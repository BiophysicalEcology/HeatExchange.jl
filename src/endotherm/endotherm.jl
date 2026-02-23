"""
    SolveMetabolicRateOptions

Options controlling the endotherm metabolic rate solver.

# Fields
- `respire`: Whether to include respiration in the heat balance (default: true)
- `temperature_tolerance`: Convergence tolerance for simultaneous temperature solution (default: 1e-3 K)
- `resp_tolerance`: Relative convergence tolerance for respiration root finding (default: 1e-5)
"""
Base.@kwdef struct SolveMetabolicRateOptions{RE,ST,BT} <: AbstractModelParameters
    respire::RE = Param(true)
    temperature_tolerance::ST = Param(1e-3u"K")
    resp_tolerance::BT = Param(1e-5)
end

"""
    solve_metabolic_rate(o::Organism, e, skin_temperature, insulation_temperature)

Solve for the metabolic rate that balances heat production with heat loss for an endotherm.

This is the main entry point for endotherm heat balance calculations. It computes heat
fluxes for dorsal and ventral body surfaces separately, then uses root-finding to
determine the metabolic rate that satisfies the respiratory heat balance.

# Arguments
- `o::Organism`: Organism with body geometry and traits
- `e`: Environment containing `environment_pars` and `environment_vars`
- `skin_temperature`: Initial guess for skin temperature
- `insulation_temperature`: Initial guess for insulation surface temperature

# Returns
NamedTuple with:
- `thermoregulation`: Temperature and conductivity outputs
- `morphology`: Body geometry and areas
- `energy_fluxes`: Heat flux components (solar, longwave, convection, etc.)
- `mass_fluxes`: Water and gas exchange rates
"""
function solve_metabolic_rate(o::Organism, e, skin_temperature, insulation_temperature)
    environment_pars = stripparams(e.environment_pars)
    environment_vars = e.environment_vars
    opts = options(o)
    ins_pars = insulation_pars(o)
    external_conduction = conduction_pars_external(o)
    internal_conduction = conduction_pars_internal(o)
    #conv = convection_pars(o)
    rad_pars = radiation_pars(o)
    evap_pars = evaporation_pars(o)
    hyd_pars = hydraulic_pars(o)
    resp_pars = respiration_pars(o)
    metab_pars = metabolism_pars(o)

    temps_out = Vector{NamedTuple}(undef, 2) # TODO preallocate
    respiration_out = Vector{NamedTuple}(undef, 1) # TODO preallocate
    geometry_pars = nothing
    temperature_tolerance = opts.temperature_tolerance

    avg_insulation_temp = insulation_temperature * 0.7 + skin_temperature * 0.3
    insulation = insulation_properties(ins_pars, avg_insulation_temp, rad_pars.ventral_fraction)
    fibres = insulation.fibres
    # if no insulation, reset bare_skin_fraction if necessary
    if insulation.insulation_test <= 0.0u"m" && evap_pars.bare_skin_fraction < 1.0
        evap_temp = EvaporationParameters(;
            skin_wetness=evap_pars.skin_wetness,
            insulation_wetness=evap_pars.insulation_wetness,
            eye_fraction=evap_pars.eye_fraction,
            bare_skin_fraction=1.0,
            insulation_fraction=evap_pars.insulation_fraction,
        )
        evap_pars = evap_temp
    end

    fat = Fat(internal_conduction.fat_fraction, internal_conduction.fat_density)

    # correct sky_factor for vegetation overhead
    vegetation_factor = rad_pars.sky_view_factor * environment_vars.shade
    sky_factor = rad_pars.sky_view_factor - vegetation_factor
    ground_factor = 1 - sky_factor - vegetation_factor

    area_silhouette = silhouette_area(o.body, rad_pars.solar_orientation)
    total_area = BiophysicalGeometry.total_area(o.body)
    area_skin = skin_area(o.body)
    area_conduction = total_area * external_conduction.conduction_fraction # not used but initialised for output
    area_evaporation = evaporation_area(o.body)
    area_convection = total_area * (1 - external_conduction.conduction_fraction)
    absorptivities = Absorptivities(rad_pars, environment_pars)
    view_factors_solar = ViewFactors(sky_factor, ground_factor, 0.0, 0.0)
    solar_conditions = SolarConditions(environment_vars)
    (; solar_flux, direct_flux, solar_sky_flux, solar_substrate_flux) = solar(
        o.body,
        absorptivities,
        view_factors_solar,
        solar_conditions,
        area_silhouette,
        area_conduction,
    )
    ventral_flux = solar_substrate_flux

    # set infrared environment
    vegetation_temperature = environment_vars.reference_air_temperature # assume vegetation casting shade is at reference (e.g. 1.2m or 2m) air temperature (deg C)
    bush_factor_ref = rad_pars.bush_view_factor # nearby bush
    sky_factor_ref = sky_factor # sky
    ground_factor_ref = ground_factor # ground
    vegetation_factor_ref = vegetation_factor # vegetation

    for (side_idx, side) in enumerate((:dorsal, :ventral))

        # Calculating solar intensity entering insulation. This will depend on whether we are calculating
        # the insulation temperature for the dorsal side or the ventral side. The dorsal side will have
        # solar inputs from the direct beam hitting the silhouette area as well as diffuse solar scattered
        # from the sky. The ventral side will have diffuse solar scattered off the substrate.
        # Resetting config factors and solar depending on whether the dorsal side (side=1) or
        # ventral side (side=2) is being estimated.
        side_sky_factor, side_ground_factor, side_bush_factor, side_vegetation_factor = if solar_flux > 0.0u"W"
            if side == :dorsal
                solar_flux = 2.0 * direct_flux + ((solar_sky_flux / sky_factor_ref) * sky_factor_ref * 2.0) # direct x 2 because assuming sun in both directions
                (sky_factor_ref * 2.0, 0.0, 0.0, vegetation_factor_ref * 2.0)
            else
                solar_flux =
                    (ventral_flux / (1.0 - sky_factor_ref - vegetation_factor_ref)) *
                    (1.0 - (2.0 * external_conduction.conduction_fraction)) # unadjust by config factor imposed in SOLAR_ENDO
                (0.0, ground_factor_ref * 2.0, bush_factor_ref * 2.0, 0.0)
            end
        else
            solar_flux = 0.0u"W"
            if side == :dorsal
                (sky_factor_ref * 2.0, 0.0, 0.0, vegetation_factor_ref * 2.0)
            else
                (0.0, ground_factor_ref * 2.0, bush_factor_ref * 2.0, 0.0)
            end
        end

        if side == :dorsal
            ϵ_body = rad_pars.body_emissivity_dorsal
        else
            ϵ_body = rad_pars.body_emissivity_ventral
        end

        # set fur depth and conductivity
        side_fibres = getproperty(fibres, side)
        fur = Fur(
            side_fibres.depth,
            side_fibres.diameter,
            side_fibres.density,
        )
        geometry_pars = Body(o.body.shape, CompositeInsulation(fur, fat))
        body_total_area = BiophysicalGeometry.total_area(geometry_pars)
        r_skin = skin_radius(geometry_pars) # body radius (including fat), m
        r_insulation = r_skin + fur.thickness # body radius including fur, m
        if geometry_pars.shape isa Cylinder || geometry_pars.shape isa Sphere
            r_compressed = r_skin + ins_pars.depth_compressed
        else
            r_compressed = r_insulation # Note that this value is never used if conduction not being modeled, but need to have a value for the calculations
        end

        # Calculating substrate_conductance: Qcond = substrate_conductance*(Tskin-Tsub)
        substrate_conductance = if side == :ventral # doing ventral side, add conduction
            area_conduction = body_total_area * external_conduction.conduction_fraction * 2
            (area_conduction * environment_vars.substrate_conductivity) / environment_pars.conduction_depth # assume conduction happens from 2.5 cm depth
        else  # doing dorsal side, no conduction. No need to adjust areas used for convection.
            0.0u"W/K"
        end

        # package up inputs
        geometry_vars = GeometryVariables(;
            side,
            substrate_conductance,
            ventral_fraction=rad_pars.ventral_fraction,
            conduction_fraction=external_conduction.conduction_fraction,
            longwave_depth_fraction=ins_pars.longwave_depth_fraction,
        )
        view_factors = ViewFactors(side_sky_factor, side_ground_factor, side_bush_factor, side_vegetation_factor)
        atmos = AtmosphericConditions(environment_vars)
        temperature = EnvironmentTemperatures(
            environment_vars.air_temperature, environment_vars.sky_temperature, environment_vars.ground_temperature, vegetation_temperature, environment_vars.bush_temperature, environment_vars.substrate_temperature
        )
        env_packed = (;
            temperature,
            view_factors,
            atmos,
            fluid=environment_pars.fluid,
            solar_flux,
            gas_fractions=environment_pars.gas_fractions,
            convection_enhancement=environment_pars.convection_enhancement,
        )
        traits = (;
            core_temperature=metab_pars.core_temperature,
            flesh_conductivity=internal_conduction.flesh_conductivity,
            fat_conductivity=internal_conduction.fat_conductivity,
            ϵ_body,
            skin_wetness=evap_pars.skin_wetness,
            insulation_wetness=evap_pars.insulation_wetness,
            bare_skin_fraction=evap_pars.bare_skin_fraction,
            eye_fraction=evap_pars.eye_fraction,
        )

        insulation = insulation_properties(
            ins_pars, insulation_temperature * 0.7 + skin_temperature * 0.3, rad_pars.ventral_fraction
        )
        # call solve_temperatures
        temps_out[side_idx] = solve_temperatures(;
            body=geometry_pars,
            insulation_pars=ins_pars,
            insulation,
            geometry_vars,
            environment_vars=env_packed,
            traits,
            temperature_tolerance=opts.temperature_tolerance,
            skin_temperature,
            insulation_temperature,
        )

        skin_temperature = temps_out[side_idx].skin_temperature # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
        insulation_temperature = temps_out[side_idx].insulation_temperature # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
        temperature_tolerance = temps_out[side_idx].tolerance # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
    end

    max_skin_temperature = max(temps_out[1].skin_temperature, temps_out[2].skin_temperature)

    # zbrent and respfun

    # Now compute a weighted mean heat generation for all the parts/components = (dorsal value *(sky_factor+vegetation_factor))+(ventral value*ground_factor)
    gen_d = temps_out[1].fluxes.net_generated
    gen_v = temps_out[2].fluxes.net_generated
    dmult = sky_factor_ref + vegetation_factor_ref
    vmult = 1 - dmult # assume that reflectivity of veg below equals reflectivity of soil so vmult left as 1 - dmult
    x = gen_d * dmult + gen_v * vmult # weighted estimate of metabolic heat generation
    flux_sum = x

    # lung temperature and temperature of exhaled air
    skin_temperature = (temps_out[1].skin_temperature + temps_out[2].skin_temperature) * 0.5 # TODO weight it by dorsal/ventral fractions?
    insulation_temperature = (temps_out[1].insulation_temperature + temps_out[2].insulation_temperature) * 0.5
    lung_temperature = (metab_pars.core_temperature + skin_temperature) * 0.5 # average of skin and core
    exit_air_temperature = min(environment_vars.air_temperature + resp_pars.exhaled_temperature_offset, lung_temperature) # temperature of exhaled air, deg C

    if opts.respire
        # now guess for metabolic rate that balances the heat budget via root-finder ZBRENT
        minimum_flux = metab_pars.metabolic_flux
        flux_m1 = metab_pars.metabolic_flux * (-2.0)
        flux_m2 = metab_pars.metabolic_flux * 10.0
        if max_skin_temperature >= metab_pars.core_temperature
            flux_m2 = metab_pars.metabolic_flux * 1.01
        end
        resp_atmos = AtmosphericConditions(environment_vars)
        f =
            x -> respiration(
                MetabolicRates(; metabolic=x, sum=flux_sum, minimum=minimum_flux),
                resp_pars,
                resp_atmos,
                geometry_pars.shape.mass,
                lung_temperature,
                environment_vars.air_temperature;
                gas_fractions=environment_pars.gas_fractions,
                O2conversion=Kleiber1961(),
            ).balance

        generated_flux = zbrent(
            f,
            ustrip(u"W", flux_m1),
            ustrip(u"W", flux_m2),
            opts.resp_tolerance * ustrip(u"W", metab_pars.metabolic_flux),
        )

        respiration_out = respiration(
            MetabolicRates(; metabolic=generated_flux, sum=flux_sum, minimum=minimum_flux),
            resp_pars,
            resp_atmos,
            geometry_pars.shape.mass,
            lung_temperature,
            environment_vars.air_temperature;
            gas_fractions=environment_pars.gas_fractions,
            O2conversion=Kleiber1961(),
        )
        generated_flux = respiration_out.generated_flux # net_generated_flux
    else
        generated_flux = flux_sum
    end

    # collate output

    # solve_temperatures outputs
    insulation_temperature_dorsal = temps_out[1].insulation_temperature   # feather/fur-air interface temp (°C)
    skin_temperature_dorsal = temps_out[1].skin_temperature   # average skin temp
    dorsal_convection_flux = temps_out[1].fluxes.convection   # convection (W)
    dorsal_conduction_flux = temps_out[1].fluxes.conduction   # conduction (W)
    dorsal_skin_evaporation_flux = temps_out[1].fluxes.skin_evaporation   # cutaneous evaporation (W)
    dorsal_longwave_flux = temps_out[1].fluxes.longwave   # radiative loss (W)
    dorsal_solar_flux = temps_out[1].fluxes.solar   # solar (W)
    dorsal_insulation_evaporation_flux = temps_out[1].fluxes.insulation_evaporation  # fur evaporation (W)
    ntry_dorsal = temps_out[1].ntry  # attempts
    success_dorsal = temps_out[1].success  # success?
    insulation_conductivity_dorsal = temps_out[1].insulation_conductivity

    insulation_temperature_ventral = temps_out[2].insulation_temperature   # feather/fur-air interface temp (°C)
    skin_temperature_ventral = temps_out[2].skin_temperature   # average skin temp
    ventral_convection_flux = temps_out[2].fluxes.convection   # convection (W)
    ventral_conduction_flux = temps_out[2].fluxes.conduction   # conduction (W)
    ventral_skin_evaporation_flux = temps_out[2].fluxes.skin_evaporation   # cutaneous evaporation (W)
    ventral_longwave_flux = temps_out[2].fluxes.longwave   # radiative loss (W)
    ventral_solar_flux = temps_out[2].fluxes.solar   # solar (W)
    ventral_insulation_evaporation_flux = temps_out[2].fluxes.insulation_evaporation  # fur evaporation (W)
    ntry_ventral = temps_out[2].ntry  # attempts
    success_ventral = temps_out[2].success  # success?
    insulation_conductivity_ventral = temps_out[2].insulation_conductivity
    # respiration outputs
    if opts.respire
        balance = respiration_out.balance
        respiration_flux = respiration_out.respiration_flux
        respiration_mass = respiration_out.respiration_mass
        air_flow = respiration_out.air_flow
        oxygen_flow_stp = respiration_out.oxygen_flow_stp
        molar_fluxes_in = respiration_out.molar_fluxes_in
        molar_fluxes_out = respiration_out.molar_fluxes_out
    else
        balance = nothing
        respiration_flux = 0.0u"W"
        respiration_mass = nothing
        air_flow = nothing
        oxygen_flow_stp = nothing
        molar_fluxes_in = nothing
        molar_fluxes_out = nothing
    end
    # evaporation calculations
    latent_heat_vaporisation = enthalpy_of_vaporisation(environment_vars.air_temperature)
    m_sweat = u"g/hr"((dorsal_skin_evaporation_flux + ventral_skin_evaporation_flux) * 0.5 / latent_heat_vaporisation)
    if opts.respire
        m_evap = u"g/hr"(respiration_mass + m_sweat)
    else
        m_evap = u"g/hr"(m_sweat)
    end

    # geometric outputs
    fur = Fur(fibres.average.depth, fibres.average.diameter, fibres.average.density)
    geometry_pars = Body(o.body.shape, CompositeInsulation(fur, fat))

    fat_mass = geometry_pars.shape.mass * fat.fraction
    volume = geometry_pars.geometry.volume
    volume_flesh = flesh_volume(geometry_pars)
    total_area = BiophysicalGeometry.total_area(geometry_pars)
    area_skin = skin_area(geometry_pars)
    area_silhouette = silhouette_area(geometry_pars, rad_pars.solar_orientation)

    # radiation outputs

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    if insulation.insulation_test > 0.0u"m"
        surface_temperature_dorsal = insulation_temperature_dorsal
        surface_temperature_ventral = insulation_temperature_ventral
        area_radiant_dorsal = total_area
        area_radiant_ventral = total_area * (1 - external_conduction.conduction_fraction)
    else
        surface_temperature_dorsal = skin_temperature_dorsal
        surface_temperature_ventral = skin_temperature_ventral
        area_radiant_dorsal = area_skin
        area_radiant_ventral = area_skin #* (1 - external_conduction.conduction_fraction) TODO make conduction occur when no insulation
    end
    # infrared dorsal
    dorsal_radiation_out_flux =
        2 * sky_factor_ref * σ * rad_pars.body_emissivity_dorsal * area_radiant_dorsal * surface_temperature_dorsal^4
    dorsal_radiation_in_flux = -dorsal_longwave_flux + dorsal_radiation_out_flux

    # infrared ventral
    ventral_radiation_out_flux =
        2 * ground_factor_ref * σ * rad_pars.body_emissivity_ventral * area_radiant_ventral * surface_temperature_ventral^4
    ventral_radiation_in_flux = -ventral_longwave_flux + ventral_radiation_out_flux
    # energy flows
    solar_flux = dorsal_solar_flux * dmult + ventral_solar_flux * vmult
    longwave_in_flux = dorsal_radiation_in_flux * dmult + ventral_radiation_in_flux * vmult
    evaporation_flux =
        dorsal_skin_evaporation_flux * dmult +
        ventral_skin_evaporation_flux * vmult +
        dorsal_insulation_evaporation_flux * dmult +
        ventral_insulation_evaporation_flux * vmult +
        respiration_flux
    longwave_out_flux = dorsal_radiation_out_flux * dmult + ventral_radiation_out_flux * vmult
    convection_flux = dorsal_convection_flux * dmult + ventral_convection_flux * vmult
    conduction_flux = dorsal_conduction_flux * dmult + ventral_conduction_flux * vmult

    insulation = insulation_properties(
        ins_pars, insulation_temperature * 0.7 + skin_temperature * 0.3, rad_pars.ventral_fraction
    )
    insulation_conductivity_effective = insulation.conductivities.average
    insulation_conductivity_compressed = insulation.conductivity_compressed
    skin_temperature = skin_temperature_dorsal * dmult + skin_temperature_ventral * vmult
    insulation_temperature = insulation_temperature_dorsal * dmult + insulation_temperature_ventral * vmult
    if o.body.shape isa Sphere
        shape_b = 1.0
    else
        shape_b = o.body.shape.b
    end
    thermoregulation = (;
        metab_pars.core_temperature,
        skin_temperature,
        insulation_temperature,
        lung_temperature,
        skin_temperature_dorsal,
        skin_temperature_ventral,
        insulation_temperature_dorsal,
        insulation_temperature_ventral,
        shape_b,
        pant=resp_pars.pant,
        skin_wetness=evap_pars.skin_wetness,
        flesh_conductivity=internal_conduction.flesh_conductivity,
        insulation_conductivity_effective,
        insulation_conductivity_dorsal,
        insulation_conductivity_ventral,
        insulation_conductivity_compressed,
        insulation_depth_dorsal=ins_pars.dorsal.depth,
        insulation_depth_ventral=ins_pars.ventral.depth,
        metab_pars.q10,
    )

    morphology = (;
        total_area,
        area_skin,
        area_evaporation,
        area_convection,
        area_conduction=area_conduction / 2,
        area_silhouette,
        sky_view_factor=sky_factor_ref,
        ground_view_factor=ground_factor_ref,
        volume,
        volume_flesh,
        characteristic_dimension=geometry_pars.geometry.characteristic_dimension,
        fat_mass,
        geometry_pars.geometry.length...,
    )

    energy_fluxes = (;
        solar_flux,
        longwave_in_flux,
        generated_flux,
        evaporation_flux,
        longwave_out_flux,
        convection_flux,
        conduction_flux,
        balance,
        ntry=max(ntry_dorsal, ntry_ventral),
        success=all((success_dorsal, success_ventral)),
    )

    mass_fluxes = (;
        air_flow,
        oxygen_flow_stp,
        m_evap,
        respiration_mass,
        m_sweat,
        molar_fluxes_in,
        molar_fluxes_out,
    )
    return (; thermoregulation, morphology, energy_fluxes, mass_fluxes)
end

function bracket_root(f, a, b; factor=2, maxiter=20)
    fa = f(a)
    fb = f(b)
    for _ in 1:maxiter
        fa * fb < 0u"W^2" && return (a, b)
        a -= factor * (b - a)
        b += factor * (b - a)
        fa = f(a)
        fb = f(b)
    end
    error("Failed to bracket root")
end

function zbrent(f, a, b, tol; maxiter=300, eps=3e-8)
    qa = ustrip(u"W", f(a * u"W"))
    qb = ustrip(u"W", f(b * u"W"))

    c = 0.0
    e = 0.0
    d = 0.0

    qc = qb

    @inbounds for i in 1:maxiter
        if qb * qc > 0.0
            c = a
            qc = qa
            d = b - a
            e = d
        end
        if abs(qc) < abs(qb)
            a = b
            b = c
            c = a
            qa = qb
            qb = qc
            qc = qa
        end

        tol1 = 2eps * abs(b) + tol/2
        xm = (c - b) / 2

        if abs(xm) ≤ tol1 || qb == 0.0
            return b * u"W"
        end
        if abs(qb) ≤ tol1 && i > 1
            return b * u"W"
        end
        if abs(e) ≥ tol1 && abs(qa) > abs(qb)
            s = qb / qa
            if a == c
                p = 2xm * s
                q = 1 - s
            else
                q = qa / qc
                r = qb / qc
                p = s * (2xm*q*(q - r) - (b - a)*(r - 1))
                q = (q - 1)*(r - 1)*(s - 1)
            end

            if p > 0
                q = -q
            end
            p = abs(p)

            if 2p < min(3xm*q - abs(tol1*q), abs(e*q))
                e = d
                d = p / q
            else
                d = xm
                e = d
            end
        else
            d = xm
            e = d
        end

        a = b
        qa = qb

        b += abs(d) > tol1 ? d : sign(xm)*tol1
        qb = ustrip(u"W", f(b * u"W"))
    end

    return b * u"W"
end
