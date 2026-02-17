"""
    SolveMetabolicRateOptions

Options controlling the endotherm metabolic rate solver.

# Fields
- `respire`: Whether to include respiration in the heat balance (default: true)
- `simulsol_tolerance`: Convergence tolerance for simultaneous temperature solution (default: 1e-3 K)
- `resp_tolerance`: Relative convergence tolerance for respiration root finding (default: 1e-5)
"""
Base.@kwdef struct SolveMetabolicRateOptions{RE,ST,BT} <: AbstractModelParameters
    respire::RE = Param(true)
    simulsol_tolerance::ST = Param(1e-3u"K")
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
    e_pars = stripparams(e.environment_pars)
    e_vars = e.environment_vars
    opts = options(o)
    ins = insulation_pars(o)
    cond_ex = conduction_pars_external(o)
    cond_in = conduction_pars_internal(o)
    #conv = convection_pars(o)
    rad = radiation_pars(o)
    evap = evaporation_pars(o)
    hyd = hydraulic_pars(o)
    resp = respiration_pars(o)
    metab = metabolism_pars(o)

    simulsol_out = Vector{NamedTuple}(undef, 2) # TODO preallocate
    respiration_out = Vector{NamedTuple}(undef, 1) # TODO preallocate
    geometry_pars = nothing
    simulsol_tolerance = opts.simulsol_tolerance

    avg_insulation_temp = insulation_temperature * 0.7 + skin_temperature * 0.3
    insulation_out = insulation_properties(ins, avg_insulation_temp, rad.ventral_fraction)
    fibres = insulation_out.fibres
    # if no insulation, reset bare_skin_fraction if necessary
    if insulation_out.insulation_test <= 0.0u"m" && evap.bare_skin_fraction < 1.0
        evap_temp = EvaporationParameters(;
            skin_wetness=evap.skin_wetness,
            insulation_wetness=evap.insulation_wetness,
            eye_fraction=evap.eye_fraction,
            bare_skin_fraction=1.0,
            insulation_fraction=evap.insulation_fraction,
        )
        evap = evap_temp
    end

    fat = Fat(cond_in.fat_fraction, cond_in.fat_density)

    # correct F_sky for vegetaion overhead
    F_vegetation = rad.sky_view_factor * e_vars.shade
    F_sky = rad.sky_view_factor - F_vegetation
    F_ground = 1 - F_sky - F_vegetation

    area_silhouette = silhouette_area(o.body, rad.solar_orientation)
    area_total = total_area(o.body)
    area_skin = skin_area(o.body)
    area_conduction = area_total * cond_ex.conduction_fraction # not used but initialised for output
    area_evaporation = evaporation_area(o.body)
    area_convection = area_total * (1 - cond_ex.conduction_fraction)
    absorptivities = Absorptivities(rad, e_pars)
    view_factors_solar = ViewFactors(F_sky, F_ground, 0.0, 0.0)
    solar_conditions = SolarConditions(e_vars)
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
    vegetation_temperature = e_vars.reference_air_temperature # assume vegetation casting shade is at reference (e.g. 1.2m or 2m) air temperature (deg C)
    F_bush_ref = rad.bush_view_factor # nearby bush
    F_sky_ref = F_sky # sky
    F_ground_ref = F_ground # ground
    F_vegetation_ref = F_vegetation # vegetation

    for (side_idx, side) in enumerate((:dorsal, :ventral))

        # Calculating solar intensity entering insulation. This will depend on whether we are calculating 
        # the insulation temperature for the dorsal side or the ventral side. The dorsal side will have 
        # solar inputs from the direct beam hitting the silhouette area as well as diffuse solar scattered 
        # from the sky. The ventral side will have diffuse solar scattered off the substrate.
        # Resetting config factors and solar depending on whether the dorsal side (side=1) or 
        # ventral side (side=2) is being estimated.
        if solar_flux > 0.0u"W"
            if side == :dorsal
                F_sky = F_sky_ref * 2.0 # proportion of upward view that is sky
                F_vegetation = F_vegetation_ref * 2.0 # proportion of upward view that is vegetation (shade)
                F_ground = 0.0
                F_bush = 0.0
                solar_flux = 2.0 * direct_flux + ((solar_sky_flux / F_sky_ref) * F_sky) # direct x 2 because assuming sun in both directions, and unadjusting solar_sky_flux for config factor imposed in SOLAR_ENDO and back to new larger one in both directions
            else
                F_sky = 0.0
                F_vegetation = 0.0
                F_ground = F_ground_ref * 2.0
                F_bush = F_bush_ref * 2.0
                solar_flux =
                    (ventral_flux / (1.0 - F_sky_ref - F_vegetation_ref)) *
                    (1.0 - (2.0 * cond_ex.conduction_fraction)) # unadjust by config factor imposed in SOLAR_ENDO to have it coming in both directions, but also cutting down according to fractional area conducting to ground (in both directions)
            end
        else
            solar_flux = 0.0u"W"
            if side == :dorsal
                F_sky = F_sky_ref * 2.0
                F_vegetation = F_vegetation_ref * 2.0
                F_ground = 0.0
                F_bush = 0.0
            else
                F_sky = 0.0
                F_vegetation = 0.0
                F_ground = F_ground_ref * 2.0
                F_bush = F_bush_ref * 2
            end
        end

        if side == :dorsal
            ϵ_body = rad.body_emissivity_dorsal
        else
            ϵ_body = rad.body_emissivity_ventral
        end

        # set fur depth and conductivity
        side_fibres = getproperty(fibres, side)
        insulation = Fur(
            side_fibres.depth,
            side_fibres.diameter,
            side_fibres.density,
        )
        geometry_pars = Body(o.body.shape, CompositeInsulation(insulation, fat))
        A_total = total_area(geometry_pars)
        r_skin = skin_radius(geometry_pars) # body radius (including fat), m
        r_insulation = r_skin + insulation.thickness # body radius including fur, m
        if geometry_pars.shape isa Cylinder || geometry_pars.shape isa Sphere
            r_compressed = r_skin + ins.depth_compressed
        else
            r_compressed = r_insulation # Note that this value is never used if conduction not being modeled, but need to have a value for the calculations
        end

        # Calculating the "cd" variable: Qcond = cd(Tskin-Tsub), where cd = Conduction area*ksub/subdepth
        if side == :ventral # doing ventral side, add conduction
            area_conduction = A_total * cond_ex.conduction_fraction * 2
            cd = (area_conduction * e_vars.k_substrate) / e_pars.conduction_depth # assume conduction happens from 2.5 cm depth
        else  # doing dorsal side, no conduction. No need to adjust areas used for convection.
            area_conduction = 0.0u"m^2"
            cd = 0.0u"W/K"
        end

        # package up inputs
        geom_vars = GeometryVariables(;
            side,
            substrate_conductance=cd,
            ventral_fraction=rad.ventral_fraction,
            conduction_fraction=cond_ex.conduction_fraction,
            longwave_depth_fraction=ins.longwave_depth_fraction,
        )
        view_factors = ViewFactors(F_sky, F_ground, F_bush, F_vegetation)
        atmos = AtmosphericConditions(e_vars)
        temperature = EnvironmentTemperatures(
            e_vars.air_temperature, e_vars.sky_temperature, e_vars.ground_temperature, vegetation_temperature, e_vars.bush_temperature, e_vars.substrate_temperature
        )
        env_vars = (;
            temperature,
            view_factors,
            atmos,
            fluid=e_pars.fluid,
            solar_flux=solar_flux,
            gas_fractions=e_pars.gas_fractions,
            convection_enhancement=e_pars.convection_enhancement,
        )
        traits = (;
            core_temperature=metab.core_temperature,
            k_flesh=cond_in.flesh_conductivity,
            k_fat=cond_in.fat_conductivity,
            ϵ_body,
            skin_wetness=evap.skin_wetness,
            insulation_wetness=evap.insulation_wetness,
            bare_skin_fraction=evap.bare_skin_fraction,
            eye_fraction=evap.eye_fraction,
        )

        insulation_out = insulation_properties(
            ins, insulation_temperature * 0.7 + skin_temperature * 0.3, rad.ventral_fraction
        )
        # call simulsol
        simulsol_out[side_idx] = simulsol(;
            body=geometry_pars,
            insulation_pars=ins,
            insulation_out,
            geom_vars,
            env_vars,
            traits,
            simulsol_tolerance=opts.simulsol_tolerance,
            skin_temperature,
            insulation_temperature,
        )

        skin_temperature = simulsol_out[side_idx].skin_temperature # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
        insulation_temperature = simulsol_out[side_idx].insulation_temperature # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
        simulsol_tolerance = simulsol_out[side_idx].tolerance # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
    end

    max_skin_temperature = max(simulsol_out[1].skin_temperature, simulsol_out[2].skin_temperature)

    # zbrent and respfun

    # Now compute a weighted mean heat generation for all the parts/components = (dorsal value *(F_sky+F_vegetation))+(ventral value*F_ground)
    gen_d = simulsol_out[1].fluxes.net_generated_flux
    gen_v = simulsol_out[2].fluxes.net_generated_flux
    dmult = F_sky_ref + F_vegetation_ref
    vmult = 1 - dmult # assume that reflectivity of veg below equals reflectivity of soil so vmult left as 1 - dmult
    x = gen_d * dmult + gen_v * vmult # weighted estimate of metabolic heat generation
    flux_sum = x

    # reset configuration factors
    F_bush = F_bush_ref # nearby bush
    F_sky = F_sky_ref # sky
    F_ground = F_ground_ref # ground
    F_vegetation = F_vegetation_ref # vegetation

    # lung temperature and temperature of exhaled air
    skin_temperature = (simulsol_out[1].skin_temperature + simulsol_out[2].skin_temperature) * 0.5 # TODO weight it by dorsal/ventral fractions?
    insulation_temperature = (simulsol_out[1].insulation_temperature + simulsol_out[2].insulation_temperature) * 0.5
    lung_temperature = (metab.core_temperature + skin_temperature) * 0.5 # average of skin and core
    exit_air_temperature = min(e_vars.air_temperature + resp.exhaled_temperature_offset, lung_temperature) # temperature of exhaled air, deg C

    if opts.respire
        # now guess for metabolic rate that balances the heat budget via root-finder ZBRENT
        minimum_flux = metab.metabolic_fluxolism
        flux_m1 = metab.metabolic_fluxolism * (-2.0)
        flux_m2 = metab.metabolic_fluxolism * 10.0
        if max_skin_temperature >= metab.core_temperature
            flux_m2 = metab.metabolic_fluxolism * 1.01
        end
        resp_atmos = AtmosphericConditions(e_vars)
        f =
            x -> respiration(
                MetabolicRates(; metabolic=x, sum=flux_sum, minimum=minimum_flux),
                resp,
                resp_atmos,
                geometry_pars.shape.mass,
                lung_temperature,
                e_vars.air_temperature;
                gas_fractions=e_pars.gas_fractions,
                O2conversion=Kleiber1961(),
            ).balance

        generated_flux = zbrent(
            f,
            ustrip(u"W", flux_m1),
            ustrip(u"W", flux_m2),
            opts.resp_tolerance * ustrip(u"W", metab.metabolic_fluxolism),
        )

        respiration_out = respiration(
            MetabolicRates(; metabolic=generated_flux, sum=flux_sum, minimum=minimum_flux),
            resp,
            resp_atmos,
            geometry_pars.shape.mass,
            lung_temperature,
            e_vars.air_temperature;
            gas_fractions=e_pars.gas_fractions,
            O2conversion=Kleiber1961(),
        )
        generated_flux = respiration_out.generated_flux # net_generated_flux
    else
        generated_flux = flux_sum
    end

    # collate output

    # simusol outputs
    insulation_temperature_dorsal = simulsol_out[1].insulation_temperature   # feather/fur-air interface temp (°C)
    skin_temperature_dorsal = simulsol_out[1].skin_temperature   # average skin temp
    dorsal_convection_flux = simulsol_out[1].fluxes.convection_flux   # convection (W)
    dorsal_conduction_flux = simulsol_out[1].fluxes.conduction_flux   # conduction (W)
    dorsal_skin_evaporation_flux = simulsol_out[1].fluxes.skin_evaporation_flux   # cutaneous evaporation (W)
    dorsal_longwave_flux = simulsol_out[1].fluxes.longwave_flux   # radiative loss (W)
    dorsal_solar_flux = simulsol_out[1].fluxes.solar_flux   # solar (W)
    dorsal_insulation_evaporation_flux = simulsol_out[1].fluxes.insulation_evaporation_flux  # fur evaporation (W)
    ntry_dorsal = simulsol_out[1].ntry  # attempts
    success_dorsal = simulsol_out[1].success  # success?
    k_insulation_dorsal = simulsol_out[1].k_insulation  # fur conductivity? (same as FORTRAN)

    insulation_temperature_ventral = simulsol_out[2].insulation_temperature   # feather/fur-air interface temp (°C)
    skin_temperature_ventral = simulsol_out[2].skin_temperature   # average skin temp
    ventral_convection_flux = simulsol_out[2].fluxes.convection_flux   # convection (W)
    ventral_conduction_flux = simulsol_out[2].fluxes.conduction_flux   # conduction (W)
    ventral_skin_evaporation_flux = simulsol_out[2].fluxes.skin_evaporation_flux   # cutaneous evaporation (W)
    ventral_longwave_flux = simulsol_out[2].fluxes.longwave_flux   # radiative loss (W)
    ventral_solar_flux = simulsol_out[2].fluxes.solar_flux   # solar (W)
    ventral_insulation_evaporation_flux = simulsol_out[2].fluxes.insulation_evaporation_flux  # fur evaporation (W)
    ntry_ventral = simulsol_out[2].ntry  # attempts
    success_ventral = simulsol_out[2].success  # success?
    k_insulation_ventral = simulsol_out[2].k_insulation  # fur conductivity? (same as FORTRAN)
    # respiration outputs
    if opts.respire
        balance = respiration_out.balance
        respiration_flux = respiration_out.respiration_flux
        m_resp = respiration_out.m_resp
        V_air = respiration_out.V_air
        V_O2_STP = respiration_out.V_O2_STP
        molar_fluxes = respiration_out.molar_fluxes
        J_air_in = molar_fluxes.J_air_in
        J_air_out = molar_fluxes.J_air_out
        J_H2O_in = molar_fluxes.J_H2O_in
        J_H2O_out = molar_fluxes.J_H2O_out
        J_O2_in = molar_fluxes.J_O2_in
        J_O2_out = molar_fluxes.J_O2_out
        J_CO2_in = molar_fluxes.J_CO2_in
        J_CO2_out = molar_fluxes.J_CO2_out
        J_N2_in = molar_fluxes.J_N2_in
        J_N2_out = molar_fluxes.J_N2_out
    else
        balance = nothing
        respiration_flux = 0.0u"W"
        m_resp = nothing
        V_air = nothing
        V_O2_STP = nothing
        J_air_in = nothing
        J_air_out = nothing
        J_H2O_in = nothing
        J_H2O_out = nothing
        J_O2_in = nothing
        J_O2_out = nothing
        J_CO2_in = nothing
        J_CO2_out = nothing
        J_N2_in = nothing
        J_N2_out = nothing
    end
    # evaporation calculations
    L_v = enthalpy_of_vaporisation(e_vars.air_temperature)
    m_sweat = u"g/hr"((dorsal_skin_evaporation_flux + ventral_skin_evaporation_flux) * 0.5 / L_v)
    if opts.respire
        m_evap = u"g/hr"(m_resp + m_sweat)
    else
        m_evap = u"g/hr"(m_sweat)
    end

    # geometric outputs
    insulation = Fur(fibres.average.depth, fibres.average.diameter, fibres.average.density)
    geometry_pars = Body(o.body.shape, CompositeInsulation(insulation, fat))

    fat_mass = geometry_pars.shape.mass * fat.fraction
    volume = geometry_pars.geometry.volume
    volume_flesh = flesh_volume(geometry_pars)
    area_total = total_area(geometry_pars)
    area_skin = skin_area(geometry_pars)
    area_silhouette = silhouette_area(geometry_pars, rad.solar_orientation)

    # radiation outputs

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    if insulation_out.insulation_test > 0.0u"m"
        surface_temperature_dorsal = insulation_temperature_dorsal
        surface_temperature_ventral = insulation_temperature_ventral
        area_radiant_dorsal = area_total
        area_radiant_ventral = area_total * (1 - cond_ex.conduction_fraction)
    else
        surface_temperature_dorsal = skin_temperature_dorsal
        surface_temperature_ventral = skin_temperature_ventral
        area_radiant_dorsal = area_skin
        area_radiant_ventral = area_skin #* (1 - cond_ex.conduction_fraction) TODO make conduction occur when no insulation
    end
    # infrared dorsal
    dorsal_radiation_out_flux =
        2 * F_sky * σ * rad.body_emissivity_dorsal * area_radiant_dorsal * surface_temperature_dorsal^4
    dorsal_radiation_in_flux = -dorsal_longwave_flux + dorsal_radiation_out_flux

    # infrared ventral
    ventral_radiation_out_flux =
        2 * F_ground * σ * rad.body_emissivity_ventral * area_radiant_ventral * surface_temperature_ventral^4
    ventral_radiation_in_flux = -ventral_longwave_flux + ventral_radiation_out_flux
    # energy flows
    solar_flux = dorsal_solar_flux * dmult + ventral_solar_flux * vmult
    longwave_in_flux = dorsal_radiation_in_flux * dmult + ventral_radiation_in_flux * vmult
    evaporation_fluxoration =
        dorsal_skin_evaporation_flux * dmult +
        ventral_skin_evaporation_flux * vmult +
        dorsal_insulation_evaporation_flux * dmult +
        ventral_insulation_evaporation_flux * vmult +
        respiration_flux
    longwave_out_flux = dorsal_radiation_out_flux * dmult + ventral_radiation_out_flux * vmult
    convection_flux = dorsal_convection_flux * dmult + ventral_convection_flux * vmult
    conduction_flux = dorsal_conduction_flux * dmult + ventral_conduction_flux * vmult

    insulation_out = insulation_properties(
        ins, insulation_temperature * 0.7 + skin_temperature * 0.3, rad.ventral_fraction
    )
    k_insulation_effective = insulation_out.conductivities.average
    k_insulation_compressed = insulation_out.conductivity_compressed
    skin_temperature = skin_temperature_dorsal * dmult + skin_temperature_ventral * vmult
    insulation_temperature = insulation_temperature_dorsal * dmult + insulation_temperature_ventral * vmult
    if o.body.shape isa Sphere
        shape_b = 1.0
    else
        shape_b = o.body.shape.b
    end
    thermoregulation = (;
        metab.core_temperature,
        skin_temperature,
        insulation_temperature,
        lung_temperature,
        skin_temperature_dorsal,
        skin_temperature_ventral,
        insulation_temperature_dorsal,
        insulation_temperature_ventral,
        shape_b,
        pant=resp.pant,
        skin_wetness=evap.skin_wetness,
        k_flesh=cond_in.flesh_conductivity,
        k_insulation_effective,
        k_insulation_dorsal,
        k_insulation_ventral,
        k_insulation_compressed,
        insulation_depth_dorsal=ins.dorsal.depth,
        insulation_depth_ventral=ins.ventral.depth,
        metab.q10,
    )

    morphology = (;
        area_total,
        area_skin,
        area_evaporation,
        area_convection,
        area_conduction=area_conduction / 2,
        area_silhouette,
        sky_view_factor=F_sky,
        ground_view_factor=F_ground,
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
        evaporation_fluxoration,
        longwave_out_flux,
        convection_flux,
        conduction_flux,
        balance,
        ntry=max(ntry_dorsal, ntry_ventral),
        success=all((success_dorsal, success_ventral)),
    )

    mass_fluxes = (;
        V_air,
        V_O2_STP,
        m_evap,
        m_resp,
        m_sweat,
        J_H2O_in,
        J_H2O_out,
        J_O2_in,
        J_O2_out,
        J_CO2_in,
        J_CO2_out,
        J_N2_in,
        J_N2_out,
        J_air_in,
        J_air_out,
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
