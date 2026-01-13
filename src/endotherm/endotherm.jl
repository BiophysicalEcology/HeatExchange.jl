Base.@kwdef struct EndoModelPars{RE,ST,BT} <: AbstractModelParameters
    respire::RE = Param(true)
    simulsol_tolerance::ST = Param(1e-3u"K")
    resp_tolerance::BT = Param(1e-5u"K")
end

function solve_metabolic_rate(T_skin, T_insulation, o, e, m)
    e_pars = stripparams(e.environment_pars)
    e_vars = e.environment_vars
    ins = insulationpars(o)
    cond_ex = conductionpars_external(o)
    cond_in = conductionpars_internal(o)
    #conv = convectionpars(o)
    rad = radiationpars(o)
    evap = evaporationpars(o)
    hyd = hydraulicpars(o)
    resp = respirationpars(o)
    metab = metabolismpars(o)

    simulsol_out = Vector{NamedTuple}(undef, 2) # TODO preallocate
    respiration_out = Vector{NamedTuple}(undef, 1) # TODO preallocate
    geometry_pars = nothing
    simulsol_tolerance = m.simulsol_tolerance

    insulation_temperature = T_insulation * 0.7 + T_skin * 0.3
    insulation_out = insulation_properties(;
        insulation=ins, insulation_temperature, ventral_fraction=rad.ventral_fraction
    )
    fibre_diameters = insulation_out.fibre_diameters
    fibre_densities = insulation_out.fibre_densities
    insulation_depths = insulation_out.insulation_depths
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

    α_body_dorsal = rad.α_body_dorsal
    α_body_ventral = rad.α_body_ventral
    fat = Fat(cond_in.fat_fraction, cond_in.ρ_fat)

    # correct F_sky for vegetaion overhead
    F_vegetation = rad.F_sky * e_vars.shade
    F_sky = rad.F_sky - F_vegetation
    F_ground = 1 - F_sky - F_vegetation

    area_silhouette = silhouette_area(o.body, rad.solar_orientation)
    area_total = total_area(o.body)
    area_skin = skin_area(o.body)
    area_conduction = area_total * cond_ex.conduction_fraction # not used but initialised for output
    area_evaporation = evaporation_area(o.body)
    area_convection = area_total * (1 - cond_ex.conduction_fraction)
    (; Q_solar, Q_direct, Q_solar_sky, Q_solar_substrate) = solar(;
        α_body_dorsal,
        α_body_ventral,
        A_silhouette=area_silhouette,
        A_total=area_total,
        A_conduction=area_conduction,
        F_ground,
        F_sky,
        α_ground=e_pars.α_ground,
        shade=e_vars.shade,
        zenith_angle=e_vars.zenith_angle,
        global_radiation=e_vars.global_radiation,
        diffuse_fraction=e_vars.diffuse_fraction,
    )
    Q_ventral = Q_solar_substrate

    # set infrared environment
    T_vegetation = e_vars.T_air_reference # assume vegetation casting shade is at reference (e.g. 1.2m or 2m) air temperature (deg C)
    F_bush_ref = rad.F_bush # nearby bush
    F_sky_ref = F_sky # sky
    F_ground_ref = F_ground # ground
    F_vegetation_ref = F_vegetation # vegetation

    for side in 1:2

        # Calculating solar intensity entering insulation. This will depend on whether we are calculating 
        # the insulation temperature for the dorsal side or the ventral side. The dorsal side will have 
        # solar inputs from the direct beam hitting the silhouette area as well as diffuse solar scattered 
        # from the sky. The ventral side will have diffuse solar scattered off the substrate.
        # Resetting config factors and solar depending on whether the dorsal side (side=1) or 
        # ventral side (side=2) is being estimated.
        if Q_solar > 0.0u"W"
            if side == 1
                F_sky = F_sky_ref * 2.0 # proportion of upward view that is sky
                F_vegetation = F_vegetation_ref * 2.0 # proportion of upward view that is vegetation (shade)
                F_ground = 0.0
                F_bush = 0.0
                Q_sol = 2.0 * Q_direct + ((Q_solar_sky / F_sky_ref) * F_sky) # direct x 2 because assuming sun in both directions, and unadjusting Q_solar_sky for config factor imposed in SOLAR_ENDO and back to new larger one in both directions
            else
                F_sky = 0.0
                F_vegetation = 0.0
                F_ground = F_ground_ref * 2.0
                F_bush = F_bush_ref * 2.0
                Q_sol =
                    (Q_ventral / (1.0 - F_sky_ref - F_vegetation_ref)) *
                    (1.0 - (2.0 * cond_ex.conduction_fraction)) # unadjust by config factor imposed in SOLAR_ENDO to have it coming in both directions, but also cutting down according to fractional area conducting to ground (in both directions)
            end
        else
            Q_sol = 0.0u"W"
            if side == 1
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

        if side == 1
            ϵ_body = rad.ϵ_body_dorsal
        else
            ϵ_body = rad.ϵ_body_ventral
        end

        # set fur depth and conductivity
        insulation = Fur(
            insulation_depths[side + 1],
            fibre_diameters[side + 1],
            fibre_densities[side + 1],
        )
        if side == 1
            insulation_conductivity = ins.insulation_conductivity_dorsal
        else
            insulation_conductivity = ins.insulation_conductivity_ventral
        end
        geometry_pars = Body(o.body.shape, CompositeInsulation(insulation, fat))
        A_total = total_area(geometry_pars)
        r_skin = skin_radius(geometry_pars) # body radius (including fat), m
        r_insulation = r_skin + insulation.thickness # body radius including fur, m
        if geometry_pars.shape isa Cylinder || geometry_pars.shape isa Sphere
            r_compressed = r_skin + ins.insulation_depth_compressed
        else
            r_compressed = r_insulation # Note that this value is never used if conduction not being modeled, but need to have a value for the calculations
        end

        # Calculating the "cd" variable: Qcond = cd(Tskin-Tsub), where cd = Conduction area*ksub/subdepth
        if side == 2 # doing ventral side, add conduction
            area_conduction = A_total * cond_ex.conduction_fraction * 2
            cd = (area_conduction * e_vars.k_substrate) / e_pars.conduction_depth # assume conduction happens from 2.5 cm depth
        else  # doing dorsal side, no conduction. No need to adjust areas used for convection.
            area_conduction = 0.0u"m^2"
            cd = 0.0u"W/K"
        end

        # package up inputs
        geom_vars = (;
            side,
            cd,
            rad.ventral_fraction,
            conduction_fraction=cond_ex.conduction_fraction,
            longwave_depth_fraction=ins.longwave_depth_fraction,
        )
        view_factors = ViewFactors(F_sky, F_ground, F_bush, F_vegetation)
        atmos = AtmosphericConditions(e_vars.rh, e_vars.wind_speed, e_vars.P_atmos)
        env_vars = (;
            view_factors,
            atmos,
            fluid=e_pars.fluid,
            T_air=e_vars.T_air,
            T_ground=e_vars.T_ground,
            T_bush=e_vars.T_bush,
            T_vegetation,
            T_sky=e_vars.T_sky,
            T_substrate=e_vars.T_substrate,
            Q_solar=Q_sol,
            gasfrac=e_pars.gasfrac,
            convection_enhancement=e_pars.convection_enhancement,
        )
        traits = (;
            T_core=metab.T_core,
            k_flesh=cond_in.k_flesh,
            k_fat=cond_in.k_fat,
            ϵ_body,
            skin_wetness=evap.skin_wetness,
            insulation_wetness=evap.insulation_wetness,
            bare_skin_fraction=evap.bare_skin_fraction,
            eye_fraction=evap.eye_fraction,
            insulation_conductivity,
            #ψ_body = hyd.ψ_body,
        )

        insulation_out = insulation_properties(;
            insulation=ins,
            insulation_temperature=T_insulation * 0.7 + T_skin * 0.3,
            rad.ventral_fraction,
        )
        # call simulsol
        simulsol_out[side,] = simulsol(;
            T_skin,
            T_insulation,
            simulsol_tolerance=m.simulsol_tolerance,
            geometry_pars,
            insulation_pars=ins,
            insulation_out,
            geom_vars,
            env_vars,
            traits,
        )

        T_skin = simulsol_out[side,].T_skin # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
        T_insulation = simulsol_out[side,].T_insulation # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
        simulsol_tolerance = simulsol_out[side,].tolerance # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
    end

    T_skin_max = max(simulsol_out[1].T_skin, simulsol_out[2].T_skin)

    # zbrent and respfun

    # Now compute a weighted mean heat generation for all the parts/components = (dorsal value *(F_sky+F_vegetation))+(ventral value*F_ground)
    gen_d = simulsol_out[1].fluxes.Q_gen_net
    gen_v = simulsol_out[2].fluxes.Q_gen_net
    dmult = F_sky_ref + F_vegetation_ref
    vmult = 1 - dmult # assume that reflectivity of veg below equals reflectivity of soil so vmult left as 1 - dmult
    x = gen_d * dmult + gen_v * vmult # weighted estimate of metabolic heat generation
    Q_sum = x

    # reset configuration factors
    F_bush = F_bush_ref # nearby bush
    F_sky = F_sky_ref # sky
    F_ground = F_ground_ref # ground
    F_vegetation = F_vegetation_ref # vegetation

    # lung temperature and temperature of exhaled air
    T_skin = (simulsol_out[1].T_skin + simulsol_out[2].T_skin) * 0.5 # TODO weight it by dorsal/ventral fractions?
    T_insulation = (simulsol_out[1].T_insulation + simulsol_out[2].T_insulation) * 0.5
    T_lung = (metab.T_core + T_skin) * 0.5 # average of skin and core
    T_air_exit = min(e_vars.T_air + resp.Δ_breath, T_lung) # temperature of exhaled air, deg C

    if m.respire
        # now guess for metabolic rate that balances the heat budget via root-finder ZBRENT
        Q_min = metab.Q_metabolism
        Q_m1 = metab.Q_metabolism * (-2.0)
        Q_m2 = metab.Q_metabolism * 10.0
        if T_skin_max >= metab.T_core
            Q_m2 = metab.Q_metabolism * 1.01
        end
        f =
            x -> respiration(;
                Q_metab=x,
                Q_min,
                Q_sum,
                T_lung,
                mass=geometry_pars.shape.mass,
                rq=resp.rq,
                fO2_extract=resp.fO2_extract,
                rh_exit=resp.rh_exit,
                T_air_exit,
                pant=resp.pant,
                T_air=e_vars.T_air,
                P_atmos=e_vars.P_atmos,
                rh=e_vars.rh,
                gasfrac=e_pars.gasfrac,
                O2conversion=Kleiber1961(),
            ).balance

        # a, b = bracket_root(f, Q_m1, Q_m2)

        # Q_gen = find_zero(x -> respiration(; 
        #                         Q_metab = x,
        #                         Q_min,
        #                         Q_sum,
        #                         T_lung,
        #                         mass = geometry_pars.shape.mass,
        #                         rq = resp.rq,
        #                         fO2_extract = resp.fO2_extract,
        #                         rh_exit = resp.rh_exit,
        #                         T_air_exit,
        #                         pant = resp.pant,                                
        #                         T_air = e_vars.T_air,
        #                         P_atmos = e_vars.P_atmos,                                
        #                         rh = e_vars.rh,
        #                         fO2 = e_pars.fO2,
        #                         fN2 = e_pars.fN2,
        #                         fCO2 = e_pars.fCO2,                                
        #                         O2conversion = Kleiber1961()).balance,
        #                         (a, b),
        #                         Roots.Brent(),
        #                         atol = m.resp_tolerance * metab.Q_metabolism,
        #                         )
        Q_gen = zbrent(
            f,
            ustrip(u"W", Q_m1),
            ustrip(u"W", Q_m2),
            m.resp_tolerance * ustrip(u"W", metab.Q_metabolism),
        )

        respiration_out = respiration(;
            Q_metab=Q_gen,
            Q_min,
            Q_sum,
            T_lung,
            mass=geometry_pars.shape.mass,
            rq=resp.rq,
            fO2_extract=resp.fO2_extract,
            rh_exit=resp.rh_exit,
            T_air_exit,
            pant=resp.pant,
            T_air=e_vars.T_air,
            P_atmos=e_vars.P_atmos,
            rh=e_vars.rh,
            gasfrac=e_pars.gasfrac,
            O2conversion=Kleiber1961(),
        )
        Q_gen = respiration_out.Q_gen # Q_gen_net
    else
        Q_gen = Q_sum
    end

    # collate output

    # simusol outputs
    T_insulation_dorsal = simulsol_out[1].T_insulation   # feather/fur-air interface temp (°C)
    T_skin_dorsal = simulsol_out[1].T_skin   # average skin temp
    Q_convection_dorsal = simulsol_out[1].fluxes.Q_convection   # convection (W)
    Q_conduction_dorsal = simulsol_out[1].fluxes.Q_conduction   # conduction (W)
    Q_evap_skin_dorsal = simulsol_out[1].fluxes.Q_evap_skin   # cutaneous evaporation (W)
    Q_longwave_dorsal = simulsol_out[1].fluxes.Q_longwave   # radiative loss (W)
    Q_solar_dorsal = simulsol_out[1].fluxes.Q_solar   # solar (W)
    Q_evap_insulation_dorsal = simulsol_out[1].fluxes.Q_evap_insulation  # fur evaporation (W)
    ntry_dorsal = simulsol_out[1].ntry  # attempts
    success_dorsal = simulsol_out[1].success  # success?
    k_insulation_dorsal = simulsol_out[1].k_insulation  # fur conductivity? (same as FORTRAN)

    T_insulation_ventral = simulsol_out[2].T_insulation   # feather/fur-air interface temp (°C)
    T_skin_ventral = simulsol_out[2].T_skin   # average skin temp
    Q_convection_ventral = simulsol_out[2].fluxes.Q_convection   # convection (W)
    Q_conduction_ventral = simulsol_out[2].fluxes.Q_conduction   # conduction (W)
    Q_evap_skin_ventral = simulsol_out[2].fluxes.Q_evap_skin   # cutaneous evaporation (W)
    Q_longwave_ventral = simulsol_out[2].fluxes.Q_longwave   # radiative loss (W)
    Q_solar_ventral = simulsol_out[2].fluxes.Q_solar   # solar (W)
    Q_evap_insulation_ventral = simulsol_out[2].fluxes.Q_evap_insulation  # fur evaporation (W)
    ntry_ventral = simulsol_out[2].ntry  # attempts
    success_ventral = simulsol_out[2].success  # success?
    k_insulation_ventral = simulsol_out[2].k_insulation  # fur conductivity? (same as FORTRAN)
    # respiration outputs
    if m.respire
        balance = respiration_out.balance
        Q_resp = respiration_out.Q_resp
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
        Q_resp = 0.0u"W"
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
    L_v = enthalpy_of_vaporisation(e_vars.T_air)
    m_sweat = u"g/hr"((Q_evap_skin_dorsal + Q_evap_skin_ventral) * 0.5 / L_v)
    if m.respire
        m_evap = u"g/hr"(m_resp + m_sweat)
    else
        m_evap = u"g/hr"(m_sweat)
    end

    # geometric outputs
    insulation = Fur(insulation_depths[1], fibre_diameters[1], fibre_densities[1])
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
        T_surface_dorsal = T_insulation_dorsal
        T_surface_ventral = T_insulation_ventral
        area_radiant_dorsal = area_total
        area_radiant_ventral = area_total * (1 - cond_ex.conduction_fraction)
    else
        T_surface_dorsal = T_skin_dorsal
        T_surface_ventral = T_skin_ventral
        area_radiant_dorsal = area_skin
        area_radiant_ventral = area_skin #* (1 - cond_ex.conduction_fraction) TODO make conduction occur when no insulation
    end
    # infrared dorsal
    Q_rad_out_dorsal =
        2 * F_sky * σ * rad.ϵ_body_dorsal * area_radiant_dorsal * T_surface_dorsal^4
    Q_rad_in_dorsal = -Q_longwave_dorsal + Q_rad_out_dorsal

    # infrared ventral
    Q_rad_out_ventral =
        2 * F_ground * σ * rad.ϵ_body_ventral * area_radiant_ventral * T_surface_ventral^4
    Q_rad_in_ventral = -Q_longwave_ventral + Q_rad_out_ventral
    # energy flows
    Q_solar = Q_solar_dorsal * dmult + Q_solar_ventral * vmult
    Q_longwave_in = Q_rad_in_dorsal * dmult + Q_rad_in_ventral * vmult
    Q_evaporation =
        Q_evap_skin_dorsal * dmult +
        Q_evap_skin_ventral * vmult +
        Q_evap_insulation_dorsal * dmult +
        Q_evap_insulation_ventral * vmult +
        Q_resp
    Q_longwave_out = Q_rad_out_dorsal * dmult + Q_rad_out_ventral * vmult
    Q_convection = Q_convection_dorsal * dmult + Q_convection_ventral * vmult
    Q_conduction = Q_conduction_dorsal * dmult + Q_conduction_ventral * vmult

    insulation_out = insulation_properties(;
        insulation=ins,
        insulation_temperature=T_insulation * 0.7 + T_skin * 0.3,
        rad.ventral_fraction,
    )
    k_insulation_effective = insulation_out.effective_conductivities[1]
    k_insulation_compressed = insulation_out.insulation_conductivity_compressed
    T_skin = T_skin_dorsal * dmult + T_skin_ventral * vmult
    T_insulation = T_insulation_dorsal * dmult + T_insulation_ventral * vmult
    if o.body.shape isa Sphere
        shape_b = 1.0
    else
        shape_b = o.body.shape.b
    end
    thermoregulation = (;
        metab.T_core,
        T_skin,
        T_insulation,
        T_lung,
        T_skin_dorsal,
        T_skin_ventral,
        T_insulation_dorsal,
        T_insulation_ventral,
        shape_b,
        pant=resp.pant,
        skin_wetness=evap.skin_wetness,
        k_flesh=cond_in.k_flesh,
        k_insulation_effective,
        k_insulation_dorsal,
        k_insulation_ventral,
        k_insulation_compressed,
        ins.insulation_depth_dorsal,
        ins.insulation_depth_ventral,
        metab.q10,
    )

    morphology = (;
        area_total,
        area_skin,
        area_evaporation,
        area_convection,
        area_conduction=area_conduction / 2,
        area_silhouette,
        F_sky,
        F_ground,
        volume,
        volume_flesh,
        characteristic_dimension=geometry_pars.geometry.characteristic_dimension,
        fat_mass,
        geometry_pars.geometry.length...,
    )

    energy_fluxes = (;
        Q_solar,
        Q_longwave_in,
        Q_gen,
        Q_evaporation,
        Q_longwave_out,
        Q_convection,
        Q_conduction,
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
