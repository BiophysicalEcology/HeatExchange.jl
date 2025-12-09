Base.@kwdef struct EndoModelPars{RE,ST,BT} <: AbstractModelParameters
    respire::RE =               Param(true)
    simulsol_tolerance::ST =    Param(1e-3u"K")
    resp_tolerance::BT =        Param(1e-5u"K")
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
                        insulation=ins, insulation_temperature, 
                        ventral_fraction = rad.ventral_fraction)
    fibre_diameters = insulation_out.fibre_diameters
    fibre_densities = insulation_out.fibre_densities
    insulation_depths = insulation_out.insulation_depths

    # if no insulation, reset bare_skin_fraction if necessary
    if insulation_out.insulation_test <= 0.0u"m" && evap.bare_skin_fraction < 1.0
        evap_temp = EvaporationParameters(
            skin_wetness        = evap.skin_wetness,
            insulation_wetness  = evap.insulation_wetness,
            eye_fraction        = evap.eye_fraction,
            bare_skin_fraction  = 1.0,
            insulation_fraction = evap.insulation_fraction,
        )    
        evap = evap_temp
    end

    α_body_dorsal = rad.α_body_dorsal
    α_body_ventral = rad.α_body_ventral
    fat = Fat(cond_in.fat_fraction, cond_in.ρ_fat)

    # correct F_sky for vegetaion overhead
    F_vegetation = rad.F_sky * e_vars.shade
    F_sky = rad.F_sky - F_vegetation
    F_ground = rad.F_ground

    area_silhouette = silhouette_area(o.body, rad.solar_orientation)
    area_total = total_area(o.body)
    area_dorsal = area_total * (1 - rad.ventral_fraction)
    area_ventral = area_total * rad.ventral_fraction * (1 - cond_ex.conduction_fraction)
    area_skin = skin_area(o.body)
    area_conduction = area_total * cond_ex.conduction_fraction
    area_evaporation = evaporation_area(o.body)
    area_convection = area_total * (1 - cond_ex.conduction_fraction)

    (; Q_solar, Q_direct, Q_solar_sky, Q_solar_substrate) =
        solar(
            α_body_dorsal,
            α_body_ventral, 
            area_silhouette,
            area_dorsal, 
            area_ventral, 
            F_ground, 
            F_sky, 
            e_pars.α_ground, 
            e_vars.shade,
            e_vars.zenith_angle, 
            e_vars.global_radiation, 
            e_vars.diffuse_fraction
            )

    Q_ventral = Q_solar_substrate

    # set infrared environment
    T_vegetation = e_vars.T_air_reference # assume vegetation casting shade is at reference (e.g. 1.2m or 2m) air temperature (deg C)

    for side in 1:2

        # Calculating solar intensity entering insulation. This will depend on whether we are calculating 
        # the insulation temperature for the dorsal side or the ventral side. The dorsal side will have 
        # solar inputs from the direct beam hitting the silhouette area as well as diffuse solar scattered 
        # from the sky. The ventral side will have diffuse solar scattered off the substrate.
        # Resetting config factors and solar depending on whether the dorsal side (side=1) or 
        # ventral side (side=2) is being estimated.
        if Q_solar > 0.0u"W"
            if side == 1
                F_sky = rad.F_sky * 2.0 # proportion of upward view that is sky
                F_vegetation = rad.F_vegetation * 2.0 # proportion of upward view that is vegetation (shade)
                F_ground = 0.0
                F_bush = 0.0
                Q_sol = 2.0 * Q_direct + ((Q_solar_sky / rad.F_sky) * F_sky) # direct x 2 because assuming sun in both directions, and unadjusting Q_solar_sky for config factor imposed in SOLAR_ENDO and back to new larger one in both directions
            else
                F_sky = 0.0
                F_vegetation = 0.0
                F_ground = rad.F_ground * 2.0
                F_bush = rad.F_bush * 2.0
                Q_sol = (Q_ventral / (1.0 - rad.F_sky - rad.F_vegetation)) * (1.0 - (2.0 * cond_ex.conduction_fraction)) # unadjust by config factor imposed in SOLAR_ENDO to have it coming in both directions, but also cutting down according to fractional area conducting to ground (in both directions)
            end
        else
            Q_sol = 0.0u"W"
            if side == 1
                F_sky = rad.F_sky * 2.0
                F_vegetation = rad.F_vegetation * 2.0
                F_ground = 0.0
                F_bush = 0.0
            else
                F_sky = 0.0
                F_vegetation = 0.0
                F_ground = rad.F_ground * 2.0
                F_bush = rad.F_bush
            end
        end

        if side == 1
            ϵ_body = rad.ϵ_body_dorsal
        else
            ϵ_body = rad.ϵ_body_ventral
        end

        # set fur depth and conductivity
        # index for effective_conductivities, is the average (1), front/dorsal (2), back/ventral(3) of the body part
        if Q_solar > 0.0u"W" || (insulation_depths[2] != insulation_depths[3])
            if side == 1
                insulation = Fur(insulation_depths[2], fibre_diameters[2], fibre_densities[2])
            else
                insulation = Fur(insulation_depths[3], fibre_diameters[3], fibre_densities[3])
            end
        else
            insulation = Fur(insulation_depths[1], fibre_diameters[1], fibre_densities[1])
        end

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
            A_conduction = A_total * (cond_ex.conduction_fraction * 2)
            cd = (A_conduction * e_vars.k_substrate) / e_pars.conduction_depth # assume conduction happens from 2.5 cm depth
        else  # doing dorsal side, no conduction. No need to adjust areas used for convection.
            A_conduction = 0.0u"m^2"
            cd = 0.0u"W/K"
        end

        # package up inputs
        geom_vars = (; 
                side,
                cd,
                ventral_fraction = rad.ventral_fraction,
                conduction_fraction = cond_ex.conduction_fraction,
                longwave_depth_fraction = ins.longwave_depth_fraction)
        env_vars = (; 
                F_sky,
                F_ground,
                F_bush,
                F_vegetation,
                fluid = e_pars.fluid,
                T_air = e_vars.T_air,
                T_ground = e_vars.T_ground,
                T_bush = e_vars.T_bush,
                T_vegetation,
                T_sky = e_vars.T_sky,
                T_substrate = e_vars.T_substrate,
                rh = e_vars.rh,
                wind_speed = e_vars.wind_speed,
                P_atmos = e_vars.P_atmos,
                Q_solar = Q_sol, 
                fO2 = e_pars.fO2,
                fCO2 = e_pars.fCO2,
                fN2 = e_pars.fN2,
                convection_enhancement = e_pars.convection_enhancement,
                )
        traits = (; 
                T_core = metab.T_core,
                k_flesh = cond_in.k_flesh,
                k_fat = cond_in.k_fat,
                ϵ_body,
                skin_wetness = evap.skin_wetness,
                insulation_wetness = evap.insulation_wetness,
                bare_skin_fraction = evap.bare_skin_fraction,
                eye_fraction = evap.eye_fraction,
                insulation_conductivity,
                #ψ_body = hyd.ψ_body,
                )
        if side == 1
            insulation_side = InsulationParameters(;
            insulation_conductivity_dorsal = ins.insulation_conductivity_dorsal,
            insulation_conductivity_ventral = ins.insulation_conductivity_dorsal,
            fibre_diameter_dorsal = ins.fibre_diameter_dorsal,
            fibre_diameter_ventral = ins.fibre_diameter_dorsal,
            fibre_length_dorsal = ins.fibre_length_dorsal,
            fibre_length_ventral = ins.fibre_length_dorsal,
            insulation_depth_dorsal = ins.insulation_depth_dorsal,
            insulation_depth_ventral = ins.insulation_depth_dorsal,    
            max_insulation_depth_dorsal = ins.max_insulation_depth_dorsal,
            max_insulation_depth_ventral = ins.max_insulation_depth_dorsal,
            fibre_density_dorsal = ins.fibre_density_dorsal,
            fibre_density_ventral = ins.fibre_density_dorsal,
            insulation_reflectance_dorsal = ins.insulation_reflectance_dorsal,
            insulation_reflectance_ventral = ins.insulation_reflectance_dorsal,
            insulation_depth_compressed = ins.insulation_depth_compressed,
            fibre_conductivity = ins.fibre_conductivity,
            longwave_depth_fraction = ins.longwave_depth_fraction,
            )
        else
            insulation_side = InsulationParameters(;
            insulation_conductivity_dorsal = ins.insulation_conductivity_ventral,
            insulation_conductivity_ventral = ins.insulation_conductivity_ventral,
            fibre_diameter_dorsal = ins.fibre_diameter_ventral,
            fibre_diameter_ventral = ins.fibre_diameter_ventral,
            fibre_length_dorsal = ins.fibre_length_ventral,
            fibre_length_ventral = ins.fibre_length_ventral,
            insulation_depth_dorsal = ins.insulation_depth_ventral,
            insulation_depth_ventral = ins.insulation_depth_ventral,    
            max_insulation_depth_dorsal = ins.max_insulation_depth_ventral,
            max_insulation_depth_ventral = ins.max_insulation_depth_ventral,
            fibre_density_dorsal = ins.fibre_density_ventral,
            fibre_density_ventral = ins.fibre_density_ventral,
            insulation_reflectance_dorsal = ins.insulation_reflectance_ventral,
            insulation_reflectance_ventral = ins.insulation_reflectance_ventral,
            insulation_depth_compressed = ins.insulation_depth_compressed,
            fibre_conductivity = ins.fibre_conductivity,
            longwave_depth_fraction = ins.longwave_depth_fraction,
            )
        end            
        insulation_out = insulation_properties(; 
                        insulation = insulation_side, insulation_temperature = T_insulation, 
                        ventral_fraction = rad.ventral_fraction)

        # call simulsol
        simulsol_out[side,] = simulsol(; 
                                T_skin, 
                                T_insulation, 
                                simulsol_tolerance = m.simulsol_tolerance,
                                geometry_pars,
                                insulation_pars = insulation_side, 
                                insulation_out,
                                geom_vars,
                                env_vars,
                                traits)
        
        T_skin = simulsol_out[side,].T_skin # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
        T_insulation = simulsol_out[side,].T_insulation # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
        simulsol_tolerance = simulsol_out[side,].tolerance # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)

    end

    T_skin_max = max(simulsol_out[1].T_skin, simulsol_out[2].T_skin)

    # zbrent and respfun

    # Now compute a weighted mean heat generation for all the parts/components = (dorsal value *(F_sky+F_vegetation))+(ventral value*F_ground)
    gen_d = simulsol_out[1].Q_gen_net
    gen_v = simulsol_out[2].Q_gen_net
    dmult = rad.F_sky + rad.F_vegetation
    vmult = 1 - dmult # assume that reflectivity of veg below equals reflectivity of soil so vmult left as 1 - dmult
    x = gen_d * dmult + gen_v * vmult # weighted estimate of metabolic heat generation
    Q_sum = x

    # reset configuration factors
    F_bush = rad.F_bush # nearby bush
    F_sky = rad.F_sky # sky
    F_ground = rad.F_ground # ground
    F_vegetation = rad.F_vegetation # vegetation

    # lung temperature and temperature of exhaled air
    T_skin = (simulsol_out[1].T_skin + simulsol_out[2].T_skin) * 0.5 # TODO weight it by dorsal/ventral fractions?
    T_insulation = (simulsol_out[1].T_insulation + simulsol_out[2].T_insulation) * 0.5
    T_lung = (metab.T_core + T_skin) * 0.5 # average of skin and core
    T_air_exit = min(e_vars.T_air + resp.Δ_breath, T_lung) # temperature of exhaled air, deg C

    if m.respire
        # now guess for metabolic rate that balances the heat budget while allowing metabolic rate
        # to remain at or above Q_basal, via root-finder ZBRENT
        Q_min = metab.Q_metabolism
        Q_m1 = metab.Q_metabolism * (-2.)
        Q_m2 = metab.Q_metabolism * 10.
        if T_skin_max >= metab.T_core
            Q_m2 = metab.Q_metabolism * 1.01
        end

        Q_gen = find_zero(x -> respiration(; 
                                Q_metab = x,
                                Q_min,
                                Q_sum,
                                T_lung,
                                mass = geometry_pars.shape.mass,
                                rq = resp.rq,
                                fO2_extract = resp.fO2_extract,
                                rh_exit = resp.rh_exit,
                                T_air_exit,
                                pant = resp.pant,                                
                                T_air = e_vars.T_air,
                                P_atmos = e_vars.P_atmos,                                
                                rh = e_vars.rh,
                                fO2 = e_pars.fO2,
                                fN2 = e_pars.fN2,
                                fCO2 = e_pars.fCO2,                                
                                O2conversion = Kleiber1961()).balance,
                                (Q_m1, Q_m2),
                                Roots.Brent(),
                                atol = m.resp_tolerance * metab.Q_metabolism,
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
            fO2=e_pars.fO2,
            fN2=e_pars.fN2,
            fCO2=e_pars.fCO2,
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
    Q_convection_dorsal = simulsol_out[1].Q_convection   # convection (W)
    Q_conduction_dorsal = simulsol_out[1].Q_conduction   # conduction (W)
    Q_evap_skin_dorsal = simulsol_out[1].Q_evap_skin   # cutaneous evaporation (W)
    Q_longwave_dorsal = simulsol_out[1].Q_longwave   # radiative loss (W)
    Q_solar_dorsal = simulsol_out[1].Q_solar   # solar (W)
    Q_evap_insulation_dorsal = simulsol_out[1].Q_evap_insulation  # fur evaporation (W)
    ntry_dorsal = simulsol_out[1].ntry  # attempts
    success_dorsal = simulsol_out[1].success  # success?
    k_insulation_dorsal = simulsol_out[1].k_insulation  # fur conductivity? (same as FORTRAN)

    T_insulation_ventral = simulsol_out[2].T_insulation   # feather/fur-air interface temp (°C)
    T_skin_ventral = simulsol_out[2].T_skin   # average skin temp
    Q_convection_ventral = simulsol_out[2].Q_convection   # convection (W)
    Q_conduction_ventral = simulsol_out[2].Q_conduction   # conduction (W)
    Q_evap_skin_ventral = simulsol_out[2].Q_evap_skin   # cutaneous evaporation (W)
    Q_longwave_ventral = simulsol_out[2].Q_longwave   # radiative loss (W)
    Q_solar_ventral = simulsol_out[2].Q_solar   # solar (W)
    Q_evap_insulation_ventral = simulsol_out[2].Q_evap_insulation  # fur evaporation (W)
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
        J_air_in = respiration_out.J_air_in
        J_air_out = respiration_out.J_air_out
        J_H2O_in = respiration_out.J_H2O_in
        J_H2O_out = respiration_out.J_H2O_out
        J_O2_in = respiration_out.J_O2_in
        J_O2_out = respiration_out.J_O2_out
        J_CO2_in = respiration_out.J_CO2_in
        J_CO2_out = respiration_out.J_CO2_out
        J_N2_in = respiration_out.J_N2_in
        J_N2_out = respiration_out.J_N2_out
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
    fat_volume = fat_mass / fat.density
    volume = geometry_pars.geometry.volume
    flesh_volume = volume - fat_volume
    area_total = total_area(geometry_pars)
    area_skin = skin_area(geometry_pars)
    radius_insulation = insulation_radius(geometry_pars)
    characteristic_dimension = 2 * radius_insulation
    @show characteristic_dimension
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
    Q_rad_out_dorsal = 2 * rad.F_sky * σ * rad.ϵ_body_dorsal * area_radiant_dorsal * T_surface_dorsal^4
    Q_rad_in_dorsal = -Q_longwave_dorsal + Q_rad_out_dorsal

    # infrared ventral
    Q_rad_out_ventral = 2 * rad.F_ground * σ * rad.ϵ_body_ventral * area_radiant_ventral * T_surface_ventral^4
    Q_rad_in_ventral = -Q_longwave_ventral + Q_rad_out_ventral
    # energy flows
    Q_solar = Q_solar_dorsal * dmult + Q_solar_ventral * vmult
    Q_longwave_in = Q_rad_in_dorsal * dmult + Q_rad_in_ventral * vmult
    Q_evaporation = Q_evap_skin_dorsal * dmult + Q_evap_skin_ventral * vmult +
                    Q_evap_insulation_dorsal * dmult + Q_evap_insulation_ventral * vmult + Q_resp
    Q_longwave_out = Q_rad_out_dorsal * dmult + Q_rad_out_ventral * vmult
    Q_convection = Q_convection_dorsal * dmult + Q_convection_ventral * vmult
    Q_conduction = Q_conduction_dorsal * dmult + Q_conduction_ventral * vmult

    insulation_out = insulation_properties(; insulation=ins,
        insulation_temperature=T_insulation * 0.7 + T_skin * 0.3, rad.ventral_fraction)
    k_insulation_effective = insulation_out.effective_conductivities[1]
    k_insulation_compressed = insulation_out.insulation_conductivity_compressed
    T_skin = T_skin_dorsal * dmult + T_skin_ventral * vmult
    T_insulation = T_insulation_dorsal * dmult + T_insulation_ventral * vmult

    thermoregulation = (;
        metab.T_core, T_skin, T_insulation, T_lung, T_skin_dorsal, T_skin_ventral, T_insulation_dorsal, 
        T_insulation_ventral, shape_b=o.body.shape.b, pant=resp.pant, skin_wetness=evap.skin_wetness, 
        k_flesh=cond_in.k_flesh, k_insulation_effective, k_insulation_dorsal, k_insulation_ventral, 
        k_insulation_compressed, ins.insulation_depth_dorsal, ins.insulation_depth_ventral, metab.q10
    )

    morphology = (;
        area_total, area_skin, area_evaporation, area_convection, area_conduction, area_silhouette,
        F_sky, F_ground, volume, flesh_volume, 
        # TODO make it so that this is the proper characteristic dimension computed via the geometry functions
        #characteristic_dimension = geometry_pars.geometry.characteristic_dimension,
        characteristic_dimension,
     fat_mass, geometry_pars.geometry.length...
    )

    energy_fluxes = (;
        Q_solar, Q_longwave_in, Q_gen, Q_evaporation, Q_longwave_out, Q_convection, Q_conduction, balance,
        ntry=max(ntry_dorsal, ntry_ventral), success=all((success_dorsal, success_ventral))
    )

    mass_fluxes = (; V_air, V_O2_STP, m_evap, m_resp, m_sweat, J_H2O_in, J_H2O_out, J_O2_in, J_O2_out,
        J_CO2_in, J_CO2_out, J_N2_in, J_N2_out, J_air_in, J_air_out
    )
    return (; thermoregulation, morphology, energy_fluxes, mass_fluxes)
end