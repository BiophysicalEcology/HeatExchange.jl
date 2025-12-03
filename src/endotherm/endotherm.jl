Base.@kwdef struct EndoModelPars{TM,TR,RE,TP,ST,BT} <: AbstractModelParameters
    thermoregulation_mode::TM = Param(1)
    thermoregulate::TR =        Param(true)
    respire::RE =               Param(trie)
    torpor::TP =                Param(false)
    simulsol_tolerance::ST =    Param(1e-3u"K")
    resp_tolerance::BT =        Param(1e-5u"K")
end

function endotherm(; model_pars, bodyshape, body_pars, integument_pars, insulation_pars, 
        physio_pars, thermoreg_pars, thermoreg_vars, environmental_pars, organism_vars, environmental_vars)
        shade = 0.0 # TODO how to deal with shade?
    # unpack parameters and variables
    (; thermoregulation_mode, thermoregulate, respire, torpor, simulsol_tolerance, resp_tolerance) = model_pars
    (; mass, ρ_flesh, ρ_fat, shape_b, shape_c) = body_pars
    (; α_body_dorsal, α_body_ventral, ϵ_body_dorsal, ϵ_body_ventral, F_sky,
        F_ground, F_vegetation, F_bush, eye_fraction, bare_skin_fraction,
        ventral_fraction) = integument_pars
    (; insulation_conductivity_dorsal, insulation_conductivity_ventral, 
        fibre_diameter_dorsal, fibre_diameter_ventral, fibre_length_dorsal,
        fibre_length_ventral, max_insulation_depth_dorsal, max_insulation_depth_ventral, 
        fibre_density_dorsal, fibre_density_ventral, 
        insulation_reflectance_dorsal, insulation_reflectance_ventral, 
        insulation_depth_compressed, fibre_conductivity, longwave_depth_fraction) = insulation_pars
    (; Q_minimum, q10, k_fat, fO2_extract, rq, Δ_breath, rh_exit) = physio_pars
    (; insulation_step, shape_b_step, shape_b_max, T_core_max, T_core_min, 
        T_core_step, k_flesh_step, k_flesh_max, pant_step, pant_multiplier, pant_max,
        skin_wetness_step, skin_wetness_max) = thermoreg_pars
    (; insulation_depth_dorsal, insulation_depth_ventral, fat_fraction, conduction_fraction, 
        k_flesh, T_core_target, pant, skin_wetness, insulation_wetness) = thermoreg_vars    
    (; α_ground, ϵ_ground, ϵ_sky, elevation, fluid, fN2, fO2, fCO2, convection_enhancement) = 
        environmental_pars
    (; T_core, T_skin, T_insulation, T_lung, ψ_org) = organism_vars
    (; T_air, T_air_reference, T_sky, T_ground, T_substrate, T_bush, T_vegetation, rh, 
        wind_speed, P_atmos, zenith_angle, k_substrate, global_radiation, 
        diffuse_fraction) = environmental_vars

    # check if starting piloerect
    if insulation_step > 0.0
        insulation_depth_dorsal = max_insulation_depth_dorsal
        insulation_depth_ventral = max_insulation_depth_ventral
        thermoreg_vars.insulation_depth_dorsal = insulation_depth_dorsal
        thermoreg_vars.insulation_depth_ventral = insulation_depth_ventral
    end
    
    ϵ_body = ϵ_body_dorsal # TODO use both dorsal and ventral
    Q_minimum_ref = Q_minimum
    T_core_ref = T_core_target
    insulation_depth_dorsal_ref = insulation_depth_dorsal
    insulation_depth_ventral_ref = insulation_depth_ventral
    T_air_reference = T_air
    F_ground_reference = F_ground
    F_sky_reference = F_sky
    F_vegetation_reference = F_vegetation
    F_bush_reference = F_bush
    Q_gen = 0.0u"W"
    pant_cost = 0.0u"W"
    shape_b_last = shape_b
    k_flesh_last = k_flesh
    T_core_last = T_core
    pant_last = pant
    skin_wetness_last = skin_wetness
    simulsol_out = Vector{NamedTuple}(undef, 2) # TODO preallocate
    respiration_out = Vector{NamedTuple}(undef, 1) # TODO preallocate
    geometry_pars = nothing
    dmult = nothing
    vmult = nothing
    T_lung = nothing
    area_silhouette = nothing
    area_total = nothing
    area_skin = nothing
    area_convection = nothing
    area_evaporation = nothing
    area_conduction = nothing
    fat = nothing
    fur = nothing
    insulation_out = nothing

    while Q_gen < Q_minimum * 0.995
        insulation_temperature = T_insulation * 0.7 + T_skin * 0.3
        insulation_out = insulation_properties(; insulation=insulation_pars, insulation_temperature, 
            ventral_fraction, insulation_depth_dorsal, insulation_depth_ventral)
        fibre_diameters = insulation_out.fibre_diameters # fur diameter array, mean, dorsal, ventral (m)
        fibre_densities = insulation_out.fibre_densities # fur density array, mean, dorsal, ventral (1/m2)
        insulation_depths = insulation_out.insulation_depths # fur depth array, mean, dorsal, ventral (m)

        fat = Fat(fat_fraction, ρ_fat)
        fur = Fur(insulation_depths[1], fibre_diameters[1], fibre_densities[1])
        composite_insulation = (fat, fur)

        # geometric properties
        geometry_pars = Body(bodyshape, CompositeInsulation(fur, fat))

        α_body_dorsal = 1 - insulation_reflectance_dorsal #solar absorptivity of dorsal fur (fractional, 0-1)
        α_body_ventral = 1 - insulation_reflectance_ventral # solar absorptivity of ventral fur (fractional, 0-1)

        # correct F_sky for vegetaion overhead
        F_vegetation = F_sky_reference * shade
        F_sky = F_sky_reference - F_vegetation
        F_ground = F_ground_reference

        area_silhouette = silhouette_area(geometry_pars, thermoreg_vars.solar_orientation)        
        area_total = get_total_area(geometry_pars) # total area, m2
        area_skin = get_skin_area(geometry_pars) # total area, m2
        area_conduction = area_total * conduction_fraction # area of skin for convection/evaporation (total skin area - hair area), m2
        area_evaporation = get_evaporation_area(geometry_pars)
        area_convection = area_total * (1 - conduction_fraction)

        (; Q_solar, Q_direct, Q_solar_sky, Q_solar_substrate) = 
            solar(α_body_dorsal, α_body_ventral, area_silhouette, 
            area_total, area_conduction, F_ground, F_sky, α_ground, shade, 
            zenith_angle, global_radiation, diffuse_fraction)

        #Q_dorsal = Q_direct + Q_solar_sky
        Q_ventral = Q_solar_substrate
        #Q_diffuse = Q_solar_sky + Q_solar_substrate

        # # convection
        # (; Q_conv, hc, hd, Sh,
        #     Q_free, Q_forc,
        #     hc_free, hc_forc,
        #     Sh_free, Sh_forc,
        #     hd_free, hd_forc) = convection(geometry_pars, area_convection, T_air, T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)

        # set infrared environment
        T_vegetation = T_air_reference # assume vegetation casting shade is at reference (e.g. 1.2m or 2m) air temperature (deg C)

        for side in 1:2

            # Calculating solar intensity entering insulation. This will depend on whether we are calculating the fur
            # temperature for the dorsal side or the ventral side. The dorsal side will have solar inputs from
            # the direct beam hitting the silhouette area as well as diffuse solar scattered from the sky.
            # The ventral side will have diffuse solar scattered off the substrate.
            # Resetting config factors and solar depending on whether the dorsal side (side=1) or ventral side (side=2) is being estimated.
            if Q_solar > 0.0u"W"
                if side == 1
                    F_sky = F_sky_reference * 2.0 # proportion of upward view that is sky
                    F_vegetation = F_vegetation_reference * 2.0 # proportion of upward view that is vegetation (shade)
                    F_ground = 0.0
                    F_bush = 0.0
                    Q_solar = 2.0 * Q_direct + ((Q_solar_sky / F_sky_reference) * F_sky) # direct x 2 because assuming sun in both directions, and unadjusting Q_solar_sky for config factor imposed in SOLAR_ENDO and back to new larger one in both directions
                else
                    F_sky = 0.0
                    F_vegetation = 0.0
                    F_ground = F_ground_reference * 2.0
                    F_bush = F_bush_reference * 2.0
                    Q_solar = (Q_ventral / (1.0 - F_sky_reference - F_vegetation_reference)) * (1.0 - (2.0 * conduction_fraction)) # unadjust by config factor imposed in SOLAR_ENDO to have it coming in both directions, but also cutting down according to fractional area conducting to ground (in both directions)
                end
            else
                Q_solar = 0.0u"W"
                if side == 1
                    F_sky = F_sky_reference * 2.0
                    F_vegetation = F_vegetation_reference * 2.0
                    F_ground = 0.0
                    F_bush = 0.0
                else
                    F_sky = 0.0
                    F_vegetation = 0.0
                    F_ground = F_ground_reference * 2.0
                    F_bush = F_bush_reference
                end
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
                insulation_conductivity = insulation_conductivity_dorsal
            else
                insulation_conductivity = insulation_conductivity_ventral                
            end
            composite_insulation = (fat, fur)
            geometry_pars = Body(bodyshape, CompositeInsulation(insulation, fat))
            # volume, accounting for any subcutaneous fat
            #volume = geometry_pars.geometry.volume - fat_fraction * mass / ρ_fat

            A_total = get_total_area(geometry_pars)
            r_skin = get_r_skin(geometry_pars) # body radius (including fat), m
            #r_flesh = get_r_flesh(geometry_pars) # body radius flesh only (no fat), m
            r_insulation = r_skin + insulation.thickness # body radius including fur, m
            #diameter = 2 * r_insulation # diameter, m
            #r_radiation = r_skin + (longwave_depth_fraction * insulation.thickness) # effective radiation radius, m
            if geometry_pars.shape isa Cylinder || geometry_pars.shape isa Sphere
                r_compressed = r_skin + insulation_depth_compressed
            else
                r_compressed = r_insulation # Note that this value is never used if conduction not being modeled, but need to have a value for the calculations
            end

            #LEN = ALENTH # length, m

            # Calculating the "cd" variable: Qcond = cd(Tskin-Tsub), where cd = Conduction area*ksub/subdepth
            if side == 2 # doing ventral side, add conduction
                A_conduction = A_total * (conduction_fraction * 2)
                cd = (A_conduction * k_substrate) / 0.025u"m" # assume conduction happens from 2.5 cm depth
            else  # doing dorsal side, no conduction. No need to adjust areas used for convection.
                A_conduction = 0.0u"m^2"
                cd = 0.0u"W/K"
            end

            # package up inputs
            geom_vars = (; side, cd, conduction_fraction, longwave_depth_fraction)
            env_vars = (; fluid, T_air, T_ground, T_bush, T_vegetation, T_sky, T_substrate,
                rh, wind_speed, P_atmos, F_sky, F_ground, F_bush, F_vegetation, Q_solar, fO2, fCO2, fN2, 
                convection_enhancement)
            traits = (; T_core, k_flesh, k_fat, ϵ_body, skin_wetness,
                insulation_wetness, bare_skin_fraction, eye_fraction, insulation_conductivity, 
                insulation_depth_dorsal, insulation_depth_ventral)
            simulsol_out[side,] = simulsol(; T_skin, T_insulation, simulsol_tolerance, geometry_pars, 
                insulation_pars, insulation_out, geom_vars, env_vars, traits)
            T_skin = simulsol_out[side,].T_skin # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
            T_insulation = simulsol_out[side,].T_insulation # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)
            simulsol_tolerance = simulsol_out[side,].tolerance # TODO check if connecting this value across runs per side is a good idea (happens in Fortran)

        end

        T_skin_max = max(simulsol_out[1].T_skin, simulsol_out[2].T_skin)

        # zbrent and respfun

        # Now compute a weighted mean heat generation for all the parts/components = (dorsal value *(F_sky+F_vegetation))+(ventral value*F_ground)
        gen_d = simulsol_out[1].Q_gen_net
        gen_v = simulsol_out[2].Q_gen_net
        dmult = F_sky_reference + F_vegetation_reference
        vmult = 1 - dmult # assume that reflectivity of veg below equals reflectivity of soil so vmult left as 1 - dmult
        x = gen_d * dmult + gen_v * vmult # weighted estimate of metabolic heat generation
        Q_sum = x

        # reset configuration factors
        F_bush = F_bush_reference # nearby bush
        F_sky = F_sky_reference # sky
        F_ground = F_ground_reference # ground
        F_vegetation = F_vegetation_reference # vegetation

        # lung temperature and temperature of exhaled air
        T_skin = (simulsol_out[1].T_skin + simulsol_out[2].T_skin) * 0.5 # TODO weight it by dorsal/ventral fractions?
        T_insulation = (simulsol_out[1].T_insulation + simulsol_out[2].T_insulation) * 0.5
        T_lung = (T_core + T_skin) * 0.5 # average of skin and core
        T_air_exit = min(T_air + Δ_breath, T_lung) # temperature of exhaled air, deg C

        if respire
            # now guess for metabolic rate that balances the heat budget while allowing metabolic rate
            # to remain at or above Q_basal, via root-finder ZBRENT
            Q_min = Q_minimum
            Q_m1 = Q_minimum * (-2.)
            Q_m2 = Q_minimum * 10.
            if T_skin_max >= T_core
                Q_m2 = Q_minimum * 1.01
            end

            Q_gen = find_zero(x -> respiration(; Q_metab = x, T_air, fO2, 
                fN2, fCO2, P_atmos, Q_min, rq, T_lung, mass, fO2_extract, rh, rh_exit, 
                T_air_exit, pant, Q_sum, O2conversion=Kleiber1961()).balance, (Q_m1, Q_m2), Roots.Brent(), atol = resp_tolerance * Q_minimum)

            respiration_out = respiration(; Q_metab = Q_gen, T_air, fO2, 
                fN2, fCO2, P_atmos, Q_min, rq, T_lung, mass, fO2_extract, rh, rh_exit, 
                T_air_exit, pant, Q_sum, O2conversion=Kleiber1961())

            Q_gen = respiration_out.Q_gen # Q_gen_net
        else
            Q_gen = Q_sum
        end
        # store last-step values
        shape_b_last = shape_b
        k_flesh_last = k_flesh
        T_core_last = T_core
        pant_last = pant
        skin_wetness_last = skin_wetness

        if thermoregulate

            if (insulation_depth_dorsal > insulation_depth_dorsal_ref) &&
               (insulation_depth_ventral > insulation_depth_ventral_ref)
                insulation_depth_dorsal = max(insulation_depth_dorsal_ref, 
                    insulation_depth_dorsal - insulation_step * fibre_length_dorsal)
                insulation_depth_ventral = max(insulation_depth_ventral_ref, 
                    insulation_depth_ventral - insulation_step * fibre_length_ventral)
                thermoreg_vars.insulation_depth_dorsal = insulation_depth_dorsal
                thermoreg_vars.insulation_depth_ventral = insulation_depth_ventral    
            else
                insulation_depth_dorsal = insulation_depth_dorsal_ref
                insulation_depth_ventral = insulation_depth_ventral_ref
                thermoreg_vars.insulation_depth_dorsal = insulation_depth_dorsal
                thermoreg_vars.insulation_depth_ventral = insulation_depth_ventral
                if shape_b < shape_b_max
                    shape_b += shape_b_step
                    geometry_pars.shape.b = shape_b
                else
                    shape_b = shape_b_max
                    geometry_pars.shape.b = shape_b

                    if k_flesh < k_flesh_max
                        k_flesh += k_flesh_step
                        thermoreg_vars.k_flesh = k_flesh
                    else
                        k_flesh = k_flesh_max
                        thermoreg_vars.k_flesh = k_flesh

                        if T_core < T_core_max
                            T_core += T_core_step
                            thermoreg_vars.T_core_target = T_core

                            q10mult = q10^((ustrip(u"K", (T_core - T_core_ref))) / 10)

                            if (thermoregulation_mode >= 2) && (pant < pant_max)
                                pant += pant_step
                                thermoreg_vars.pant = pant
                                pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                            (pant_multiplier - 1) * Q_minimum_ref

                                if thermoregulation_mode == 3
                                    skin_wetness += skin_wetness_step
                                    thermoreg_vars.skin_wetness = skin_wetness
                                    if skin_wetness > skin_wetness_max
                                        skin_wetness = skin_wetness_max
                                        thermoreg_vars.skin_wetness = skin_wetness 
                                    end
                                end
                            end

                            Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                        else
                            T_core = T_core_max
                            thermoreg_vars.T_core_target = T_core                            
                            q10mult = q10^((ustrip(u"K", (T_core - T_core_ref))) / 10)

                            if pant < pant_max
                                pant += pant_step
                                thermoreg_vars.pant = pant
                                pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                            (pant_multiplier - 1) * Q_minimum_ref
                                Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                                if thermoregulation_mode == 3
                                    skin_wetness += skin_wetness_step
                                    #thermoreg_vars.skin_wetness = skin_wetness
                                    if skin_wetness > skin_wetness_max
                                        skin_wetness = skin_wetness_max
                                    end
                                end
                            else
                                pant = pant_max
                                #thermoreg_vars.pant = pant
                                pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                            (pant_multiplier - 1) * Q_minimum_ref

                                Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                                skin_wetness += skin_wetness_step
                                #thermoreg_vars.skin_wetness = skin_wetness
                                if (skin_wetness > skin_wetness_max) || (skin_wetness_step <= 0)
                                    skin_wetness = skin_wetness_max
                                    #thermoreg_vars.skin_wetness = skin_wetness
                                    return
                                end
                            end
                        end
                    end
                end
            end
        else
            break
        end
    end

    # collate output

    # simusol outputs
    T_insulation_dorsal      = simulsol_out[1].T_insulation   # feather/fur-air interface temp (°C)
    T_skin_dorsal            = simulsol_out[1].T_skin   # average skin temp
    Q_convection_dorsal      = simulsol_out[1].Q_convection   # convection (W)
    Q_conduction_dorsal      = simulsol_out[1].Q_conduction   # conduction (W)
    Q_evap_skin_dorsal       = simulsol_out[1].Q_evap_skin   # cutaneous evaporation (W)
    Q_longwave_dorsal        = simulsol_out[1].Q_longwave   # radiative loss (W)
    Q_solar_dorsal           = simulsol_out[1].Q_solar   # solar (W)
    Q_evap_insulation_dorsal = simulsol_out[1].Q_evap_insulation  # fur evaporation (W)
    ntry_dorsal              = simulsol_out[1].ntry  # attempts
    success_dorsal           = simulsol_out[1].success  # success?
    k_insulation_dorsal      = simulsol_out[1].k_insulation  # fur conductivity? (same as FORTRAN)

    T_insulation_ventral      = simulsol_out[2].T_insulation   # feather/fur-air interface temp (°C)
    T_skin_ventral            = simulsol_out[2].T_skin   # average skin temp
    Q_convection_ventral      = simulsol_out[2].Q_convection   # convection (W)
    Q_conduction_ventral      = simulsol_out[2].Q_conduction   # conduction (W)
    Q_evap_skin_ventral       = simulsol_out[2].Q_evap_skin   # cutaneous evaporation (W)
    Q_longwave_ventral        = simulsol_out[2].Q_longwave   # radiative loss (W)
    Q_solar_ventral           = simulsol_out[2].Q_solar   # solar (W)
    Q_evap_insulation_ventral = simulsol_out[2].Q_evap_insulation  # fur evaporation (W)
    ntry_ventral              = simulsol_out[2].ntry  # attempts
    success_ventral           = simulsol_out[2].success  # success?
    k_insulation_ventral      = simulsol_out[2].k_insulation  # fur conductivity? (same as FORTRAN)

    # respiration outputs
    if respire
        balance   = respiration_out.balance
        Q_resp    = respiration_out.Q_resp
        m_resp    = respiration_out.m_resp
        V_air     = respiration_out.V_air
        V_O2_STP  = respiration_out.V_O2_STP
        J_air_in  = respiration_out.J_air_in
        J_air_out = respiration_out.J_air_out
        J_H2O_in  = respiration_out.J_H2O_in
        J_H2O_out = respiration_out.J_H2O_out
        J_O2_in   = respiration_out.J_O2_in
        J_O2_out  = respiration_out.J_O2_out
        J_CO2_in  = respiration_out.J_CO2_in
        J_CO2_out = respiration_out.J_CO2_out
        J_N2_in   = respiration_out.J_N2_in
        J_N2_out  = respiration_out.J_N2_out
    else
        balance   = nothing
        Q_resp    = 0.0u"W"
        m_resp    = nothing
        V_air     = nothing
        V_O2_STP  = nothing
        J_air_in  = nothing
        J_air_out = nothing
        J_H2O_in  = nothing
        J_H2O_out = nothing
        J_O2_in   = nothing
        J_O2_out  = nothing
        J_CO2_in  = nothing
        J_CO2_out = nothing
        J_N2_in   = nothing
        J_N2_out  = nothing
    end
    # evaporation calculations
    L_v   = enthalpy_of_vaporisation(T_air)
    m_sweat  = u"g/hr"((Q_evap_skin_dorsal + Q_evap_skin_ventral) * 0.5 / L_v)
    if respire
        m_evap   = u"g/hr"(m_resp + m_sweat)
    else 
        m_evap = u"g/hr"(m_sweat)
    end

    # geometric outputs
    fat_mass = geometry_pars.shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = geometry_pars.geometry.volume
    flesh_volume = volume - fat_volume
    area_total = get_total_area(geometry_pars)
    area_skin = get_skin_area(geometry_pars)

    # radiation outputs

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    if insulation_out.insulation_test > 0.0u"m"
        T_surface_dorsal = T_insulation_dorsal
        T_surface_ventral = T_insulation_ventral
        area_radiant_dorsal = area_total
        area_radiant_ventral = area_total * (1 - conduction_fraction)
    else
        T_surface_dorsal = T_skin_dorsal
        T_surface_ventral = T_skin_ventral
        area_radiant_dorsal = area_skin
        area_radiant_ventral = area_skin #* (1 - conduction_fraction) TODO make conduction occur when no insulation
    end
    # infrared dorsal
    Q_rad_out_dorsal = 2 * F_sky * σ * ϵ_body_dorsal * area_radiant_dorsal * T_surface_dorsal^4
    Q_rad_in_dorsal  = -Q_longwave_dorsal + Q_rad_out_dorsal

    # infrared ventral
    Q_rad_out_ventral = 2 * F_ground * σ * ϵ_body_ventral * area_radiant_ventral * T_surface_ventral^4
    Q_rad_in_ventral  = -Q_longwave_ventral + Q_rad_out_ventral
    # energy flows
    Q_solar = Q_solar_dorsal * dmult + Q_solar_ventral * vmult
    Q_longwave_in  = Q_rad_in_dorsal * dmult + Q_rad_in_ventral * vmult
    Q_evaporation  = Q_evap_skin_dorsal * dmult + Q_evap_skin_ventral * vmult + 
        Q_evap_insulation_dorsal * dmult + Q_evap_insulation_ventral * vmult + Q_resp
    Q_longwave_out = Q_rad_out_dorsal * dmult + Q_rad_out_ventral * vmult
    Q_convection  = Q_convection_dorsal * dmult + Q_convection_ventral * vmult
    Q_conduction  = Q_conduction_dorsal * dmult + Q_conduction_ventral * vmult
    @show Q_rad_out_dorsal, Q_rad_in_dorsal, Q_rad_out_ventral, Q_rad_in_ventral, Q_longwave_in, Q_longwave_out

    insulation_out = insulation_properties(; insulation=insulation_pars, 
        insulation_temperature = T_insulation * 0.7 + T_skin * 0.3, ventral_fraction, 
        insulation_depth_dorsal, insulation_depth_ventral)
    k_insulation_effective = insulation_out.effective_conductivities[1]
    k_insulation_compressed = insulation_out.insulation_conductivity_compressed
    T_skin = T_skin_dorsal * dmult + T_skin_ventral * vmult
    T_insulation = T_insulation_dorsal * dmult + T_insulation_ventral * vmult

    thermoregulation = (;
        T_core=T_core_last, T_skin, T_insulation, T_lung, T_skin_dorsal, T_skin_ventral, T_insulation_dorsal, T_insulation_ventral,
        shape_b=shape_b_last, pant = pant_last, skin_wetness = skin_wetness_last, k_flesh = k_flesh_last, k_insulation_effective,
        k_insulation_dorsal, k_insulation_ventral, k_insulation_compressed, insulation_depth_dorsal, insulation_depth_ventral, q10
    )

    morphology = (;
        area_total, area_skin, area_evaporation, area_convection, area_conduction, area_silhouette,
        F_sky, F_ground, volume, flesh_volume, 
        characteristic_dimension = geometry_pars.geometry.characteristic_dimension, 
        fat_mass, geometry_pars.geometry.length...
    )

    energy_fluxes = (;
        Q_solar, Q_longwave_in, Q_gen, Q_evaporation, Q_longwave_out, Q_convection, Q_conduction, balance,
        ntry = max(ntry_dorsal, ntry_ventral), success = all((success_dorsal, success_ventral))
    )

    mass_fluxes = (; V_air, V_O2_STP, m_evap, m_resp, m_sweat, J_H2O_in, J_H2O_out, J_O2_in, J_O2_out,
        J_CO2_in, J_CO2_out, J_N2_in, J_N2_out, J_air_in, J_air_out
    )
    return(; thermoregulation, morphology, energy_fluxes, mass_fluxes)
end

