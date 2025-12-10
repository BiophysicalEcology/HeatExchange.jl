function simulsol(; T_skin, T_insulation, simulsol_tolerance, geometry_pars, insulation_pars,
    insulation_out, geom_vars, env_vars, traits)

    insulation_test = u"m"(insulation_out.insulation_test)
    success = true

    if insulation_test > 0.0u"m"
        return (; T_insulation, T_skin, Q_convection, Q_conduction, Q_gen_net, Q_evap_skin, 
            Q_longwave, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation,
            Q_rad_ground, Q_evap_insulation, k_insulation, success, ntry) = solve_with_insulation!(T_skin, T_insulation,
            geometry_pars, insulation_pars, insulation_out, geom_vars, env_vars, traits,
            simulsol_tolerance)
    else
        return (; T_insulation, T_skin, Q_convection, Q_conduction, Q_gen_net, Q_evap_skin,
            Q_longwave, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation,
            Q_rad_ground, Q_evap_insulation, k_insulation, success, ntry) = solve_without_insulation!(T_skin, T_insulation,
            geometry_pars, insulation_pars, insulation_out, geom_vars, env_vars, traits,
            simulsol_tolerance)        
    end
end

function solve_without_insulation!(T_skin, T_insulation,
    geometry_pars, insulation_pars, insulation_out, geom_vars, env_vars, traits,
    simulsol_tolerance
)
    (; side, cd, conduction_fraction, longwave_depth_fraction) = geom_vars
    (; fluid, T_air, T_bush, T_vegetation, T_ground, T_sky, T_substrate, rh, wind_speed, P_atmos, 
        F_sky, F_ground, F_bush, F_vegetation, Q_solar, fO2, fCO2, fN2, convection_enhancement) = env_vars
    (; T_core, k_flesh, k_fat, ϵ_body, skin_wetness, bare_skin_fraction,
        eye_fraction) = traits
    tolerance = simulsol_tolerance
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    ntry = 0

    area_total = total_area(geometry_pars)
    area_evaporation = area_total
    area_convection = area_total * (1 - conduction_fraction)
    volume = geometry_pars.geometry.volume
    r_skin = skin_radius(geometry_pars)
    Q_evap_insulation = 0.0u"W"
    Q_conduction = 0.0u"W"

    while true
        ntry += 1
        for i in 1:20

            # Evaporative heat loss
            (; hc, hd, hd_free) = convection(; 
                body = geometry_pars,
                area = area_convection,
                T_air,
                T_surface = T_skin,
                wind_speed,
                P_atmos,
                fluid,
                fO2,
                fCO2,
                fN2,
                convection_enhancement,
                )
            Q_evap_skin = evaporation(; 
                T_surface = T_skin,
                wetness = skin_wetness,
                area = area_evaporation,
                hd,
                hd_free,
                eye_fraction,
                bare_fraction=bare_skin_fraction,
                T_air, 
                rh,
                P_atmos,
                fO2,
                fCO2,
                fN2).Q_evap

            # Q_rad variables for radiant exchange
            Q_rad1 = area_convection * (F_sky * 4.0 * ϵ_body * σ * ((T_skin + T_sky) / 2)^3)
            Q_rad2 = area_convection * (F_bush * 4.0 * ϵ_body * σ * ((T_skin + T_bush) / 2)^3)
            Q_rad3 = area_convection * (F_vegetation * 4.0 * ϵ_body * σ * ((T_skin + T_vegetation) / 2)^3)
            Q_rad4 = area_convection * (F_ground * 4.0 * ϵ_body * σ * ((T_skin + T_ground) / 2)^3)
            T_skin1 = ((4.0 * k_flesh * volume) / (r_skin^2) * T_core) - Q_evap_skin +
                      hc * area_convection * T_air + Q_solar
            T_skin2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground
            T_skin3 = ((4.0 * k_flesh * volume) / (r_skin^2)) + hc * area_convection + Q_rad1 +
                      Q_rad2 + Q_rad3 + Q_rad4

            T_skin_calc = (T_skin1 + T_skin2) / T_skin3

            Q_rad_sky = Q_rad1 * (T_skin_calc - T_sky)
            Q_rad_bush = Q_rad2 * (T_skin_calc - T_bush)
            Q_rad_vegetation = Q_rad3 * (T_skin_calc - T_vegetation)
            Q_rad_ground = Q_rad4 * (T_skin_calc - T_ground)

            Q_longwave = Q_rad_sky + Q_rad_bush + Q_rad_vegetation + Q_rad_ground
            Q_convection = hc * area_convection * (T_skin_calc - T_air)
            #Q_env = Q_longwave + Q_convection - Q_solar

            ΔT_skin = abs(T_skin - T_skin_calc)

            if ΔT_skin < tolerance
                #TODO why is this not shape-specific?
                Q_gen_net = (4 * k_flesh * volume / r_skin^2) * (T_core - T_skin_calc)
                success = true
                return (; T_insulation, T_skin=T_skin_calc, Q_convection, Q_conduction, Q_gen_net,
                            Q_evap_skin, Q_longwave, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, 
                            Q_rad_ground, Q_evap_insulation, k_insulation=nothing, tolerance, success, ntry)
            else
                T_skin = T_skin_calc
                T_insulation = T_skin_calc
                ntry += 1
                if ntry == 101
                    if tolerance <= 0.001u"K"
                        tolerance = 0.01u"K"
                        ntry = 0
                    else
                        success = false
                        Q_gen_net = 0.0u"W"
                        return (; T_insulation, T_skin=T_skin_calc, Q_convection, Q_conduction, Q_gen_net,
                            Q_evap_skin, Q_longwave, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, Q_rad_ground,
                            Q_evap_insulation, k_insulation=nothing, tolerance, success, ntry)
                    end
                end
            end
        end
    end
end

function solve_with_insulation!(T_skin, T_insulation,
    geometry_pars, insulation_pars, insulation_out, geom_vars, env_vars, traits,
    simulsol_tolerance
)
    (; side, cd, ventral_fraction, conduction_fraction, longwave_depth_fraction) = geom_vars
    (; fluid, T_air, T_substrate, T_bush, T_vegetation, T_ground, T_sky, rh, wind_speed, P_atmos, 
        F_sky, F_ground, F_bush, F_vegetation, Q_solar, fO2, fCO2, fN2, convection_enhancement) = env_vars
    (; T_core, k_flesh, k_fat, ϵ_body, skin_wetness, insulation_wetness, bare_skin_fraction,
        eye_fraction, insulation_conductivity) = traits
    
    tolerance = simulsol_tolerance

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    area_evaporation = evaporation_area(geometry_pars)
    area_total = total_area(geometry_pars)
    area_convection = area_total * (1 - conduction_fraction)
    insulation_test = insulation_out.insulation_test

    ntry = 0
    solct = 0
    solution_procedure = 1
    success = true
    Q_gen_net = 0.0u"W"

    while ntry < 20
        ntry += 1
        for i in 1:20
            # Evaporative heat loss
            # first from the skin
            (; hc, hd, hd_free) = convection(;
                body=geometry_pars,
                area=area_convection,
                T_air,
                T_surface=T_insulation,
                wind_speed,
                P_atmos,
                fluid,
                fO2,
                fCO2,
                fN2,
                convection_enhancement)
                
            Q_evap_skin = evaporation(;
                T_surface = T_skin,
                wetness = skin_wetness,
                area=area_evaporation,
                hd,
                hd_free,
                eye_fraction,
                bare_fraction = bare_skin_fraction,
                T_air,
                rh,
                P_atmos, 
                fO2,
                fCO2,
                fN2).Q_evap
            # second from insulation
            if insulation_wetness > 0 && insulation_test > 0.0u"m"
                Q_evap_insulation = evaporation(; 
                    T_surface = T_insulation,
                    wetness = insulation_wetness,
                    area = area_convection,
                    hd,
                    hd_free = hd,
                    eye_fraction = 0.0,
                    bare_fraction = 1.0, #
                    T_air, 
                    rh,
                    P_atmos,
                    fO2, 
                    fCO2,
                    fN2).Q_evap
            else
                Q_evap_insulation = 0.0u"W"
            end
            
            # Radiation properties
            (; effective_conductivities, absorption_coefficients) = insulation_out =
                insulation_properties(;
                    insulation = insulation_pars,
                    insulation_temperature = T_insulation * 0.7 + T_skin * 0.3,
                    ventral_fraction,
                )

            absorption_coefficient = absorption_coefficients[side + 1]
            k_eff = effective_conductivities[side + 1]

            # Effective insulation conductivity
            if !isnothing(insulation_conductivity)
                k_insulation = insulation_conductivity
            else
                T_rad_approx = T_skin * (1 - longwave_depth_fraction) +
                               T_insulation * longwave_depth_fraction
                k_rad = (16 * σ * T_rad_approx^3) / (3 * absorption_coefficient)
                k_insulation = k_eff + k_rad
            end

            (; T_radiant, T_ins_compressed, cd1, cd2, cd3, dv1, dv2, dv3, dv4) =
                radiant_temperature(; 
                    body = geometry_pars,
                    insulation = insulation_out,
                    insulation_pars,
                    Q_evap = Q_evap_skin,
                    T_core,
                    T_skin,
                    T_substrate,
                    T_insulation, 
                    k_flesh, 
                    k_fat,
                    k_insulation,
                    cd,
                    longwave_depth_fraction, 
                    conduction_fraction)

            # Radiative heat fluxes
            Q_rad1 = area_convection * F_sky * 4 * ϵ_body * σ * ((T_radiant + T_sky) / 2)^3
            Q_rad2 = area_convection * F_bush * 4 * ϵ_body * σ * ((T_radiant + T_bush) / 2)^3
            Q_rad3 = area_convection * F_vegetation * 4 * ϵ_body * σ * ((T_radiant + T_vegetation) / 2)^3
            Q_rad4 = area_convection * F_ground * 4 * ϵ_body * σ * ((T_radiant + T_ground) / 2)^3

            if conduction_fraction < 1
                # These calculations are for when there is less than 100% conduction.
                # The term Q_evap_insulation is included for heat lost due to evaporation from
                # the insulation surface
                (; T_insulation_calc, T_radiant2) =
                    insulation_radiant_temperature(;
                        body = geometry_pars,
                        insulation=insulation_out,
                        insulation_pars,
                        T_core,
                        T_ins_compressed, 
                        T_air,
                        T_sky,
                        T_ground, 
                        T_vegetation, 
                        T_bush,
                        T_substrate, 
                        area_convection,
                        hc,
                        cd,
                        k_insulation,
                        Q_solar,
                        Q_evap_insulation,
                        Q_rad1, Q_rad2, Q_rad3, Q_rad4,
                        cd1, cd2, cd3, dv1, dv2, dv3, dv4,
                        longwave_depth_fraction,
                        conduction_fraction)

                Q_rad_sky = Q_rad1 * (T_radiant2 - T_sky)
                Q_rad_bush = Q_rad2 * (T_radiant2 - T_bush)
                Q_rad_vegetation = Q_rad3 * (T_radiant2 - T_vegetation)
                Q_rad_ground = Q_rad4 * (T_radiant2 - T_ground)
                Q_longwave = Q_rad_sky + Q_rad_bush + Q_rad_vegetation + Q_rad_ground
                Q_convection = hc * area_convection * (T_insulation_calc - T_air)
                Q_conduction = cd * (T_ins_compressed - T_substrate)
            else
                (; T_ins_compressed) =
                    compressed_radiant_temperature(;
                        body = geometry_pars,
                        insulation = insulation_out,
                        insulation_pars,
                        k_flesh,
                        k_fat,
                        T_core,
                        T_substrate,
                        cd,
                        )
                Q_rad_sky = 0.0u"W"
                Q_rad_bush = 0.0u"W"
                Q_rad_vegetation = 0.0u"W"
                Q_rad_ground = 0.0u"W"        
                Q_longwave = 0.0u"W"
                Q_convection = 0.0u"W"
                Q_evap_insulation = 0.0u"W"
                Q_solar = 0.0u"W"
                Q_conduction = cd * (T_ins_compressed - T_substrate)
                T_insulation_calc = T_ins_compressed
            end
            Q_env = Q_longwave + Q_convection + Q_conduction + Q_evap_insulation - Q_solar
            T_skin_mean, T_skin_calc1 =
                mean_skin_temperature(; 
                    body = geometry_pars,
                    insulation=insulation_out,
                    insulation_pars,
                    Q_env,
                    Q_evap_skin,
                    k_flesh,
                    k_fat,
                    T_core,
                    T_insulation_calc,
                    T_ins_compressed,
                    cd1,
                    cd2,
                    cd3,
                    conduction_fraction,
                    )

            ΔT_insulation = abs(T_insulation - T_insulation_calc)
            ΔT_skin = abs(T_skin - T_skin_mean)

            # first convergence test (ΔT_insulation)
            if ΔT_insulation < tolerance
                # Next check T_skin convergence
                if ΔT_skin < tolerance
                    Q_gen_net = net_metabolic_heat(; body = geometry_pars, T_core, T_skin, k_flesh, k_fat)
                    success = true
                    return return (; T_insulation, T_skin=T_skin_mean, Q_convection, Q_conduction,
                        Q_gen_net, Q_evap_skin, Q_longwave, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation,
                        Q_rad_ground, Q_evap_insulation, k_insulation, tolerance, success, ntry)
                else
                    # Not converged, restart iteration
                    if ntry < 20
                        T_skin = T_skin_calc1
                        continue
                    else
                        success_ref[] = 0
                        Q_gen_net = compute_qgennet(
                            ipt, tc, tskcalcav,
                            rflesh, rskin, ak1, ak2, vol,
                            ssqg, bs, bg
                        )
                        return (; T_insulation, T_skin_mean, Q_convection, Q_conduction, Q_gen_net,
                            Q_evap_skin, Q_longwave, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, Q_rad_ground,
                            Q_evap_insulation, k_insulation, tolerance, success, ntry)
                    end
                end

            else
                # No ΔT_insulation convergence → update T_insulation
                T_insulation = update_T_insulation!(T_insulation, T_insulation_calc, ΔT_insulation, solution_procedure)
            end
            # update T_skin
            T_skin = T_skin_mean
            solct += 1

            # fallback if stuck
            if solct ≥ 100
                if solution_procedure != 3
                    solct = 0
                    solution_procedure += 1
                else
                    # Didn't converge → relax tolerance or fail
                    if tolerance <= 0.001u"K"
                        tolerance = 0.01u"K"
                        solct = 0
                        solution_procedure = 1
                    else
                        success = false
                        Q_gen_net = 0.0u"W"
                        return (; T_insulation, T_skin, Q_convection, Q_conduction, Q_gen_net,
                            Q_evap_skin, Q_longwave, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, Q_rad_ground,
                            Q_evap_insulation, k_insulation, tolerance, success, ntry)
                    end
                end
            end
        end
    end
    return (; T_insulation, T_skin, Q_convection, Q_conduction, Q_gen_net,
        Q_evap_skin, Q_longwave, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, Q_rad_ground,
        Q_evap_insulation, k_insulation, tolerance, success, ntry)
end

function update_T_insulation!(T_insulation, T_insulation_calc, ΔT_insulation, solution_procedure)
    if solution_procedure == 1
        # first solution procedure: set T_insulation guess to the calculated T_insulation
        T_insulation = T_insulation_calc

    else
        if solution_procedure == 2
            # second solution procedure: set T_insulation to the average of previous and calculated
            T_insulation = (T_insulation_calc + T_insulation) / 2

        else
            # final solution procedure: incrementally adjust T_insulation
            if (T_insulation - T_insulation_calc) < 0.0u"K"
                # T_insulation < T_insulation_calc → increase T_insulation
                if ΔT_insulation > 3.5u"K"
                    T_insulation = T_insulation + 0.5u"K"
                end
                if (ΔT_insulation > 1.0u"K") && (ΔT_insulation < 3.5u"K")
                    T_insulation = T_insulation + 0.05u"K"
                end
                if (ΔT_insulation > 0.1u"K") && (ΔT_insulation < 1.0u"K")
                    T_insulation = T_insulation + 0.05u"K"
                end
                if (ΔT_insulation > 0.01u"K") && (ΔT_insulation < 0.1u"K")
                    T_insulation = T_insulation + 0.005u"K"
                end
                if (ΔT_insulation > 0.0u"K") && (ΔT_insulation < 0.01u"K")
                    T_insulation = T_insulation + 0.0001u"K"
                end
                if (ΔT_insulation > 0.0u"K") && (ΔT_insulation < 0.001u"K")
                    T_insulation = T_insulation + 0.00001u"K"
                end

            else
                # T_insulation > T_insulation_calc → decrease T_insulation
                if ΔT_insulation > 3.5u"K"
                    T_insulation = T_insulation - 0.5u"K"
                end
                if (ΔT_insulation > 1.0u"K") && (ΔT_insulation < 3.5u"K")
                    T_insulation = T_insulation - 0.05u"K"
                end
                if (ΔT_insulation > 0.1u"K") && (ΔT_insulation < 1.0u"K")
                    T_insulation = T_insulation - 0.05u"K"
                end
                if (ΔT_insulation > 0.01u"K") && (ΔT_insulation < 0.1u"K")
                    T_insulation = T_insulation - 0.005u"K"
                end
                if (ΔT_insulation > 0.001u"K") && (ΔT_insulation < 0.01u"K")
                    T_insulation = T_insulation - 0.0001u"K"
                end
                if (ΔT_insulation > 0.0u"K") && (ΔT_insulation < 0.001u"K")
                    T_insulation = T_insulation - 0.00001u"K"
                end
            end
        end
    end
    return T_insulation
end