function simulsol(;
    T_skin,
    T_insulation,
    simulsol_tolerance,
    geometry_pars::AbstractBody,
    insulation_pars::InsulationParameters,
    insulation_out::InsulationOutput,
    geom_vars::GeometryVariables,
    env_vars::NamedTuple,
    traits::NamedTuple,
)
    insulation_test = u"m"(insulation_out.insulation_test)
    success = true

    if insulation_test > 0.0u"m"
        return solve_with_insulation!(
            T_skin,
            T_insulation,
            geometry_pars,
            insulation_pars,
            insulation_out,
            geom_vars,
            env_vars,
            traits,
            simulsol_tolerance,
        )
    else
        return solve_without_insulation!(
            T_skin,
            T_insulation,
            geometry_pars,
            geom_vars,
            env_vars,
            traits,
            simulsol_tolerance,
        )
    end
end

function solve_without_insulation!(
    T_skin, T_insulation, geometry_pars::AbstractBody, geom_vars::GeometryVariables, env_vars::NamedTuple, traits::NamedTuple, simulsol_tolerance
)
    #(; conduction_fraction) = geom_vars
    (;
        view_factors,
        atmos,
        fluid,
        T_air,
        T_bush,
        T_vegetation,
        T_ground,
        T_sky,
        Q_solar,
        gasfrac,
        convection_enhancement,
    ) = env_vars
    (; F_sky, F_ground, F_bush, F_vegetation) = view_factors
    (; rh, wind_speed, P_atmos) = atmos
    (; T_core, k_flesh, ϵ_body, skin_wetness, bare_skin_fraction, eye_fraction) = traits
    tolerance = simulsol_tolerance
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    ntry = 0
    volume = flesh_volume(geometry_pars)
    area_total = total_area(geometry_pars)
    area_evaporation = area_total
    area_convection = area_total #* (1 - conduction_fraction)
    r_skin = skin_radius(geometry_pars)
    Q_evap_insulation = 0.0u"W"
    Q_conduction = 0.0u"W"

    while true
        ntry += 1
        for i in 1:20

            # Evaporative heat loss
            (; hc, hd, hd_free) = convection(;
                body=geometry_pars,
                area=area_convection,
                T_air,
                T_surface=T_skin,
                wind_speed,
                P_atmos,
                fluid,
                gasfrac,
                convection_enhancement,
            )
            evap_pars_local = EvaporationParameters(;
                skin_wetness,
                eye_fraction,
                bare_skin_fraction,
            )
            transfer_local = TransferCoefficients(; heat=hc, mass=hd, mass_free=hd_free)
            atmos_local = AtmosphericConditions(rh, wind_speed, P_atmos)
            Q_evap_skin = evaporation(
                T_skin,
                T_air,
                area_evaporation,
                evap_pars_local,
                transfer_local,
                atmos_local;
                gasfrac,
            ).Q_evap

            # Q_rad variables for radiant exchange
            Q_rad1 = area_convection * (F_sky * 4.0 * ϵ_body * σ * ((T_skin + T_sky) / 2)^3)
            Q_rad2 =
                area_convection * (F_bush * 4.0 * ϵ_body * σ * ((T_skin + T_bush) / 2)^3)
            Q_rad3 =
                area_convection *
                (F_vegetation * 4.0 * ϵ_body * σ * ((T_skin + T_vegetation) / 2)^3)
            Q_rad4 =
                area_convection *
                (F_ground * 4.0 * ϵ_body * σ * ((T_skin + T_ground) / 2)^3)
            T_skin1 =
                ((4.0 * k_flesh * volume) / (r_skin^2) * T_core) - Q_evap_skin +
                hc * area_convection * T_air +
                Q_solar
            T_skin2 =
                Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground
            T_skin3 =
                ((4.0 * k_flesh * volume) / (r_skin^2)) +
                hc * area_convection +
                Q_rad1 +
                Q_rad2 +
                Q_rad3 +
                Q_rad4

            T_skin_calc = (T_skin1 + T_skin2) / T_skin3

            Q_rad_sky = Q_rad1 * (T_skin_calc - T_sky)
            Q_rad_bush = Q_rad2 * (T_skin_calc - T_bush)
            Q_rad_vegetation = Q_rad3 * (T_skin_calc - T_vegetation)
            Q_rad_ground = Q_rad4 * (T_skin_calc - T_ground)

            Q_longwave = Q_rad_sky + Q_rad_bush + Q_rad_vegetation + Q_rad_ground
            Q_convection = hc * area_convection * (T_skin_calc - T_air)

            # Build fluxes (Q_gen_net updated on success)
            fluxes = HeatFluxes(
                Q_convection,
                Q_conduction,
                0.0u"W",  # Q_gen_net placeholder
                Q_evap_skin,
                Q_evap_insulation,
                Q_longwave,
                Q_solar,
                Q_rad_sky,
                Q_rad_bush,
                Q_rad_vegetation,
                Q_rad_ground,
            )

            ΔT_skin = abs(T_skin - T_skin_calc)

            if ΔT_skin < tolerance
                #TODO why is this not shape-specific?
                Q_gen_net = (4 * k_flesh * volume / r_skin^2) * (T_core - T_skin_calc)
                fluxes = setproperties(fluxes; Q_gen_net)
                return (;
                    T_insulation,
                    T_skin=T_skin_calc,
                    fluxes,
                    k_insulation=nothing,
                    tolerance,
                    success=true,
                    ntry,
                )
            else
                T_skin = T_skin_calc
                T_insulation = T_skin_calc
                ntry += 1
                if ntry == 101
                    if tolerance <= 0.001u"K"
                        tolerance = 0.01u"K"
                        ntry = 0
                    else
                        return (;
                            T_insulation,
                            T_skin=T_skin_calc,
                            fluxes,
                            k_insulation=nothing,
                            tolerance,
                            success=false,
                            ntry,
                        )
                    end
                end
            end
        end
    end
end

function solve_with_insulation!(
    T_skin,
    T_insulation,
    geometry_pars::AbstractBody,
    ins::InsulationParameters,
    insulation_out::InsulationOutput,
    geom_vars::GeometryVariables,
    env_vars::NamedTuple,
    traits::NamedTuple,
    simulsol_tolerance,
)
    (; side, substrate_conductance, ventral_fraction, conduction_fraction, longwave_depth_fraction) = geom_vars
    cd = substrate_conductance  # short name for math
    (;
        view_factors,
        atmos,
        fluid,
        T_air,
        T_substrate,
        T_bush,
        T_vegetation,
        T_ground,
        T_sky,
        Q_solar,
        gasfrac,
        convection_enhancement,
    ) = env_vars
    (; F_sky, F_ground, F_bush, F_vegetation) = view_factors
    (; rh, wind_speed, P_atmos) = atmos
    (;
        T_core,
        k_flesh,
        k_fat,
        ϵ_body,
        skin_wetness,
        insulation_wetness,
        bare_skin_fraction,
        eye_fraction,
        insulation_conductivity,
    ) = traits

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
                gasfrac,
                convection_enhancement,
            )
            evap_pars_skin = EvaporationParameters(;
                skin_wetness,
                eye_fraction,
                bare_skin_fraction,
            )
            transfer_skin = TransferCoefficients(; heat=hc, mass=hd, mass_free=hd_free)
            atmos_skin = AtmosphericConditions(rh, wind_speed, P_atmos)
            Q_evap_skin = evaporation(
                T_skin,
                T_air,
                area_evaporation,
                evap_pars_skin,
                transfer_skin,
                atmos_skin;
                gasfrac,
            ).Q_evap
            # second from insulation
            if insulation_wetness > 0 && insulation_test > 0.0u"m"
                evap_pars_ins = EvaporationParameters(;
                    skin_wetness=insulation_wetness,
                    eye_fraction=0.0,
                    bare_skin_fraction=1.0,
                )
                transfer_ins = TransferCoefficients(; heat=hc, mass=hd, mass_free=hd)
                Q_evap_insulation = evaporation(
                    T_insulation,
                    T_air,
                    area_convection,
                    evap_pars_ins,
                    transfer_ins,
                    atmos_skin;
                    gasfrac,
                ).Q_evap
            else
                Q_evap_insulation = 0.0u"W"
            end
            if side == 1
                insulation_side = InsulationParameters(;
                    insulation_conductivity_dorsal=ins.insulation_conductivity_dorsal,
                    insulation_conductivity_ventral=ins.insulation_conductivity_dorsal,
                    fibre_diameter_dorsal=ins.fibre_diameter_dorsal,
                    fibre_diameter_ventral=ins.fibre_diameter_dorsal,
                    fibre_length_dorsal=ins.fibre_length_dorsal,
                    fibre_length_ventral=ins.fibre_length_dorsal,
                    insulation_depth_dorsal=insulation_out.insulation_depths[2],
                    insulation_depth_ventral=insulation_out.insulation_depths[2],
                    fibre_density_dorsal=ins.fibre_density_dorsal,
                    fibre_density_ventral=ins.fibre_density_dorsal,
                    insulation_reflectance_dorsal=ins.insulation_reflectance_dorsal,
                    insulation_reflectance_ventral=ins.insulation_reflectance_dorsal,
                    insulation_depth_compressed=insulation_out.insulation_depths[2],
                    fibre_conductivity=ins.fibre_conductivity,
                    longwave_depth_fraction=ins.longwave_depth_fraction,
                )
            else
                insulation_side = InsulationParameters(;
                    insulation_conductivity_dorsal=ins.insulation_conductivity_ventral,
                    insulation_conductivity_ventral=ins.insulation_conductivity_ventral,
                    fibre_diameter_dorsal=ins.fibre_diameter_ventral,
                    fibre_diameter_ventral=ins.fibre_diameter_ventral,
                    fibre_length_dorsal=ins.fibre_length_ventral,
                    fibre_length_ventral=ins.fibre_length_ventral,
                    insulation_depth_dorsal=insulation_out.insulation_depths[3],
                    insulation_depth_ventral=insulation_out.insulation_depths[3],
                    fibre_density_dorsal=ins.fibre_density_ventral,
                    fibre_density_ventral=ins.fibre_density_ventral,
                    insulation_reflectance_dorsal=ins.insulation_reflectance_ventral,
                    insulation_reflectance_ventral=ins.insulation_reflectance_ventral,
                    insulation_depth_compressed=insulation_out.insulation_depths[3], # TODO this should be ins.insulation_depth_compressed here and in Fortran
                    fibre_conductivity=ins.fibre_conductivity,
                    longwave_depth_fraction=ins.longwave_depth_fraction,
                )
            end

            # recompute insulation radiation properties
            (; effective_conductivities, insulation_conductivity_compressed) = insulation_properties(;
                insulation=insulation_side,
                insulation_temperature=T_insulation * 0.7 + T_skin * 0.3,
                ventral_fraction,
            )
            absorption_coefficient = insulation_out.absorption_coefficients[side + 1]
            k_eff = effective_conductivities[side + 1]
            # update thermally sensitive insulation parameters for current skin/insulation temperature
            insulation_out = setproperties(insulation_out;
                insulation_conductivity_compressed,
                effective_conductivities,
            )
            # Effective insulation conductivity
            if !isnothing(insulation_conductivity)
                k_insulation = insulation_conductivity
            else
                T_rad_approx =
                    T_skin * (1 - longwave_depth_fraction) +
                    T_insulation * longwave_depth_fraction
                k_rad = (16 * σ * T_rad_approx^3) / (3 * absorption_coefficient)
                k_insulation = k_eff + k_rad
            end
            ks = ThermalConductivities(k_flesh, k_fat, k_insulation)
            org_temps = OrganismTemperatures(T_core, T_skin, T_insulation)
            (; T_radiant, T_ins_compressed, cds, dvs) = radiant_temperature(;
                body=geometry_pars,
                insulation=insulation_out,
                insulation_pars=ins,
                Q_evap=Q_evap_skin,
                org_temps,
                T_substrate,
                ks,
                cd,
                longwave_depth_fraction,
                conduction_fraction,
                side,
            )
            # Radiative heat fluxes
            Q_rad1 = area_convection * F_sky * 4 * ϵ_body * σ * ((T_radiant + T_sky) / 2)^3
            Q_rad2 =
                area_convection * F_bush * 4 * ϵ_body * σ * ((T_radiant + T_bush) / 2)^3
            Q_rad3 =
                area_convection *
                F_vegetation *
                4 *
                ϵ_body *
                σ *
                ((T_radiant + T_vegetation) / 2)^3
            Q_rad4 =
                area_convection * F_ground * 4 * ϵ_body * σ * ((T_radiant + T_ground) / 2)^3
            Q_rads = RadiationCoeffs(Q_rad1, Q_rad2, Q_rad3, Q_rad4)
            env_temps = EnvironmentTemperatures(
                T_air, T_sky, T_ground, T_vegetation, T_bush, T_substrate
            )
            if conduction_fraction < 1
                # These calculations are for when there is less than 100% conduction.
                # The term Q_evap_insulation is included for heat lost due to evaporation from
                # the insulation surface
                (; T_insulation_calc, T_radiant2) = insulation_radiant_temperature(;
                    body=geometry_pars,
                    insulation=insulation_out,
                    insulation_pars=ins,
                    T_core,
                    T_ins_compressed,
                    env_temps,
                    area_convection,
                    hc,
                    cd,
                    k_insulation,
                    Q_solar,
                    Q_evap_insulation,
                    Q_rads,
                    cds,
                    dvs,
                    longwave_depth_fraction,
                    conduction_fraction,
                    side,
                )

                Q_rad_sky = Q_rad1 * (T_radiant2 - T_sky)
                Q_rad_bush = Q_rad2 * (T_radiant2 - T_bush)
                Q_rad_vegetation = Q_rad3 * (T_radiant2 - T_vegetation)
                Q_rad_ground = Q_rad4 * (T_radiant2 - T_ground)
                Q_longwave = Q_rad_sky + Q_rad_bush + Q_rad_vegetation + Q_rad_ground
                Q_convection = hc * area_convection * (T_insulation_calc - T_air)
                Q_conduction = u"W"(cd * (T_ins_compressed - T_substrate))
            else
                (; T_ins_compressed) = compressed_radiant_temperature(;
                    body=geometry_pars,
                    insulation=insulation_out,
                    insulation_pars=ins,
                    ks,
                    T_core,
                    T_substrate,
                    cd,
                    side,
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
            T_skin_mean, T_skin_calc1 = mean_skin_temperature(;
                body=geometry_pars,
                insulation=insulation_out,
                insulation_pars=ins,
                Q_env,
                Q_evap_skin,
                ks,
                T_core,
                T_insulation_calc,
                T_ins_compressed,
                cds,
                conduction_fraction,
            )

            # Build fluxes (Q_gen_net updated on success)
            fluxes = HeatFluxes(
                Q_convection,
                Q_conduction,
                Q_gen_net,  # placeholder, updated below
                Q_evap_skin,
                Q_evap_insulation,
                Q_longwave,
                Q_solar,
                Q_rad_sky,
                Q_rad_bush,
                Q_rad_vegetation,
                Q_rad_ground,
            )

            ΔT_insulation = abs(T_insulation - T_insulation_calc)
            ΔT_skin = abs(T_skin - T_skin_mean)

            # first convergence test (ΔT_insulation)
            if ΔT_insulation < tolerance
                # Next check T_skin convergence
                if ΔT_skin < tolerance
                    Q_gen_net = net_metabolic_heat(; body=geometry_pars, T_core, T_skin, ks)
                    fluxes = setproperties(fluxes; Q_gen_net)
                    return (;
                        T_insulation,
                        T_skin=T_skin_mean,
                        fluxes,
                        k_insulation,
                        tolerance,
                        success=true,
                        ntry,
                    )
                else
                    # Not converged, restart iteration
                    if ntry < 20
                        T_skin = T_skin_calc1
                        continue
                    else
                        Q_gen_net = net_metabolic_heat(; body=geometry_pars, T_core, T_skin, ks)
                        fluxes = setproperties(fluxes; Q_gen_net)
                        return (;
                            T_insulation,
                            T_skin=T_skin_mean,
                            fluxes,
                            k_insulation,
                            tolerance,
                            success=false,
                            ntry,
                        )
                    end
                end

            else
                # No ΔT_insulation convergence → update T_insulation
                T_insulation = update_T_insulation!(
                    T_insulation, T_insulation_calc, ΔT_insulation, solution_procedure
                )
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
                        return (;
                            T_insulation,
                            T_skin,
                            fluxes,
                            k_insulation,
                            tolerance,
                            success=false,
                            ntry,
                        )
                    end
                end
            end
        end
    end
    return (; T_insulation, T_skin, fluxes, k_insulation, tolerance, success, ntry)
end

function update_T_insulation!(
    T_insulation, T_insulation_calc, ΔT_insulation, solution_procedure
)
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
