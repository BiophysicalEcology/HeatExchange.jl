"""
    simulsol(; body, insulation_pars, insulation_out, geom_vars, env_vars, traits, simulsol_tolerance, skin_temperature, insulation_temperature)

Simultaneously solve for skin and insulation surface temperatures.

Iteratively finds temperatures that satisfy heat balance equations for an endotherm,
accounting for convection, radiation, evaporation, and conduction through insulation.

# Keywords
- `body::AbstractBody`: Body geometry
- `insulation_pars::InsulationParameters`: Insulation parameters
- `insulation_out::InsulationProperties`: Computed insulation properties
- `geom_vars::GeometryVariables`: Geometric variables (side, conduction_fraction, etc.)
- `env_vars::NamedTuple`: Environmental variables (temperatures, wind, humidity, etc.)
- `traits::NamedTuple`: Organism traits (core_temperature, conductivities, emissivity, etc.)
- `simulsol_tolerance`: Convergence tolerance for temperature iteration
- `skin_temperature`: Initial guess for skin temperature
- `insulation_temperature`: Initial guess for insulation surface temperature

# Returns
NamedTuple with:
- `insulation_temperature`: Converged insulation surface temperature
- `skin_temperature`: Converged mean skin temperature
- `fluxes::HeatFluxes`: Heat flux components
- `k_insulation`: Effective insulation conductivity
- `tolerance`: Final tolerance used
- `success`: Whether convergence was achieved
- `ntry`: Number of iterations
"""
function simulsol(;
    body::AbstractBody,
    insulation_pars::InsulationParameters,
    insulation_out::InsulationProperties,
    geom_vars::GeometryVariables,
    env_vars::NamedTuple,
    traits::NamedTuple,
    simulsol_tolerance,
    skin_temperature,
    insulation_temperature,
)
    insulation_test = u"m"(insulation_out.insulation_test)
    success = true

    if insulation_test > 0.0u"m"
        return solve_with_insulation!(
            body,
            insulation_pars,
            insulation_out,
            geom_vars,
            env_vars,
            traits,
            simulsol_tolerance,
            skin_temperature,
            insulation_temperature,
        )
    else
        return solve_without_insulation!(
            body,
            geom_vars,
            env_vars,
            traits,
            simulsol_tolerance,
            skin_temperature,
            insulation_temperature,
        )
    end
end

function solve_without_insulation!(
    body::AbstractBody, geom_vars::GeometryVariables, env_vars::NamedTuple, traits::NamedTuple, simulsol_tolerance, skin_temperature, insulation_temperature
)
    (;
        temperature,
        view_factors,
        atmos,
        fluid,
        solar_flux,
        gas_fractions,
        convection_enhancement,
    ) = env_vars
    air_temperature = temperature.air
    sky_temperature = temperature.sky
    ground_temperature = temperature.ground
    vegetation_temperature = temperature.vegetation
    bush_temperature = temperature.bush
    F_sky = view_factors.sky
    F_ground = view_factors.ground
    F_bush = view_factors.bush
    F_vegetation = view_factors.vegetation
    (; relative_humidity, wind_speed, atmospheric_pressure) = atmos
    (; core_temperature, k_flesh, ϵ_body, skin_wetness, bare_skin_fraction, eye_fraction) = traits
    tolerance = simulsol_tolerance
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    ntry = 0
    volume = flesh_volume(body)
    area_total = total_area(body)
    area_evaporation = area_total
    area_convection = area_total #* (1 - conduction_fraction)
    r_skin = skin_radius(body)
    insulation_evaporation_flux = 0.0u"W"
    conduction_flux = 0.0u"W"

    while true
        ntry += 1
        for i in 1:20

            # Evaporative heat loss
            (; hc, hd, hd_free) = convection(;
                body=body,
                area=area_convection,
                air_temperature,
                surface_temperature=skin_temperature,
                wind_speed,
                atmospheric_pressure,
                fluid,
                gas_fractions,
                convection_enhancement,
            )
            evap_pars_local = EvaporationParameters(;
                skin_wetness,
                eye_fraction,
                bare_skin_fraction,
            )
            transfer_local = TransferCoefficients(; heat=hc, mass=hd, mass_free=hd_free)
            atmos_local = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
            skin_evaporation_flux = evaporation(
                evap_pars_local,
                transfer_local,
                atmos_local,
                area_evaporation,
                skin_temperature,
                air_temperature;
                gas_fractions,
            ).evaporation_flux

            # Q_rad variables for radiant exchange
            sky_radiation_coeff = area_convection * (F_sky * 4.0 * ϵ_body * σ * ((skin_temperature + sky_temperature) / 2)^3)
            bush_radiation_coeff =
                area_convection * (F_bush * 4.0 * ϵ_body * σ * ((skin_temperature + bush_temperature) / 2)^3)
            vegetation_radiation_coeff =
                area_convection *
                (F_vegetation * 4.0 * ϵ_body * σ * ((skin_temperature + vegetation_temperature) / 2)^3)
            ground_radiation_coeff =
                area_convection *
                (F_ground * 4.0 * ϵ_body * σ * ((skin_temperature + ground_temperature) / 2)^3)
            skin_temperature1 =
                ((4.0 * k_flesh * volume) / (r_skin^2) * core_temperature) - skin_evaporation_flux +
                hc * area_convection * air_temperature +
                solar_flux
            skin_temperature2 =
                sky_radiation_coeff * sky_temperature + bush_radiation_coeff * bush_temperature + vegetation_radiation_coeff * vegetation_temperature + ground_radiation_coeff * ground_temperature
            skin_temperature3 =
                ((4.0 * k_flesh * volume) / (r_skin^2)) +
                hc * area_convection +
                sky_radiation_coeff +
                bush_radiation_coeff +
                vegetation_radiation_coeff +
                ground_radiation_coeff

            skin_temperature_calc = (skin_temperature1 + skin_temperature2) / skin_temperature3

            sky_radiation_flux = sky_radiation_coeff * (skin_temperature_calc - sky_temperature)
            bush_radiation_flux = bush_radiation_coeff * (skin_temperature_calc - bush_temperature)
            vegetation_radiation_flux = vegetation_radiation_coeff * (skin_temperature_calc - vegetation_temperature)
            ground_radiation_flux = ground_radiation_coeff * (skin_temperature_calc - ground_temperature)

            longwave_flux = sky_radiation_flux + bush_radiation_flux + vegetation_radiation_flux + ground_radiation_flux
            convection_flux = hc * area_convection * (skin_temperature_calc - air_temperature)

            # Build fluxes (net_generated_flux updated on success)
            fluxes = HeatFluxes(
                convection_flux,
                conduction_flux,
                0.0u"W",  # net_generated_flux placeholder
                skin_evaporation_flux,
                insulation_evaporation_flux,
                longwave_flux,
                solar_flux,
                sky_radiation_flux,
                bush_radiation_flux,
                vegetation_radiation_flux,
                ground_radiation_flux,
            )

            Δskin_temperature = abs(skin_temperature - skin_temperature_calc)

            if Δskin_temperature < tolerance
                #TODO why is this not shape-specific?
                net_generated_flux = (4 * k_flesh * volume / r_skin^2) * (core_temperature - skin_temperature_calc)
                fluxes = setproperties(fluxes; net_generated_flux)
                return (;
                    insulation_temperature,
                    skin_temperature=skin_temperature_calc,
                    fluxes,
                    k_insulation=nothing,
                    tolerance,
                    success=true,
                    ntry,
                )
            else
                skin_temperature = skin_temperature_calc
                insulation_temperature = skin_temperature_calc
                ntry += 1
                if ntry == 101
                    if tolerance <= 0.001u"K"
                        tolerance = 0.01u"K"
                        ntry = 0
                    else
                        return (;
                            insulation_temperature,
                            skin_temperature=skin_temperature_calc,
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
    body::AbstractBody,
    ins::InsulationParameters,
    insulation_out::InsulationProperties,
    geom_vars::GeometryVariables,
    env_vars::NamedTuple,
    traits::NamedTuple,
    simulsol_tolerance,
    skin_temperature,
    insulation_temperature,
)
    (; side, substrate_conductance, ventral_fraction, conduction_fraction, longwave_depth_fraction) = geom_vars
    cd = substrate_conductance  # short name for math
    (;
        temperature,
        view_factors,
        atmos,
        fluid,
        solar_flux,
        gas_fractions,
        convection_enhancement,
    ) = env_vars
    air_temperature = temperature.air
    sky_temperature = temperature.sky
    ground_temperature = temperature.ground
    vegetation_temperature = temperature.vegetation
    bush_temperature = temperature.bush
    substrate_temperature = temperature.substrate
    F_sky = view_factors.sky
    F_ground = view_factors.ground
    F_bush = view_factors.bush
    F_vegetation = view_factors.vegetation
    (; relative_humidity, wind_speed, atmospheric_pressure) = atmos
    (;
        core_temperature,
        k_flesh,
        k_fat,
        ϵ_body,
        skin_wetness,
        insulation_wetness,
        bare_skin_fraction,
        eye_fraction,
    ) = traits

    tolerance = simulsol_tolerance

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    area_evaporation = evaporation_area(body)
    area_total = total_area(body)
    area_convection = area_total * (1 - conduction_fraction)
    insulation_test = insulation_out.insulation_test

    ntry = 0
    solct = 0
    solution_procedure = 1
    success = true
    net_generated_flux = 0.0u"W"

    while ntry < 20
        ntry += 1
        for i in 1:20
            # Evaporative heat loss
            # first from the skin
            (; hc, hd, hd_free) = convection(;
                body=body,
                area=area_convection,
                air_temperature,
                surface_temperature=insulation_temperature,
                wind_speed,
                atmospheric_pressure,
                fluid,
                gas_fractions,
                convection_enhancement,
            )
            evap_pars_skin = EvaporationParameters(;
                skin_wetness,
                eye_fraction,
                bare_skin_fraction,
            )
            transfer_skin = TransferCoefficients(; heat=hc, mass=hd, mass_free=hd_free)
            atmos_skin = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
            skin_evaporation_flux = evaporation(
                evap_pars_skin,
                transfer_skin,
                atmos_skin,
                area_evaporation,
                skin_temperature,
                air_temperature;
                gas_fractions,
            ).evaporation_flux
            # second from insulation
            if insulation_wetness > 0 && insulation_test > 0.0u"m"
                evap_pars_ins = EvaporationParameters(;
                    skin_wetness=insulation_wetness,
                    eye_fraction=0.0,
                    bare_skin_fraction=1.0,
                )
                transfer_ins = TransferCoefficients(; heat=hc, mass=hd, mass_free=hd)
                insulation_evaporation_flux = evaporation(
                    evap_pars_ins,
                    transfer_ins,
                    atmos_skin,
                    area_convection,
                    insulation_temperature,
                    air_temperature;
                    gas_fractions,
                ).evaporation_flux
            else
                insulation_evaporation_flux = 0.0u"W"
            end
            # Recompute insulation thermal properties for current temperatures
            insulation_temp = insulation_temperature * 0.7 + skin_temperature * 0.3
            air_k = dry_air_properties(insulation_temp).thermal_conductivity
            side_depth = getproperty(insulation_out.fibres, side).depth
            side_fibres = setproperties(getproperty(ins, side); depth=side_depth)
            side_thermal = insulation_thermal_conductivity(side_fibres, air_k)

            # Update conductivities for the current side
            conductivities = if side == :dorsal
                setproperties(insulation_out.conductivities; dorsal=side_thermal.effective_conductivity)
            else
                setproperties(insulation_out.conductivities; ventral=side_thermal.effective_conductivity)
            end
            absorption_coefficient = getproperty(insulation_out.absorption_coefficients, side)
            k_eff = getproperty(conductivities, side)
            # update thermally sensitive insulation parameters for current skin/insulation temperature
            insulation_out = setproperties(insulation_out;
                conductivity_compressed=side_thermal.effective_conductivity,
                conductivities,
            )
            # Effective insulation conductivity
            approx_radiant_temperature =
                skin_temperature * (1 - longwave_depth_fraction) +
                insulation_temperature * longwave_depth_fraction
            k_rad = (16 * σ * approx_radiant_temperature^3) / (3 * absorption_coefficient)
            k_insulation = k_eff + k_rad
            ks = ThermalConductivities(k_flesh, k_fat, k_insulation)
            org_temps = OrganismTemperatures(core_temperature, skin_temperature, insulation_temperature)
            (; radiant_temperature, compressed_insulation_temperature, cds, dvs) = radiant_temperature(;
                body=body,
                insulation=insulation_out,
                insulation_pars=ins,
                org_temps,
                ks,
                side,
                cd,
                longwave_depth_fraction,
                conduction_fraction,
                evaporation_flux=skin_evaporation_flux,
                substrate_temperature,
            )
            # Radiative heat fluxes
            sky_radiation_coeff = area_convection * F_sky * 4 * ϵ_body * σ * ((radiant_temperature + sky_temperature) / 2)^3
            bush_radiation_coeff =
                area_convection * F_bush * 4 * ϵ_body * σ * ((radiant_temperature + bush_temperature) / 2)^3
            vegetation_radiation_coeff =
                area_convection *
                F_vegetation *
                4 *
                ϵ_body *
                σ *
                ((radiant_temperature + vegetation_temperature) / 2)^3
            ground_radiation_coeff =
                area_convection * F_ground * 4 * ϵ_body * σ * ((radiant_temperature + ground_temperature) / 2)^3
            radiation_coeffs = RadiationCoeffs(sky_radiation_coeff, bush_radiation_coeff, vegetation_radiation_coeff, ground_radiation_coeff)
            env_temps = EnvironmentTemperatures(
                air_temperature, sky_temperature, ground_temperature, vegetation_temperature, bush_temperature, substrate_temperature
            )
            if conduction_fraction < 1
                # These calculations are for when there is less than 100% conduction.
                # The term insulation_evaporation_flux is included for heat lost due to evaporation from
                # the insulation surface
                (; insulation_temperature_calc, radiant_temperature2) = insulation_radiant_temperature(;
                    body=body,
                    insulation=insulation_out,
                    insulation_pars=ins,
                    env_temps,
                    radiation_coeffs,
                    cds,
                    dvs,
                    side,
                    area_convection,
                    hc,
                    cd,
                    k_insulation,
                    longwave_depth_fraction,
                    conduction_fraction,
                    solar_flux,
                    insulation_evaporation_flux,
                    core_temperature,
                    compressed_insulation_temperature,
                )

                sky_radiation_flux = sky_radiation_coeff * (radiant_temperature2 - sky_temperature)
                bush_radiation_flux = bush_radiation_coeff * (radiant_temperature2 - bush_temperature)
                vegetation_radiation_flux = vegetation_radiation_coeff * (radiant_temperature2 - vegetation_temperature)
                ground_radiation_flux = ground_radiation_coeff * (radiant_temperature2 - ground_temperature)
                longwave_flux = sky_radiation_flux + bush_radiation_flux + vegetation_radiation_flux + ground_radiation_flux
                convection_flux = hc * area_convection * (insulation_temperature_calc - air_temperature)
                conduction_flux = u"W"(cd * (compressed_insulation_temperature - substrate_temperature))
            else
                (; compressed_insulation_temperature) = compressed_radiant_temperature(;
                    body=body,
                    insulation=insulation_out,
                    insulation_pars=ins,
                    ks,
                    side,
                    cd,
                    core_temperature,
                    substrate_temperature,
                )
                sky_radiation_flux = 0.0u"W"
                bush_radiation_flux = 0.0u"W"
                vegetation_radiation_flux = 0.0u"W"
                ground_radiation_flux = 0.0u"W"
                longwave_flux = 0.0u"W"
                convection_flux = 0.0u"W"
                insulation_evaporation_flux = 0.0u"W"
                solar_flux = 0.0u"W"
                conduction_flux = cd * (compressed_insulation_temperature - substrate_temperature)
                insulation_temperature_calc = compressed_insulation_temperature
            end
            environment_flux = longwave_flux + convection_flux + conduction_flux + insulation_evaporation_flux - solar_flux
            skin_temperature_mean, skin_temperature_calc1 = mean_skin_temperature(;
                body=body,
                insulation=insulation_out,
                insulation_pars=ins,
                ks,
                cds,
                conduction_fraction,
                environment_flux,
                skin_evaporation_flux,
                core_temperature,
                insulation_temperature_calc,
                compressed_insulation_temperature,
            )

            # Build fluxes (net_generated_flux updated on success)
            fluxes = HeatFluxes(
                convection_flux,
                conduction_flux,
                net_generated_flux,  # placeholder, updated below
                skin_evaporation_flux,
                insulation_evaporation_flux,
                longwave_flux,
                solar_flux,
                sky_radiation_flux,
                bush_radiation_flux,
                vegetation_radiation_flux,
                ground_radiation_flux,
            )

            Δinsulation_temperature = abs(insulation_temperature - insulation_temperature_calc)
            Δskin_temperature = abs(skin_temperature - skin_temperature_mean)

            # first convergence test (Δinsulation_temperature)
            if Δinsulation_temperature < tolerance
                # Next check skin_temperature convergence
                if Δskin_temperature < tolerance
                    net_generated_flux = net_metabolic_heat(; body=body, ks, core_temperature, skin_temperature)
                    fluxes = setproperties(fluxes; net_generated_flux)
                    return (;
                        insulation_temperature,
                        skin_temperature=skin_temperature_mean,
                        fluxes,
                        k_insulation,
                        tolerance,
                        success=true,
                        ntry,
                    )
                else
                    # Not converged, restart iteration
                    if ntry < 20
                        skin_temperature = skin_temperature_calc1
                        continue
                    else
                        net_generated_flux = net_metabolic_heat(; body=body, ks, core_temperature, skin_temperature)
                        fluxes = setproperties(fluxes; net_generated_flux)
                        return (;
                            insulation_temperature,
                            skin_temperature=skin_temperature_mean,
                            fluxes,
                            k_insulation,
                            tolerance,
                            success=false,
                            ntry,
                        )
                    end
                end

            else
                # No Δinsulation_temperature convergence → update insulation_temperature
                insulation_temperature = update_insulation_temperature!(
                    insulation_temperature, insulation_temperature_calc, Δinsulation_temperature, solution_procedure
                )
            end
            # update skin_temperature
            skin_temperature = skin_temperature_mean
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
                            insulation_temperature,
                            skin_temperature,
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
    return (; insulation_temperature, skin_temperature, fluxes, k_insulation, tolerance, success, ntry)
end

function update_insulation_temperature!(
    insulation_temperature, insulation_temperature_calc, Δinsulation_temperature, solution_procedure
)
    if solution_procedure == 1
        # first solution procedure: set insulation_temperature guess to the calculated insulation_temperature
        insulation_temperature = insulation_temperature_calc

    else
        if solution_procedure == 2
            # second solution procedure: set insulation_temperature to the average of previous and calculated
            insulation_temperature = (insulation_temperature_calc + insulation_temperature) / 2

        else
            # final solution procedure: incrementally adjust insulation_temperature
            if (insulation_temperature - insulation_temperature_calc) < 0.0u"K"
                # insulation_temperature < insulation_temperature_calc → increase insulation_temperature
                if Δinsulation_temperature > 3.5u"K"
                    insulation_temperature = insulation_temperature + 0.5u"K"
                end
                if (Δinsulation_temperature > 1.0u"K") && (Δinsulation_temperature < 3.5u"K")
                    insulation_temperature = insulation_temperature + 0.05u"K"
                end
                if (Δinsulation_temperature > 0.1u"K") && (Δinsulation_temperature < 1.0u"K")
                    insulation_temperature = insulation_temperature + 0.05u"K"
                end
                if (Δinsulation_temperature > 0.01u"K") && (Δinsulation_temperature < 0.1u"K")
                    insulation_temperature = insulation_temperature + 0.005u"K"
                end
                if (Δinsulation_temperature > 0.0u"K") && (Δinsulation_temperature < 0.01u"K")
                    insulation_temperature = insulation_temperature + 0.0001u"K"
                end
                if (Δinsulation_temperature > 0.0u"K") && (Δinsulation_temperature < 0.001u"K")
                    insulation_temperature = insulation_temperature + 0.00001u"K"
                end

            else
                # insulation_temperature > insulation_temperature_calc → decrease insulation_temperature
                if Δinsulation_temperature > 3.5u"K"
                    insulation_temperature = insulation_temperature - 0.5u"K"
                end
                if (Δinsulation_temperature > 1.0u"K") && (Δinsulation_temperature < 3.5u"K")
                    insulation_temperature = insulation_temperature - 0.05u"K"
                end
                if (Δinsulation_temperature > 0.1u"K") && (Δinsulation_temperature < 1.0u"K")
                    insulation_temperature = insulation_temperature - 0.05u"K"
                end
                if (Δinsulation_temperature > 0.01u"K") && (Δinsulation_temperature < 0.1u"K")
                    insulation_temperature = insulation_temperature - 0.005u"K"
                end
                if (Δinsulation_temperature > 0.001u"K") && (Δinsulation_temperature < 0.01u"K")
                    insulation_temperature = insulation_temperature - 0.0001u"K"
                end
                if (Δinsulation_temperature > 0.0u"K") && (Δinsulation_temperature < 0.001u"K")
                    insulation_temperature = insulation_temperature - 0.00001u"K"
                end
            end
        end
    end
    return insulation_temperature
end
