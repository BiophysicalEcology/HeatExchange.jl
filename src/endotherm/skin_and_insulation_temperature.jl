"""
    solve_temperatures(; body, insulation_pars, insulation, geom_vars, env_vars, traits, temperature_tolerance, skin_temperature, insulation_temperature)

Simultaneously solve for skin and insulation surface temperatures.

Iteratively finds temperatures that satisfy heat balance equations for an endotherm,
accounting for convection, radiation, evaporation, and conduction through insulation.

# Keywords
- `body::AbstractBody`: Body geometry
- `insulation_pars::InsulationParameters`: Insulation parameters
- `insulation::InsulationProperties`: Computed insulation properties
- `geom_vars::GeometryVariables`: Geometric variables (side, conduction_fraction, etc.)
- `env_vars::NamedTuple`: Environmental variables (temperatures, wind, humidity, etc.)
- `traits::NamedTuple`: Organism traits (core_temperature, conductivities, emissivity, etc.)
- `temperature_tolerance`: Convergence tolerance for temperature iteration
- `skin_temperature`: Initial guess for skin temperature
- `insulation_temperature`: Initial guess for insulation surface temperature

# Returns
NamedTuple with:
- `insulation_temperature`: Converged insulation surface temperature
- `skin_temperature`: Converged mean skin temperature
- `fluxes::HeatFluxes`: Heat flux components
- `insulation_conductivity`: Effective insulation conductivity
- `tolerance`: Final tolerance used
- `success`: Whether convergence was achieved
- `ntry`: Number of iterations
"""
function solve_temperatures(;
    body::AbstractBody,
    insulation_pars::InsulationParameters,
    insulation::InsulationProperties,
    geom_vars::GeometryVariables,
    env_vars::NamedTuple,
    traits::NamedTuple,
    temperature_tolerance,
    skin_temperature,
    insulation_temperature,
)
    insulation_test = u"m"(insulation.insulation_test)
    success = true

    if insulation_test > 0.0u"m"
        return solve_with_insulation!(
            body,
            insulation_pars,
            insulation,
            geom_vars,
            env_vars,
            traits,
            temperature_tolerance,
            skin_temperature,
            insulation_temperature,
        )
    else
        return solve_without_insulation!(
            body,
            geom_vars,
            env_vars,
            traits,
            temperature_tolerance,
            skin_temperature,
            insulation_temperature,
        )
    end
end

function solve_without_insulation!(
    body::AbstractBody, geom_vars::GeometryVariables, env_vars::NamedTuple, traits::NamedTuple, temperature_tolerance, skin_temperature, insulation_temperature
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
    (; relative_humidity, wind_speed, atmospheric_pressure) = atmos
    (; core_temperature, flesh_conductivity, ϵ_body, skin_wetness, bare_skin_fraction, eye_fraction) = traits
    tolerance = temperature_tolerance
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
            conv = convection(;
                body,
                area=area_convection,
                air_temperature,
                surface_temperature=skin_temperature,
                wind_speed,
                atmospheric_pressure,
                fluid,
                gas_fractions,
                convection_enhancement,
            )
            heat_transfer_coefficient = conv.heat.combined
            evap_pars_local = EvaporationParameters(;
                skin_wetness,
                eye_fraction,
                bare_skin_fraction,
            )
            atmos_local = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
            skin_evaporation_flux = evaporation(
                evap_pars_local,
                conv.mass,
                atmos_local,
                area_evaporation,
                skin_temperature,
                air_temperature;
                gas_fractions,
            ).evaporation_flux

            # Q_rad variables for radiant exchange
            sky_radiation_coeff = area_convection * (view_factors.sky * 4.0 * ϵ_body * σ * ((skin_temperature + sky_temperature) / 2)^3)
            bush_radiation_coeff =
                area_convection * (view_factors.bush * 4.0 * ϵ_body * σ * ((skin_temperature + bush_temperature) / 2)^3)
            vegetation_radiation_coeff =
                area_convection *
                (view_factors.vegetation * 4.0 * ϵ_body * σ * ((skin_temperature + vegetation_temperature) / 2)^3)
            ground_radiation_coeff =
                area_convection *
                (view_factors.ground * 4.0 * ϵ_body * σ * ((skin_temperature + ground_temperature) / 2)^3)
            skin_temperature1 =
                ((4.0 * flesh_conductivity * volume) / (r_skin^2) * core_temperature) - skin_evaporation_flux +
                heat_transfer_coefficient * area_convection * air_temperature +
                solar_flux
            skin_temperature2 =
                sky_radiation_coeff * sky_temperature + bush_radiation_coeff * bush_temperature + vegetation_radiation_coeff * vegetation_temperature + ground_radiation_coeff * ground_temperature
            skin_temperature3 =
                ((4.0 * flesh_conductivity * volume) / (r_skin^2)) +
                heat_transfer_coefficient * area_convection +
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
            convection_flux = heat_transfer_coefficient * area_convection * (skin_temperature_calc - air_temperature)

            # Build fluxes (net_generated updated on success)
            fluxes = HeatFluxes(
                convection_flux,
                conduction_flux,
                0.0u"W",  # net_generated placeholder
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
                net_generated = (4 * flesh_conductivity * volume / r_skin^2) * (core_temperature - skin_temperature_calc)
                fluxes = setproperties(fluxes; net_generated)
                return (;
                    insulation_temperature,
                    skin_temperature=skin_temperature_calc,
                    fluxes,
                    insulation_conductivity=nothing,
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
                            insulation_conductivity=nothing,
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
    insulation_pars::InsulationParameters,
    insulation::InsulationProperties,
    geom_vars::GeometryVariables,
    env_vars::NamedTuple,
    traits::NamedTuple,
    temperature_tolerance,
    skin_temperature,
    insulation_temperature,
)
    (; side, substrate_conductance, ventral_fraction, conduction_fraction, longwave_depth_fraction) = geom_vars
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
    (; relative_humidity, wind_speed, atmospheric_pressure) = atmos
    (;
        core_temperature,
        flesh_conductivity,
        fat_conductivity,
        ϵ_body,
        skin_wetness,
        insulation_wetness,
        bare_skin_fraction,
        eye_fraction,
    ) = traits

    tolerance = temperature_tolerance

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    area_evaporation = evaporation_area(body)
    area_total = total_area(body)
    area_convection = area_total * (1 - conduction_fraction)
    insulation_test = insulation.insulation_test

    ntry = 0
    solct = 0
    solution_procedure = 1
    success = true
    net_generated = 0.0u"W"

    while ntry < 20
        ntry += 1
        for i in 1:20
            # Evaporative heat loss
            # first from the skin
            conv = convection(;
                body,
                area=area_convection,
                air_temperature,
                surface_temperature=insulation_temperature,
                wind_speed,
                atmospheric_pressure,
                fluid,
                gas_fractions,
                convection_enhancement,
            )
            heat_transfer_coefficient = conv.heat.combined
            evap_pars_skin = EvaporationParameters(;
                skin_wetness,
                eye_fraction,
                bare_skin_fraction,
            )
            atmos_skin = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
            skin_evaporation_flux = evaporation(
                evap_pars_skin,
                conv.mass,
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
                # Use combined for free (insulation uses combined mass transfer for all)
                mass_ins = TransferCoefficients(conv.mass.combined, conv.mass.combined, conv.mass.forced)
                insulation_evaporation_flux = evaporation(
                    evap_pars_ins,
                    mass_ins,
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
            side_depth = getproperty(insulation.fibres, side).depth
            side_fibres = setproperties(getproperty(insulation_pars, side); depth=side_depth)
            side_thermal = insulation_thermal_conductivity(side_fibres, air_k)

            # Update conductivities for the current side
            side_conductivities = if side == :dorsal
                setproperties(insulation.conductivities; dorsal=side_thermal.effective_conductivity)
            else
                setproperties(insulation.conductivities; ventral=side_thermal.effective_conductivity)
            end
            absorption_coefficient = getproperty(insulation.absorption_coefficients, side)
            effective_conductivity = getproperty(side_conductivities, side)
            # update thermally sensitive insulation parameters for current skin/insulation temperature
            insulation = setproperties(insulation;
                conductivity_compressed=side_thermal.effective_conductivity,
                conductivities=side_conductivities,
            )
            # Effective insulation conductivity
            approx_radiant_temperature =
                skin_temperature * (1 - longwave_depth_fraction) +
                insulation_temperature * longwave_depth_fraction
            radiative_conductivity = (16 * σ * approx_radiant_temperature^3) / (3 * absorption_coefficient)
            insulation_conductivity = effective_conductivity + radiative_conductivity
            conductivities = ThermalConductivities(flesh_conductivity, fat_conductivity, insulation_conductivity)
            org_temps = OrganismTemperatures(core_temperature, skin_temperature, insulation_temperature)
            radiant_temp_result = radiant_temperature(;
                body,
                insulation,
                insulation_pars,
                org_temps,
                conductivities,
                side,
                substrate_conductance,
                longwave_depth_fraction,
                conduction_fraction,
                evaporation_flux=skin_evaporation_flux,
                substrate_temperature,
            )
            calculated_radiant_temperature = radiant_temp_result.radiant_temperature
            compressed_insulation_temperature = radiant_temp_result.compressed_insulation_temperature
            conductances = radiant_temp_result.conductances
            divisors = radiant_temp_result.divisors
            # Radiative heat fluxes
            sky_radiation_coeff = area_convection * view_factors.sky * 4 * ϵ_body * σ * ((calculated_radiant_temperature + sky_temperature) / 2)^3
            bush_radiation_coeff =
                area_convection * view_factors.bush * 4 * ϵ_body * σ * ((calculated_radiant_temperature + bush_temperature) / 2)^3
            vegetation_radiation_coeff =
                area_convection *
                view_factors.vegetation *
                4 *
                ϵ_body *
                σ *
                ((calculated_radiant_temperature + vegetation_temperature) / 2)^3
            ground_radiation_coeff =
                area_convection * view_factors.ground * 4 * ϵ_body * σ * ((calculated_radiant_temperature + ground_temperature) / 2)^3
            radiation_coeffs = RadiationCoeffs(sky_radiation_coeff, bush_radiation_coeff, vegetation_radiation_coeff, ground_radiation_coeff)
            env_temps = EnvironmentTemperatures(
                air_temperature, sky_temperature, ground_temperature, vegetation_temperature, bush_temperature, substrate_temperature
            )
            if conduction_fraction < 1
                # These calculations are for when there is less than 100% conduction.
                # The term insulation_evaporation_flux is included for heat lost due to evaporation from
                # the insulation surface
                insulation_temp_result = insulation_radiant_temperature(;
                    body,
                    insulation,
                    insulation_pars,
                    env_temps,
                    coeffs=radiation_coeffs,
                    conductances,
                    divisors,
                    side,
                    area_convection,
                    heat_transfer_coefficient,
                    substrate_conductance,
                    insulation_conductivity,
                    longwave_depth_fraction,
                    conduction_fraction,
                    solar_flux,
                    insulation_evaporation_flux,
                    core_temperature,
                    compressed_insulation_temperature,
                )
                insulation_temperature_calc = insulation_temp_result.calculated_insulation_temperature
                radiant_temperature2 = insulation_temp_result.radiant_temperature2

                sky_radiation_flux = sky_radiation_coeff * (radiant_temperature2 - sky_temperature)
                bush_radiation_flux = bush_radiation_coeff * (radiant_temperature2 - bush_temperature)
                vegetation_radiation_flux = vegetation_radiation_coeff * (radiant_temperature2 - vegetation_temperature)
                ground_radiation_flux = ground_radiation_coeff * (radiant_temperature2 - ground_temperature)
                longwave_flux = sky_radiation_flux + bush_radiation_flux + vegetation_radiation_flux + ground_radiation_flux
                convection_flux = heat_transfer_coefficient * area_convection * (insulation_temperature_calc - air_temperature)
                conduction_flux = u"W"(substrate_conductance * (compressed_insulation_temperature - substrate_temperature))
            else
                (; compressed_insulation_temperature) = compressed_radiant_temperature(;
                    body,
                    insulation,
                    insulation_pars,
                    conductivities,
                    side,
                    substrate_conductance,
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
                conduction_flux = substrate_conductance * (compressed_insulation_temperature - substrate_temperature)
                insulation_temperature_calc = compressed_insulation_temperature
            end
            environment_flux = longwave_flux + convection_flux + conduction_flux + insulation_evaporation_flux - solar_flux
            skin_temperature_mean, skin_temperature_calc1 = mean_skin_temperature(;
                body,
                insulation,
                insulation_pars,
                conductivities,
                conductances,
                conduction_fraction,
                environment_flux,
                skin_evaporation_flux,
                core_temperature,
                calculated_insulation_temperature=insulation_temperature_calc,
                compressed_insulation_temperature,
            )

            # Build fluxes (net_generated updated on success)
            fluxes = HeatFluxes(
                convection_flux,
                conduction_flux,
                net_generated,  # placeholder, updated below
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
                    net_generated = net_metabolic_heat(; body, conductivities, core_temperature, skin_temperature)
                    fluxes = setproperties(fluxes; net_generated)
                    return (;
                        insulation_temperature,
                        skin_temperature=skin_temperature_mean,
                        fluxes,
                        insulation_conductivity,
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
                        net_generated = net_metabolic_heat(; body, conductivities, core_temperature, skin_temperature)
                        fluxes = setproperties(fluxes; net_generated)
                        return (;
                            insulation_temperature,
                            skin_temperature=skin_temperature_mean,
                            fluxes,
                            insulation_conductivity,
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
                            insulation_conductivity,
                            tolerance,
                            success=false,
                            ntry,
                        )
                    end
                end
            end
        end
    end
    return (; insulation_temperature, skin_temperature, fluxes, insulation_conductivity, tolerance, success, ntry)
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
