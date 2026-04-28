function _insulation_evaporation(conv, atmos, area_convection, insulation_temperature, air_temperature,
                                  insulation_wetness, insulation_test; gas_fractions)
    insulation_wetness > 0 && insulation_test > 0.0u"m" || return 0.0u"W"
    evap_pars_ins = AnimalEvaporationParameters(;
        skin_wetness = insulation_wetness,
        eye_fraction = 0.0,
        bare_skin_fraction = 1.0,
    )
    mass_ins = TransferCoefficients(conv.mass_transfer_coefficient.combined, conv.mass_transfer_coefficient.combined, conv.mass_transfer_coefficient.forced)
    evaporation(evap_pars_ins, mass_ins, atmos, area_convection, insulation_temperature, air_temperature;
                gas_fractions).evaporation_heat_flow
end

function _insulation_conductivity(insulation, insulation_pars, side, insulation_temperature, skin_temperature, longwave_depth_fraction, σ)
    insulation_temp_mean = insulation_temperature * 0.7 + skin_temperature * 0.3
    air_k = dry_air_properties(insulation_temp_mean).thermal_conductivity
    side_depth = getproperty(insulation.fibres, side).depth
    side_fibres = setproperties(getproperty(insulation_pars, side); depth = side_depth)
    effective_conductivity = insulation_thermal_conductivity(side_fibres, air_k).effective_conductivity
    absorption_coefficient = getproperty(insulation.absorption_coefficients, side)
    approx_radiant_temperature = skin_temperature * (1 - longwave_depth_fraction) + insulation_temperature * longwave_depth_fraction
    radiative_conductivity = (16 * σ * approx_radiant_temperature^3) / (3 * absorption_coefficient)
    return (; insulation_conductivity = effective_conductivity + radiative_conductivity, effective_conductivity)
end

function _radiation_coefficients(area, view_factors, ϵ, σ, radiant_temperature, env_temps)
    T = env_temps
    F = view_factors
    coeff(vf, t) = area * vf * 4 * ϵ * σ * ((radiant_temperature + t) / 2)^3
    RadiationCoeffs(coeff(F.sky, T.sky), coeff(F.bush, T.bush),
                    coeff(F.vegetation, T.vegetation), coeff(F.ground, T.ground))
end

"""
    solve_temperatures(; body, insulation_pars, insulation, geometry_vars, environment_vars, traits, temperature_tolerance, skin_temperature, insulation_temperature)

Simultaneously solve for skin and insulation surface temperatures.

Iteratively finds temperatures that satisfy heat balance equations for an endotherm,
accounting for convection, radiation, evaporation, and conduction through insulation.

# Keywords
- `body::AbstractBody`: Body geometry
- `insulation_pars::InsulationParameters`: Insulation parameters
- `insulation::InsulationProperties`: Computed insulation properties
- `geometry_vars::GeometryVariables`: Geometric variables (side, conduction_fraction, etc.)
- `environment_vars::NamedTuple`: Environmental variables (temperatures, wind, humidity, etc.)
- `traits::NamedTuple`: Organism traits (core_temperature, conductivities, emissivity, etc.)
- `temperature_tolerance`: Convergence tolerance for temperature iteration
- `skin_temperature`: Initial guess for skin temperature
- `insulation_temperature`: Initial guess for insulation surface temperature

# Returns
NamedTuple with:
- `insulation_temperature`: Converged insulation surface temperature
- `skin_temperature`: Converged mean skin temperature
- `flows::HeatFlows`: Heat flow components
- `insulation_conductivity`: Effective insulation conductivity
- `tolerance`: Final tolerance used
- `success`: Whether convergence was achieved
- `ntry`: Number of iterations
"""
function solve_temperatures(;
    body::AbstractBody,
    insulation_pars::InsulationParameters,
    insulation::InsulationProperties,
    geometry_vars::GeometryVariables,
    environment_vars::NamedTuple,
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
            geometry_vars,
            environment_vars,
            traits,
            temperature_tolerance,
            skin_temperature,
            insulation_temperature,
        )
    else
        return solve_without_insulation!(
            body,
            geometry_vars,
            environment_vars,
            traits,
            temperature_tolerance,
            skin_temperature,
            insulation_temperature,
        )
    end
end

function solve_without_insulation!(
    body::AbstractBody, geometry_vars::GeometryVariables, environment_vars::NamedTuple, traits::NamedTuple, temperature_tolerance, skin_temperature, insulation_temperature
)
    (;
        temperature,
        view_factors,
        atmos,
        fluid,
        solar_flow,
        gas_fractions,
        convection_enhancement,
    ) = environment_vars
    T = temperature
    F = view_factors
    air_temperature = T.air
    (; relative_humidity, wind_speed, atmospheric_pressure) = atmos
    (; core_temperature, flesh_conductivity, ϵ_body, skin_wetness, bare_skin_fraction, eye_fraction) = traits
    tolerance = temperature_tolerance
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    ntry = 0
    volume = flesh_volume(body)
    total_area = BiophysicalGeometry.total_area(body)
    area_evaporation = total_area
    area_convection = total_area #* (1 - conduction_fraction)
    r_skin = skin_radius(body)
    insulation_evaporation_heat_flow = 0.0u"W"
    conduction_flow = 0.0u"W"

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
            heat_transfer_coefficient = conv.heat_transfer_coefficient.combined
            evap_pars_local = AnimalEvaporationParameters(;
                skin_wetness,
                eye_fraction,
                bare_skin_fraction,
            )
            atmos_local = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
            skin_evaporation_flow = evaporation(
                evap_pars_local,
                conv.mass_transfer_coefficient,
                atmos_local,
                area_evaporation,
                skin_temperature,
                air_temperature;
                gas_fractions,
            ).evaporation_heat_flow

            # Radiation coefficients for radiant exchange
            _rc = _radiation_coefficients(area_convection, F, ϵ_body, σ, skin_temperature, T)
            sky_radiation_coeff        = _rc.sky
            bush_radiation_coeff       = _rc.bush
            vegetation_radiation_coeff = _rc.vegetation
            ground_radiation_coeff     = _rc.ground
            skin_temperature1 =
                ((4.0 * flesh_conductivity * volume) / (r_skin^2) * core_temperature) - skin_evaporation_flow +
                heat_transfer_coefficient * area_convection * T.air +
                solar_flow
            skin_temperature2 =
                sky_radiation_coeff * T.sky + bush_radiation_coeff * T.bush + vegetation_radiation_coeff * T.vegetation + ground_radiation_coeff * T.ground
            skin_temperature3 =
                ((4.0 * flesh_conductivity * volume) / (r_skin^2)) +
                heat_transfer_coefficient * area_convection +
                sky_radiation_coeff +
                bush_radiation_coeff +
                vegetation_radiation_coeff +
                ground_radiation_coeff

            skin_temperature_calc = (skin_temperature1 + skin_temperature2) / skin_temperature3

            sky_radiation_flow = sky_radiation_coeff * (skin_temperature_calc - T.sky)
            bush_radiation_flow = bush_radiation_coeff * (skin_temperature_calc - T.bush)
            vegetation_radiation_flow = vegetation_radiation_coeff * (skin_temperature_calc - T.vegetation)
            ground_radiation_flow = ground_radiation_coeff * (skin_temperature_calc - T.ground)

            longwave_flow = sky_radiation_flow + bush_radiation_flow + vegetation_radiation_flow + ground_radiation_flow
            convection_flow = heat_transfer_coefficient * area_convection * (skin_temperature_calc - T.air)

            # Build flows (net_metabolic updated on success)
            flows = HeatFlows(
                convection_flow,
                conduction_flow,
                0.0u"W",  # net_metabolic placeholder
                skin_evaporation_flow,
                insulation_evaporation_heat_flow,
                longwave_flow,
                solar_flow,
                sky_radiation_flow,
                bush_radiation_flow,
                vegetation_radiation_flow,
                ground_radiation_flow,
            )

            Δskin_temperature = abs(skin_temperature - skin_temperature_calc)

            if Δskin_temperature < tolerance
                #TODO why is this not shape-specific?
                net_metabolic = (4 * flesh_conductivity * volume / r_skin^2) * (core_temperature - skin_temperature_calc)
                flows = setproperties(flows; net_metabolic)
                return (;
                    insulation_temperature,
                    skin_temperature=skin_temperature_calc,
                    flows,
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
                            flows,
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
    geometry_vars::GeometryVariables,
    environment_vars::NamedTuple,
    traits::NamedTuple,
    temperature_tolerance,
    skin_temperature,
    insulation_temperature,
)
    (; side, conductance_coefficient, ventral_fraction, conduction_fraction, longwave_depth_fraction) = geometry_vars
    (;
        temperature,
        view_factors,
        atmos,
        fluid,
        solar_flow,
        gas_fractions,
        convection_enhancement,
    ) = environment_vars
    env_temps = temperature
    T = env_temps
    F = view_factors
    air_temperature = T.air
    substrate_temperature = T.substrate
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
    total_area = BiophysicalGeometry.total_area(body)
    area_convection = total_area * (1 - conduction_fraction)
    insulation_test = insulation.insulation_test

    ntry = 0
    solct = 0
    solution_procedure = 1
    success = true
    net_metabolic = 0.0u"W"

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
            heat_transfer_coefficient = conv.heat_transfer_coefficient.combined
            evap_pars_skin = AnimalEvaporationParameters(;
                skin_wetness,
                eye_fraction,
                bare_skin_fraction,
            )
            atmos_skin = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
            skin_evaporation_flow = evaporation(
                evap_pars_skin,
                conv.mass_transfer_coefficient,
                atmos_skin,
                area_evaporation,
                skin_temperature,
                air_temperature;
                gas_fractions,
            ).evaporation_heat_flow
            # second from insulation
            insulation_evaporation_heat_flow = _insulation_evaporation(
                conv, atmos_skin, area_convection, insulation_temperature, air_temperature,
                insulation_wetness, insulation_test; gas_fractions)
            # Recompute insulation thermal properties for current temperatures
            (; insulation_conductivity, effective_conductivity) = _insulation_conductivity(
                insulation, insulation_pars, side, insulation_temperature, skin_temperature,
                longwave_depth_fraction, σ)
            side_conductivities = if side == :dorsal
                setproperties(insulation.conductivities; dorsal=effective_conductivity)
            else
                setproperties(insulation.conductivities; ventral=effective_conductivity)
            end
            insulation = setproperties(insulation;
                conductivity_compressed=effective_conductivity,
                conductivities=side_conductivities,
            )
            conductivities = ThermalConductivities(flesh_conductivity, fat_conductivity, insulation_conductivity)
            org_temps = OrganismTemperatures(core_temperature, skin_temperature, insulation_temperature)
            radiant_temp_result = radiant_temperature(;
                body,
                insulation,
                insulation_pars,
                org_temps,
                conductivities,
                side,
                conductance_coefficient,
                longwave_depth_fraction,
                conduction_fraction,
                evaporation_flow=skin_evaporation_flow,
                substrate_temperature,
            )
            calculated_radiant_temperature = radiant_temp_result.radiant_temperature
            compressed_insulation_temperature = radiant_temp_result.compressed_insulation_temperature
            conductances = radiant_temp_result.conductances
            divisors = radiant_temp_result.divisors
            # Radiative heat flows
            radiation_coeffs = _radiation_coefficients(area_convection, F, ϵ_body, σ,
                calculated_radiant_temperature, T)
            sky_radiation_coeff        = radiation_coeffs.sky
            bush_radiation_coeff       = radiation_coeffs.bush
            vegetation_radiation_coeff = radiation_coeffs.vegetation
            ground_radiation_coeff     = radiation_coeffs.ground
            if conduction_fraction < 1
                # These calculations are for when there is less than 100% conduction.
                # The term insulation_evaporation_heat_flow is included for heat lost due to evaporation from
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
                    conductance_coefficient,
                    insulation_conductivity,
                    longwave_depth_fraction,
                    conduction_fraction,
                    solar_flow,
                    insulation_evaporation_heat_flow,
                    core_temperature,
                    compressed_insulation_temperature,
                )
                insulation_temperature_calc = insulation_temp_result.calculated_insulation_temperature
                radiant_temperature2 = insulation_temp_result.radiant_temperature2

                sky_radiation_flow = sky_radiation_coeff * (radiant_temperature2 - T.sky)
                bush_radiation_flow = bush_radiation_coeff * (radiant_temperature2 - T.bush)
                vegetation_radiation_flow = vegetation_radiation_coeff * (radiant_temperature2 - T.vegetation)
                ground_radiation_flow = ground_radiation_coeff * (radiant_temperature2 - T.ground)
                longwave_flow = sky_radiation_flow + bush_radiation_flow + vegetation_radiation_flow + ground_radiation_flow
                convection_flow = heat_transfer_coefficient * area_convection * (insulation_temperature_calc - T.air)
                conduction_flow = u"W"(conductance_coefficient * (compressed_insulation_temperature - substrate_temperature))
            else
                (; compressed_insulation_temperature) = compressed_radiant_temperature(;
                    body,
                    insulation,
                    insulation_pars,
                    conductivities,
                    side,
                    conductance_coefficient,
                    core_temperature,
                    substrate_temperature,
                )
                sky_radiation_flow = 0.0u"W"
                bush_radiation_flow = 0.0u"W"
                vegetation_radiation_flow = 0.0u"W"
                ground_radiation_flow = 0.0u"W"
                longwave_flow = 0.0u"W"
                convection_flow = 0.0u"W"
                insulation_evaporation_heat_flow = 0.0u"W"
                solar_flow = 0.0u"W"
                conduction_flow = conductance_coefficient * (compressed_insulation_temperature - substrate_temperature)
                insulation_temperature_calc = compressed_insulation_temperature
            end
            environment_flow = longwave_flow + convection_flow + conduction_flow + insulation_evaporation_heat_flow - solar_flow
            skin_temperature_mean, skin_temperature_calc1 = mean_skin_temperature(;
                body,
                insulation,
                insulation_pars,
                conductivities,
                conductances,
                conduction_fraction,
                environment_flow,
                skin_evaporation_flow,
                core_temperature,
                calculated_insulation_temperature=insulation_temperature_calc,
                compressed_insulation_temperature,
            )

            # Build flows (net_metabolic updated on success)
            flows = HeatFlows(
                convection_flow,
                conduction_flow,
                net_metabolic,  # placeholder, updated below
                skin_evaporation_flow,
                insulation_evaporation_heat_flow,
                longwave_flow,
                solar_flow,
                sky_radiation_flow,
                bush_radiation_flow,
                vegetation_radiation_flow,
                ground_radiation_flow,
            )

            Δinsulation_temperature = abs(insulation_temperature - insulation_temperature_calc)
            Δskin_temperature = abs(skin_temperature - skin_temperature_mean)

            # first convergence test (Δinsulation_temperature)
            if Δinsulation_temperature < tolerance
                # Next check skin_temperature convergence
                if Δskin_temperature < tolerance
                    net_metabolic = net_metabolic_heat(; body, conductivities, core_temperature, skin_temperature)
                    flows = setproperties(flows; net_metabolic)
                    return (;
                        insulation_temperature,
                        skin_temperature=skin_temperature_mean,
                        flows,
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
                        net_metabolic = net_metabolic_heat(; body, conductivities, core_temperature, skin_temperature)
                        flows = setproperties(flows; net_metabolic)
                        return (;
                            insulation_temperature,
                            skin_temperature=skin_temperature_mean,
                            flows,
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
                            flows,
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
    return (; insulation_temperature, skin_temperature, flows, insulation_conductivity, tolerance, success, ntry)
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
