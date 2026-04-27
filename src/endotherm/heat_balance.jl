"""
    heat_balance(T_core, T_skin, T_ins, Q_gen; ...)

Compute endotherm heat budget residuals at explicitly given temperatures and metabolic rate.

Non-iterative: all state variables are explicit inputs. Returns three residuals that
equal zero at a valid steady-state heat balance, making this function suitable as an
IPOPT/NLP constraint function differentiable via ForwardDiff.

This extends the existing `heat_balance(T_body, organism, e)` ectotherm dispatch with
a multi-argument endotherm method. Decision variables are positional arguments; fixed
parameters come from keyword arguments.

Designed to work per body side (`:dorsal` or `:ventral`) using the same pre-packed
`geometry_vars` and `environment_vars` structures that `solve_temperatures` accepts.

# Arguments
- `T_core`: Core body temperature (K) — setpoint or decision variable
- `T_skin`: Skin surface temperature (K) — decision variable (replaces iterative solve)
- `T_ins`: Insulation surface (fur-air interface) temperature (K) — decision variable
- `Q_gen`: Metabolic heat generation rate (W) — decision variable (replaces zbrent)

# Keywords
- `body::AbstractBody`: Body geometry (shape + composite insulation)
- `insulation_pars::InsulationParameters`: Insulation parameters (fibre properties, depths)
- `insulation::InsulationProperties`: Precomputed insulation properties; temperature-sensitive
  conductivities are recomputed internally from `T_skin` and `T_ins`
- `geometry_vars::GeometryVariables`: Geometric variables (side, conductance coefficient, etc.)
- `environment_vars::NamedTuple`: Packed environment (temperatures, view factors, atmos, solar)
- `traits::NamedTuple`: Fixed organism traits (fat conductivity, emissivity, evap fractions)
- `resp_pars`: Respiration parameters
- `minimum_metabolic_heat`: Minimum metabolic rate floor (W); defaults to zero
- `k_flesh`: Flesh thermal conductivity (W/m/K); overrides `traits.flesh_conductivity`
- `pant`: Panting multiplier; overrides `resp_pars.pant`
- `skin_wetness`: Skin wetness fraction; overrides `traits.skin_wetness`

# Returns
NamedTuple with heat fluxes and three residuals:
- `solar_heat_flow`, `convection_heat_flow`, `radiation_heat_flow`, `conduction_heat_flow`
- `skin_evaporation_heat_flow`, `insulation_evaporation_heat_flow`, `respiration_heat_flow`
- `net_metabolic_heat_internal`: heat conducted through flesh/fat layer
- `sky_radiation_flow`, `bush_radiation_flow`, `vegetation_radiation_flow`, `ground_radiation_flow`
- `residual_energy_balance` (W): Q_gen + Q_solar − Q_resp − Q_evap − Q_conv − Q_rad − Q_cond = 0
- `residual_internal_conduction` (W): (Q_gen − Q_resp) − net_metabolic_heat_internal = 0
- `residual_skin_temperature` (K): T_skin − T_skin_from_mean_skin_temperature = 0
"""
function heat_balance(
    T_core, T_skin, T_ins, Q_gen;
    body::AbstractBody,
    insulation_pars::InsulationParameters,
    insulation::InsulationProperties,
    geometry_vars::GeometryVariables,
    environment_vars::NamedTuple,
    traits::NamedTuple,
    resp_pars,
    minimum_metabolic_heat = zero(Q_gen),
    k_flesh = traits.flesh_conductivity,
    pant = resp_pars.pant,
    skin_wetness = traits.skin_wetness,
)
    (; side, conductance_coefficient, conduction_fraction, longwave_depth_fraction) = geometry_vars
    (;
        temperature, view_factors, atmos, fluid, solar_flow, gas_fractions, convection_enhancement,
    ) = environment_vars
    air_temperature     = temperature.air
    sky_temperature     = temperature.sky
    ground_temperature  = temperature.ground
    vegetation_temperature = temperature.vegetation
    bush_temperature    = temperature.bush
    substrate_temperature  = temperature.substrate
    (; relative_humidity, wind_speed, atmospheric_pressure) = atmos
    (; fat_conductivity, ϵ_body, insulation_wetness, bare_skin_fraction, eye_fraction) = traits

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    # Body areas
    total_area      = BiophysicalGeometry.total_area(body)
    area_evaporation = evaporation_area(body)
    area_convection  = total_area * (1 - conduction_fraction)

    # -------------------------------------------------------------------------
    # Recompute temperature-dependent insulation conductivity at current T_skin, T_ins.
    # This mirrors the per-iteration update in solve_with_insulation! and ensures
    # the function is differentiable with respect to the temperature decision variables.
    # -------------------------------------------------------------------------
    insulation_temp_mean = T_ins * 0.7 + T_skin * 0.3
    air_k = dry_air_properties(insulation_temp_mean).thermal_conductivity
    side_depth = getproperty(insulation.fibres, side).depth
    side_fibres = setproperties(getproperty(insulation_pars, side); depth = side_depth)
    side_thermal = insulation_thermal_conductivity(side_fibres, air_k)
    effective_conductivity = side_thermal.effective_conductivity

    absorption_coefficient = getproperty(insulation.absorption_coefficients, side)
    approx_radiant_temperature = T_skin * (1 - longwave_depth_fraction) + T_ins * longwave_depth_fraction
    radiative_conductivity = (16 * σ * approx_radiant_temperature^3) / (3 * absorption_coefficient)
    insulation_conductivity = effective_conductivity + radiative_conductivity

    # ThermalConductivities with k_flesh override (vasodilation decision variable)
    conductivities = ThermalConductivities(k_flesh, fat_conductivity, insulation_conductivity)

    # Update insulation struct with recomputed compressed conductivity (used in radiant_temperature)
    insulation_updated = setproperties(insulation; conductivity_compressed = effective_conductivity)

    # -------------------------------------------------------------------------
    # Convection at insulation surface (pure function of T_ins)
    # -------------------------------------------------------------------------
    conv = convection(;
        body,
        area = area_convection,
        air_temperature,
        surface_temperature = T_ins,
        wind_speed,
        atmospheric_pressure,
        fluid,
        gas_fractions,
        convection_enhancement,
    )
    heat_transfer_coefficient = conv.heat.combined

    # -------------------------------------------------------------------------
    # Skin evaporation (uses mass transfer coefficients from convection call)
    # -------------------------------------------------------------------------
    evap_pars_skin = AnimalEvaporationParameters(; skin_wetness, eye_fraction, bare_skin_fraction)
    atmos_local = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
    skin_evaporation_heat_flow = evaporation(
        evap_pars_skin,
        conv.mass,
        atmos_local,
        area_evaporation,
        T_skin,
        air_temperature;
        gas_fractions,
    ).evaporation_heat_flow

    # -------------------------------------------------------------------------
    # Radiant temperature at insulation depth (shape-dispatched pure function).
    # Provides T_rad for radiation coefficients, compressed_insulation_temperature
    # for conduction, and conductances/divisors for mean_skin_temperature.
    # -------------------------------------------------------------------------
    org_temps = OrganismTemperatures(T_core, T_skin, T_ins)
    radiant_temp_result = radiant_temperature(;
        body,
        insulation = insulation_updated,
        insulation_pars,
        org_temps,
        conductivities,
        side,
        conductance_coefficient,
        longwave_depth_fraction,
        conduction_fraction,
        evaporation_flow = skin_evaporation_heat_flow,
        substrate_temperature,
    )
    T_rad = radiant_temp_result.radiant_temperature
    compressed_insulation_temperature = radiant_temp_result.compressed_insulation_temperature
    conductances = radiant_temp_result.conductances
    divisors     = radiant_temp_result.divisors

    # -------------------------------------------------------------------------
    # Insulation surface evaporation (if insulation is wet — conditional is on
    # fixed parameters, so not inside the differentiable path)
    # -------------------------------------------------------------------------
    insulation_evaporation_heat_flow = if insulation_wetness > 0 && insulation.insulation_test > 0.0u"m"
        evap_pars_ins = AnimalEvaporationParameters(;
            skin_wetness = insulation_wetness,
            eye_fraction = 0.0,
            bare_skin_fraction = 1.0,
        )
        mass_ins = TransferCoefficients(conv.mass.combined, conv.mass.combined, conv.mass.forced)
        evaporation(
            evap_pars_ins, mass_ins, atmos_local, area_convection, T_ins, air_temperature;
            gas_fractions,
        ).evaporation_heat_flow
    else
        0.0u"W"
    end

    # -------------------------------------------------------------------------
    # Radiation exchange (coefficients linearised at T_rad, flows use T_rad)
    # -------------------------------------------------------------------------
    sky_radiation_coeff =
        area_convection * view_factors.sky * 4 * ϵ_body * σ * ((T_rad + sky_temperature) / 2)^3
    bush_radiation_coeff =
        area_convection * view_factors.bush * 4 * ϵ_body * σ * ((T_rad + bush_temperature) / 2)^3
    vegetation_radiation_coeff =
        area_convection * view_factors.vegetation * 4 * ϵ_body * σ *
        ((T_rad + vegetation_temperature) / 2)^3
    ground_radiation_coeff =
        area_convection * view_factors.ground * 4 * ϵ_body * σ * ((T_rad + ground_temperature) / 2)^3

    sky_radiation_flow        = sky_radiation_coeff        * (T_rad - sky_temperature)
    bush_radiation_flow       = bush_radiation_coeff       * (T_rad - bush_temperature)
    vegetation_radiation_flow = vegetation_radiation_coeff * (T_rad - vegetation_temperature)
    ground_radiation_flow     = ground_radiation_coeff     * (T_rad - ground_temperature)
    radiation_heat_flow =
        sky_radiation_flow + bush_radiation_flow + vegetation_radiation_flow + ground_radiation_flow

    # -------------------------------------------------------------------------
    # Convection flow (at outer insulation surface T_ins)
    # -------------------------------------------------------------------------
    convection_heat_flow = heat_transfer_coefficient * area_convection * (T_ins - air_temperature)

    # -------------------------------------------------------------------------
    # Conduction flow to substrate (via compressed insulation temperature)
    # -------------------------------------------------------------------------
    conduction_heat_flow = u"W"(conductance_coefficient * (compressed_insulation_temperature - substrate_temperature))

    # -------------------------------------------------------------------------
    # Net metabolic heat conducted through flesh and fat (existing shape-dispatched function)
    # -------------------------------------------------------------------------
    net_metabolic_heat_internal = net_metabolic_heat(;
        body,
        conductivities,
        core_temperature = T_core,
        skin_temperature = T_skin,
    )

    # -------------------------------------------------------------------------
    # Skin temperature from the insulation-side calculation.
    # Combines core-side (skin_temperature_calc1) and insulation-side
    # (skin_temperature_calc2) into a mean; both equal T_skin at convergence.
    # -------------------------------------------------------------------------
    environment_flow = radiation_heat_flow + convection_heat_flow + conduction_heat_flow +
                       insulation_evaporation_heat_flow - solar_flow
    mean_skin_temp_result = mean_skin_temperature(;
        body,
        insulation = insulation_updated,
        insulation_pars,
        conductivities,
        conductances,
        conduction_fraction,
        environment_flow,
        skin_evaporation_flow = skin_evaporation_heat_flow,
        core_temperature = T_core,
        calculated_insulation_temperature = T_ins,
        compressed_insulation_temperature,
    )
    T_skin_from_heat_balance = mean_skin_temp_result.mean_skin_temperature

    # -------------------------------------------------------------------------
    # Respiration with pant override (pure function of Q_gen and temperatures)
    # -------------------------------------------------------------------------
    lung_temperature = (T_core + T_skin) / 2
    resp_pars_effective = setproperties(resp_pars; pant)
    resp_out = respiration(
        MetabolicRates(; metabolic = Q_gen, sum = Q_gen, minimum = minimum_metabolic_heat),
        resp_pars_effective,
        atmos_local,
        body.shape.mass,
        lung_temperature,
        air_temperature;
        gas_fractions,
        O2conversion = Kleiber1961(),
    )
    respiration_heat_flow = uconvert(u"W", resp_out.respiration_heat_flow)

    # =========================================================================
    # Residuals — each equals zero at a valid steady-state heat balance
    # =========================================================================

    # Global power balance: total generation = total losses
    residual_energy_balance = Q_gen + solar_flow -
        respiration_heat_flow -
        skin_evaporation_heat_flow - insulation_evaporation_heat_flow -
        convection_heat_flow - radiation_heat_flow - conduction_heat_flow

    # Internal conduction: net metabolic heat (after respiration) must match
    # what the flesh/fat layer conducts given T_core and T_skin
    residual_internal_conduction = (Q_gen - respiration_heat_flow) - net_metabolic_heat_internal

    # Skin temperature: T_skin must be consistent with both the core-side and
    # insulation-side estimates from mean_skin_temperature (units: K)
    residual_skin_temperature = T_skin - T_skin_from_heat_balance

    return (;
        solar_heat_flow = solar_flow,
        convection_heat_flow,
        radiation_heat_flow,
        conduction_heat_flow,
        skin_evaporation_heat_flow,
        insulation_evaporation_heat_flow,
        respiration_heat_flow,
        net_metabolic_heat_internal,
        sky_radiation_flow,
        bush_radiation_flow,
        vegetation_radiation_flow,
        ground_radiation_flow,
        residual_energy_balance,
        residual_internal_conduction,
        residual_skin_temperature,
    )
end
