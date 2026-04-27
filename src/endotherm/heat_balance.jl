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

    # Recompute temperature-dependent insulation conductivity at current T_skin, T_ins.
    # Mirrors the per-iteration update in solve_with_insulation! for differentiability.
    (; insulation_conductivity, effective_conductivity) = _insulation_conductivity(
        insulation, insulation_pars, side, T_ins, T_skin, longwave_depth_fraction, σ)
    conductivities = ThermalConductivities(k_flesh, fat_conductivity, insulation_conductivity)
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

    # Insulation surface evaporation (conditional on fixed params, outside differentiable path)
    insulation_evaporation_heat_flow = _insulation_evaporation(
        conv, atmos_local, area_convection, T_ins, air_temperature,
        insulation_wetness, insulation.insulation_test; gas_fractions)

    # -------------------------------------------------------------------------
    # Radiation exchange (coefficients linearised at T_rad, flows use T_rad)
    # -------------------------------------------------------------------------
    _rc = _radiation_coefficients(area_convection, view_factors, ϵ_body, σ, T_rad,
        sky_temperature, bush_temperature, vegetation_temperature, ground_temperature)
    sky_radiation_coeff        = _rc.sky
    bush_radiation_coeff       = _rc.bush
    vegetation_radiation_coeff = _rc.vegetation
    ground_radiation_coeff     = _rc.ground

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
