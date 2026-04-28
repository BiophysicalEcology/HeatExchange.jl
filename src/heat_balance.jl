# Heat balance for all organism types.
# Outer dispatch routes through EvaluationStrategy and evaporation type.
# _radiative_convective_flows handles the shared single-body physics.

"""
    _radiative_convective_flows(surface_temperature, organism, environment_pars, environment_vars)

Private helper: compute solar, longwave, and convective heat flows for a body surface.

Returns a NamedTuple:
- `solar_flow` — absorbed solar radiation (W)
- `longwave_flow_in` — incoming longwave radiation (W)
- `longwave_flow_out` — outgoing longwave radiation (W)
- `convection_heat_flow` — convective heat loss (W)
- `convection_out` — full convection output (for mass transfer coefficients etc.)
- `convection_area` — area used for convection (m²)
- `conduction_area` — area in contact with substrate (m²)
- `solar_out` — full solar output
- `longwave_gain_out` — full longwave-in output
- `longwave_loss_out` — full longwave-out output
"""
function _radiative_convective_flows(surface_temperature, o::Organism, environment_pars, environment_vars)
    external_conduction = conduction_pars_external(o)
    rad_pars = radiation_pars(o)
    conv_pars = convection_pars(o)

    total_area = BiophysicalGeometry.total_area(o.body)
    convection_area = total_area * (1 - external_conduction.conduction_fraction)
    conduction_area = total_area * external_conduction.conduction_fraction
    silhouette_area = BiophysicalGeometry.silhouette_area(o.body, rad_pars.solar_orientation, environment_vars.zenith_angle)

    absorptivities = Absorptivities(rad_pars, environment_pars)
    emissivities = Emissivities(rad_pars, environment_pars)
    view_factors = ViewFactors(rad_pars.sky_view_factor, rad_pars.ground_view_factor, 0.0, 0.0)
    solar_conditions = SolarConditions(environment_vars)
    environmental_temps = EnvironmentTemperatures(environment_vars)

    solar_out = solar(o.body, absorptivities, view_factors, solar_conditions, silhouette_area, conduction_area)
    solar_flow = solar_out.solar_flow

    longwave_gain_out = radiation_in(
        o.body, view_factors, emissivities, environmental_temps;
        conduction_fraction=external_conduction.conduction_fraction,
    )
    longwave_flow_in = longwave_gain_out.longwave_flow_in

    longwave_loss_out = radiation_out(
        o.body, view_factors, emissivities,
        external_conduction.conduction_fraction,
        surface_temperature, surface_temperature,
    )
    longwave_flow_out = longwave_loss_out.longwave_flow_out

    convection_out = convection(;
        body=o.body,
        area=convection_area,
        air_temperature=environment_vars.air_temperature,
        surface_temperature,
        wind_speed=environment_vars.wind_speed,
        atmospheric_pressure=environment_vars.atmospheric_pressure,
        fluid=environment_pars.fluid,
        gas_fractions=environment_pars.gas_fractions,
        convection_enhancement=environment_pars.convection_enhancement,
        characteristic_dimension_formula=conv_pars.characteristic_dimension_formula,
    )
    convection_heat_flow = convection_out.convection_flow

    return (;
        solar_flow, longwave_flow_in, longwave_flow_out,
        convection_heat_flow, convection_out, convection_area, conduction_area,
        solar_out, longwave_gain_out, longwave_loss_out,
    )
end

"""
    heat_balance(core_temperature, organism::Organism, e)

Calculate heat balance for an organism (animal or leaf) at a given core temperature.

Dispatches first on `evaluation_strategy(organism)`:
- `SingleBody` — single-body physics (Naked animals and leaves); further dispatches on
  `evaporation_pars` type and insulation type.
- `MultiSided` — dorsal/ventral split for insulated animals; uses `_pack_sides` +
  `solve_temperatures` to converge skin/insulation temperatures, then returns the combined
  residual `metabolic_heat_flow − respiration_heat_flow − net_metabolic`.

# Arguments
- `core_temperature`: Core temperature to evaluate
- `organism::Organism`: Organism with body geometry and traits
- `e`: Environment containing `environment_pars` and `environment_vars`
"""
heat_balance(T, organism::Organism, e) =
    heat_balance(T, evaluation_strategy(organism), evaporation_pars(organism), organism, e)

# SingleBody: drop the strategy, dispatch on evaporation type
heat_balance(T, ::SingleBody, ep, o::Organism, e) = heat_balance(T, ep, o, e)

# Animal/ectotherm SingleBody: dispatch on insulation type
heat_balance(T, ::AnimalEvaporationParameters, o::Organism, e) =
    heat_balance(T, insulation(body(o)), o, e)


function heat_balance(core_temperature, ::Naked, o::Organism, e)
    environment_pars = stripparams(e.environment_pars)
    environment_vars = e.environment_vars
    internal_conduction = conduction_pars_internal(o)
    evap_pars = evaporation_pars(o)
    hyd_pars = hydraulic_pars(o)
    resp_pars = respiration_pars(o)
    metab_pars = metabolism_pars(o)

    # metabolism
    metabolic_heat_flow = metabolic_rate(metab_pars.model, o.body.shape.mass, core_temperature)

    # respiration
    rates = MetabolicRates(; metabolic=metabolic_heat_flow)
    atmos = AtmosphericConditions(environment_vars)
    respiration_out = respiration(
        rates,
        resp_pars,
        atmos,
        o.body.shape.mass,
        clamp(core_temperature, u"K"(1.0u"°C"), u"K"(50.0u"°C")),
        environment_vars.air_temperature;
        gas_fractions=environment_pars.gas_fractions,
    )
    respiration_heat_flow = respiration_out.respiration_heat_flow
    oxygen_flow = u"ml/hr"(Joules_to_O2(metabolic_heat_flow))

    # net metabolic heat generation → surface temperature
    net_metabolic_heat_production = metabolic_heat_flow - respiration_heat_flow
    specific_metabolic_heat_production = net_metabolic_heat_production / o.body.geometry.volume
    (; surface_temperature, lung_temperature) = surface_and_lung_temperature(;
        body=o.body,
        flesh_conductivity=internal_conduction.flesh_conductivity,
        specific_metabolic_heat_production,
        core_temperature,
    )

    # radiative + convective flows
    flows = _radiative_convective_flows(surface_temperature, o, environment_pars, environment_vars)
    (; solar_flow, longwave_flow_in, longwave_flow_out, convection_heat_flow,
       convection_out, convection_area, conduction_area,
       solar_out, longwave_gain_out, longwave_loss_out) = flows

    # conduction
    conduction_flow = conduction(;
        conduction_area,
        L=environment_pars.conduction_depth,
        surface_temperature,
        substrate_temperature=environment_vars.substrate_temperature,
        substrate_conductivity=environment_vars.substrate_conductivity,
    )

    # evaporation — mouth opens when panting
    evap_eff = if resp_pars.pant > 1
        AnimalEvaporationParameters(;
            skin_wetness        = min(1.0, evap_pars.skin_wetness + resp_pars.mouth_fraction),
            insulation_wetness  = evap_pars.insulation_wetness,
            eye_fraction        = evap_pars.eye_fraction,
            bare_skin_fraction  = evap_pars.bare_skin_fraction,
            insulation_fraction = evap_pars.insulation_fraction,
        )
    else
        evap_pars
    end
    evaporation_out = evaporation(
        evap_eff,
        convection_out.mass_transfer_coefficient,
        atmos,
        convection_area,
        surface_temperature,
        environment_vars.air_temperature;
        water_potential=hyd_pars.water_potential,
        gas_fractions=environment_pars.gas_fractions,
    )
    evaporation_heat_flow = evaporation_out.evaporation_heat_flow

    # heat balance
    heat_flow_in  = solar_flow + longwave_flow_in + metabolic_heat_flow
    heat_flow_out = longwave_flow_out + convection_heat_flow + evaporation_heat_flow +
                    respiration_heat_flow + conduction_flow
    heat_balance_val = heat_flow_in - heat_flow_out

    energy_balance = (;
        solar_flow, longwave_flow_in, metabolic_heat_flow, respiration_heat_flow,
        evaporation_heat_flow, longwave_flow_out, convection_heat_flow, conduction_flow,
        heat_balance=heat_balance_val,
    )
    mass_balance = (;
        oxygen_flow,
        respiration_mass=respiration_out.respiration_mass_flow,
        cutaneous_mass=evaporation_out.cutaneous_mass_flow,
        eye_mass=evaporation_out.eyes_mass_flow,
    )
    return (;
        heat_balance=heat_balance_val,
        core_temperature,
        surface_temperature,
        skin_temperature=surface_temperature,
        insulation_temperature=surface_temperature,
        lung_temperature,
        energy_balance,
        mass_balance,
        respiration_out,
        solar_out,
        longwave_gain_out,
        longwave_loss_out,
        convection_out,
        evaporation_out,
    )
end

# ---------------------------------------------------------------------------
# Leaf heat balance
# ---------------------------------------------------------------------------

"""
    heat_balance(core_temperature, evap_pars::LeafEvaporationParameters, organism::Organism, e)

Calculate heat balance for a leaf at a given core temperature.

The leaf is treated as a thin, isothermal surface: a high `flesh_conductivity` in the
`InternalConductionParameters` will make surface ≈ internal temperature. For thick succulent
stems (cactus), set flesh_conductivity to the appropriate tissue value.

Transpiration is computed from stomatal vapour conductances via
`evaporation(::LeafEvaporationParameters, ...)`.

Dark respiration heat is included in `metabolic_heat_flow` via the metabolic model
(e.g. `PlantDarkRespiration`); `respiration_heat_flow` in `energy_balance` is `0.0 W`
(leaves have no pulmonary respiratory heat loss).
"""
function heat_balance(core_temperature, evap_pars::LeafEvaporationParameters, o::Organism, e)
    environment_pars = stripparams(e.environment_pars)
    environment_vars = e.environment_vars
    internal_conduction = conduction_pars_internal(o)
    hyd_pars = hydraulic_pars(o)
    metab_pars = metabolism_pars(o)

    # metabolism (dark respiration or nothing)
    metabolic_heat_flow = metabolic_rate(metab_pars.model, o.body.shape.mass, core_temperature)

    # net specific metabolic heat → surface temperature
    specific_metabolic_heat_production = metabolic_heat_flow / o.body.geometry.volume
    (; surface_temperature) = surface_and_lung_temperature(;
        body=o.body,
        flesh_conductivity=internal_conduction.flesh_conductivity,
        specific_metabolic_heat_production,
        core_temperature,
    )

    # radiative + convective flows
    flows = _radiative_convective_flows(surface_temperature, o, environment_pars, environment_vars)
    (; solar_flow, longwave_flow_in, longwave_flow_out,
       convection_heat_flow, convection_out, convection_area,
       solar_out, longwave_gain_out, longwave_loss_out) = flows

    # transpiration
    atmos = AtmosphericConditions(environment_vars)
    evaporation_out = evaporation(
        evap_pars,
        convection_out.mass_transfer_coefficient,
        atmos,
        convection_area,
        surface_temperature,
        environment_vars.air_temperature;
        water_potential=hyd_pars.water_potential,
        gas_fractions=environment_pars.gas_fractions,
    )
    evaporation_heat_flow = evaporation_out.evaporation_heat_flow

    heat_balance_val = solar_flow + longwave_flow_in + metabolic_heat_flow -
                       longwave_flow_out - convection_heat_flow - evaporation_heat_flow

    energy_balance = (;
        solar_flow, longwave_flow_in, longwave_flow_out,
        convection_heat_flow, evaporation_heat_flow,
        metabolic_heat_flow,
        respiration_heat_flow = 0.0u"W",
        conduction_flow = 0.0u"W",
        heat_balance=heat_balance_val,
    )
    oxygen_consumption_rate = u"ml/hr"(Joules_to_O2(Kleiber1961(), metabolic_heat_flow, 1.0))
    mass_balance = (; transpiration_mass=evaporation_out.transpiration_mass_flow, oxygen_consumption_rate)

    return (;
        heat_balance=heat_balance_val,
        core_temperature,
        surface_temperature,
        skin_temperature=surface_temperature,
        insulation_temperature=surface_temperature,
        energy_balance,
        mass_balance,
        evaporation_out,
        solar_out,
        longwave_gain_out,
        longwave_loss_out,
        convection_out,
    )
end

# ---------------------------------------------------------------------------
# Per-side insulated heat balance (4-arg, non-iterative)
# Takes explicit core/skin/insulation temperatures and metabolic_heat_flow; pre-packed per-side structures.
# Suitable as an NLP constraint function (differentiable via ForwardDiff).
# ---------------------------------------------------------------------------

"""
    heat_balance(core_temperature, skin_temperature, insulation_temperature, metabolic_heat_flow; ...)

Compute endotherm heat budget residuals at explicitly given temperatures and metabolic rate.

Non-iterative: all state variables are explicit inputs. Returns three residuals that
equal zero at a valid steady-state heat balance, making this function suitable as an
IPOPT/NLP constraint function differentiable via ForwardDiff.

This extends the existing `heat_balance(body_temperature, organism, e)` ectotherm dispatch with
a multi-argument endotherm method. Decision variables are positional arguments; fixed
parameters come from keyword arguments.

Designed to work per body side (`:dorsal` or `:ventral`) using the same pre-packed
`geometry_vars` and `environment_vars` structures that `solve_temperatures` accepts.

# Arguments
- `core_temperature`: Core body temperature (K) — setpoint or decision variable
- `skin_temperature`: Skin surface temperature (K) — decision variable (replaces iterative solve)
- `insulation_temperature`: Insulation surface (fur-air interface) temperature (K) — decision variable
- `metabolic_heat_flow`: Metabolic heat generation rate (W) — decision variable (replaces zbrent)

# Keywords
- `body::AbstractBody`: Body geometry (shape + composite insulation)
- `insulation_pars::InsulationParameters`: Insulation parameters (fibre properties, depths)
- `insulation::InsulationProperties`: Precomputed insulation properties; temperature-sensitive
  conductivities are recomputed internally from `skin_temperature` and `insulation_temperature`
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
- `residual_energy_balance` (W): metabolic + solar − resp − evap − conv − rad − cond = 0
- `residual_internal_conduction` (W): (metabolic − resp) − net_metabolic_heat_internal = 0
- `residual_skin_temperature` (K): skin_temperature − skin_temperature_from_heat_balance = 0
"""
function heat_balance(
    core_temperature, skin_temperature, insulation_temperature, metabolic_heat_flow;
    body::AbstractBody,
    insulation_pars::InsulationParameters,
    insulation::InsulationProperties,
    geometry_vars::GeometryVariables,
    environment_vars::NamedTuple,
    traits::NamedTuple,
    resp_pars,
    minimum_metabolic_heat = zero(metabolic_heat_flow),
    k_flesh = traits.flesh_conductivity,
    pant = resp_pars.pant,
    skin_wetness = traits.skin_wetness,
)
    (; side, conductance_coefficient, conduction_fraction, longwave_depth_fraction) = geometry_vars
    (;
        temperature, view_factors, atmos, fluid, solar_flow, gas_fractions, convection_enhancement,
    ) = environment_vars
    env_temps = temperature
    T = env_temps
    F = view_factors
    air_temperature   = T.air
    substrate_temperature = T.substrate
    (; relative_humidity, wind_speed, atmospheric_pressure) = atmos
    (; fat_conductivity, ϵ_body, insulation_wetness, bare_skin_fraction, eye_fraction) = traits

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    # Body areas
    total_area      = BiophysicalGeometry.total_area(body)
    area_evaporation = evaporation_area(body)
    area_convection  = total_area * (1 - conduction_fraction)

    # Recompute temperature-dependent insulation conductivity at current temperatures.
    (; insulation_conductivity, effective_conductivity) = _insulation_conductivity(
        insulation, insulation_pars, side, insulation_temperature, skin_temperature, longwave_depth_fraction, σ)
    conductivities = ThermalConductivities(k_flesh, fat_conductivity, insulation_conductivity)
    insulation_updated = setproperties(insulation; conductivity_compressed = effective_conductivity)

    # -------------------------------------------------------------------------
    # Convection at insulation surface
    # -------------------------------------------------------------------------
    conv = convection(;
        body,
        area = area_convection,
        air_temperature,
        surface_temperature = insulation_temperature,
        wind_speed,
        atmospheric_pressure,
        fluid,
        gas_fractions,
        convection_enhancement,
    )
    heat_transfer_coefficient = conv.heat_transfer_coefficient.combined

    # -------------------------------------------------------------------------
    # Skin evaporation (uses mass transfer coefficients from convection call)
    # -------------------------------------------------------------------------
    evap_pars_skin = AnimalEvaporationParameters(; skin_wetness, eye_fraction, bare_skin_fraction)
    atmos_local = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
    skin_evaporation_heat_flow = evaporation(
        evap_pars_skin,
        conv.mass_transfer_coefficient,
        atmos_local,
        area_evaporation,
        skin_temperature,
        air_temperature;
        gas_fractions,
    ).evaporation_heat_flow

    # -------------------------------------------------------------------------
    # Radiant temperature at insulation depth (shape-dispatched pure function).
    # -------------------------------------------------------------------------
    org_temps = OrganismTemperatures(core_temperature, skin_temperature, insulation_temperature)
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
    radiant_temp = radiant_temp_result.radiant_temperature
    compressed_insulation_temperature = radiant_temp_result.compressed_insulation_temperature
    conductances = radiant_temp_result.conductances

    # Insulation surface evaporation (conditional on fixed params, outside differentiable path)
    insulation_evaporation_heat_flow = _insulation_evaporation(
        conv, atmos_local, area_convection, insulation_temperature, air_temperature,
        insulation_wetness, insulation.insulation_test; gas_fractions)

    # -------------------------------------------------------------------------
    # Radiation exchange (coefficients linearised at radiant_temp, flows use radiant_temp)
    # -------------------------------------------------------------------------
    _rc = _radiation_coefficients(area_convection, F, ϵ_body, σ, radiant_temp, T)
    sky_radiation_coeff        = _rc.sky
    bush_radiation_coeff       = _rc.bush
    vegetation_radiation_coeff = _rc.vegetation
    ground_radiation_coeff     = _rc.ground

    sky_radiation_flow        = sky_radiation_coeff        * (radiant_temp - T.sky)
    bush_radiation_flow       = bush_radiation_coeff       * (radiant_temp - T.bush)
    vegetation_radiation_flow = vegetation_radiation_coeff * (radiant_temp - T.vegetation)
    ground_radiation_flow     = ground_radiation_coeff     * (radiant_temp - T.ground)
    radiation_heat_flow =
        sky_radiation_flow + bush_radiation_flow + vegetation_radiation_flow + ground_radiation_flow

    # -------------------------------------------------------------------------
    # Convection flow (at outer insulation surface)
    # -------------------------------------------------------------------------
    convection_heat_flow = heat_transfer_coefficient * area_convection * (insulation_temperature - air_temperature)

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
        core_temperature,
        skin_temperature,
    )

    # -------------------------------------------------------------------------
    # Skin temperature from the insulation-side calculation.
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
        core_temperature,
        calculated_insulation_temperature = insulation_temperature,
        compressed_insulation_temperature,
    )
    skin_temperature_from_heat_balance = mean_skin_temp_result.mean_skin_temperature

    # -------------------------------------------------------------------------
    # Respiration with pant override
    # -------------------------------------------------------------------------
    lung_temperature = (core_temperature + skin_temperature) / 2
    resp_pars_effective = setproperties(resp_pars; pant)
    resp_out = respiration(
        MetabolicRates(; metabolic = metabolic_heat_flow, sum = metabolic_heat_flow, minimum = minimum_metabolic_heat),
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

    residual_energy_balance = metabolic_heat_flow + solar_flow -
        respiration_heat_flow -
        skin_evaporation_heat_flow - insulation_evaporation_heat_flow -
        convection_heat_flow - radiation_heat_flow - conduction_heat_flow

    residual_internal_conduction = (metabolic_heat_flow - respiration_heat_flow) - net_metabolic_heat_internal

    residual_skin_temperature = skin_temperature - skin_temperature_from_heat_balance

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
