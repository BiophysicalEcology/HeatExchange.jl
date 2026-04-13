# ectotherm / leaf heat balance

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
    )
    convection_heat_flow = convection_out.heat_flow

    return (;
        solar_flow, longwave_flow_in, longwave_flow_out,
        convection_heat_flow, convection_out, convection_area, conduction_area,
        solar_out, longwave_gain_out, longwave_loss_out,
    )
end

"""
    heat_balance(core_temperature, organism::Organism, e)

Calculate heat balance for an organism (animal or leaf) at a given core temperature.

Dispatches on `evaporation_pars(organism)` type:
- `AnimalEvaporationParameters` → ectotherm/animal heat balance
- `LeafEvaporationParameters`   → leaf heat balance (defined in `src/leaf/leaf.jl`)

# Arguments
- `core_temperature`: Core temperature to evaluate
- `organism::Organism`: Organism with body geometry and traits
- `e`: Environment containing `environment_pars` and `environment_vars`
"""
heat_balance(T, organism::Organism, e) =
    heat_balance(T, evaporation_pars(organism), organism, e)

# Animal/ectotherm dispatch
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
    T_lung_resp = clamp(core_temperature, u"K"(1.0u"°C"), u"K"(50.0u"°C"))
    rates = MetabolicRates(; metabolic=metabolic_heat_flow)
    atmos = AtmosphericConditions(environment_vars)
    respiration_out = respiration(
        rates,
        resp_pars,
        atmos,
        o.body.shape.mass,
        T_lung_resp,
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
        convection_out.mass,
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
        respiration_mass=respiration_out.respiration_mass,
        cutaneous_mass=evaporation_out.m_cut,
        eye_mass=evaporation_out.m_eyes,
    )
    return (;
        heat_balance=heat_balance_val,
        core_temperature,
        surface_temperature,
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

function heat_balance(core_temperature, ::Fur, o::Organism, e)
    # stub for insulated animals — not yet implemented
end

"""
    get_Tb(mod::Model, environment_pars, vars)

Find the equilibrium core temperature for an ectotherm.

Uses root-finding (bisection) to find the core temperature where the heat
balance equals zero.

# Arguments
- `mod::Model`: Organism model
- `environment_pars`: Environmental parameters
- `vars`: Environmental variables

# Returns
Full heat_balance output at the equilibrium core temperature.
"""
function get_Tb(mod::Model, environment_pars, vars)
    air_temperature = vars.environment.air_temperature
    core_temperature = find_zero(
        t -> heat_balance(t, mod, environment_pars, vars), (air_temperature - 40u"K", air_temperature + 100u"K"), Bisection()
    )
    heat_balance(core_temperature, mod, environment_pars, vars)
end

"""
    solve_temperature(organism, environment; T_bracket=(270.0u"K", 370.0u"K")) → T

Find the steady-state temperature of an organism (animal or leaf) by root-finding on
`heat_balance`. Dispatch (animal vs leaf) is determined automatically by
`evaporation_pars(organism)` type. Returns air temperature as fallback if root-finding fails.

Eye/stomata state is a fixed parameter of the organism — no switching here.
Temperature-triggered behavioural changes (eye opening, stomatal closure) belong in
BiophysicalBehaviour.jl.
"""
function solve_temperature(organism, environment; T_bracket=(270.0u"K", 370.0u"K"))
    lo = ustrip(u"K", T_bracket[1])
    hi = ustrip(u"K", T_bracket[2])
    try
        T_sol = zbrent(
            T -> ustrip(u"W", heat_balance(T * u"K", organism, environment).heat_balance),
            lo, hi, 1e-3,
        )
        T_sol * u"K"
    catch
        environment.environment_vars.air_temperature
    end
end
