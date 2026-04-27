# Leaf heat balance dispatch

"""
    heat_balance(leaf_temperature, evap_pars::LeafEvaporationParameters, organism::Organism, e)

Calculate heat balance for a leaf at a given leaf temperature.

The leaf is treated as a thin, isothermal surface: a high `flesh_conductivity` in the
`InternalConductionParameters` will make surface ≈ internal temperature. For thick succulent
stems (cactus), set flesh_conductivity to the appropriate tissue value.

Transpiration is computed from stomatal vapour conductances via
`evaporation(::LeafEvaporationParameters, ...)`.

# Arguments
- `leaf_temperature`: Leaf temperature to evaluate (K)
- `evap_pars::LeafEvaporationParameters`: Stomatal and cuticular conductances
- `organism::Organism`: Organism with body geometry and traits
- `e`: Environment containing `environment_pars` and `environment_vars`

# Returns
NamedTuple with:
- `heat_balance`: Net heat balance (W; zero at steady-state)
- `leaf_temperature`: Input leaf temperature (K)
- `surface_temperature`: Leaf surface temperature (K)
- `energy_balance`: NamedTuple of heat flow components
- `mass_balance`: NamedTuple with `transpiration_mass` (g/s) and `oxygen_consumption_rate` (ml/hr)
- `evaporation_out`: Full evaporation output
- `solar_out`, `longwave_gain_out`, `longwave_loss_out`, `convection_out`: Detailed outputs
"""
function heat_balance(leaf_temperature, evap_pars::LeafEvaporationParameters, o::Organism, e)
    environment_pars = stripparams(e.environment_pars)
    environment_vars = e.environment_vars
    internal_conduction = conduction_pars_internal(o)
    hyd_pars = hydraulic_pars(o)
    metab_pars = metabolism_pars(o)

    # metabolism (dark respiration or nothing)
    metabolic_heat_flow = metabolic_rate(metab_pars.model, o.body.shape.mass, leaf_temperature)

    # net specific metabolic heat → surface temperature
    specific_metabolic_heat_production = metabolic_heat_flow / o.body.geometry.volume
    (; surface_temperature) = surface_and_lung_temperature(;
        body=o.body,
        flesh_conductivity=internal_conduction.flesh_conductivity,
        specific_metabolic_heat_production,
        core_temperature=leaf_temperature,
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
        convection_out.mass,
        atmos,
        convection_area,
        surface_temperature,
        environment_vars.air_temperature;
        water_potential=hyd_pars.water_potential,
        gas_fractions=environment_pars.gas_fractions,
    )
    evaporation_heat_flow = evaporation_out.evaporation_heat_flow

    # heat balance (no respiration term for leaves — included in metabolic_heat_flow if needed)
    heat_balance_val = solar_flow + longwave_flow_in + metabolic_heat_flow -
                       longwave_flow_out - convection_heat_flow - evaporation_heat_flow

    energy_balance = (;
        solar_flow, longwave_flow_in, longwave_flow_out,
        convection_heat_flow, evaporation_heat_flow,
        metabolic_heat_flow, heat_balance=heat_balance_val,
    )
    oxygen_consumption_rate = u"ml/hr"(Joules_to_O2(Kleiber1961(), metabolic_heat_flow, 1.0))
    mass_balance = (; transpiration_mass=evaporation_out.transpiration_water_loss, oxygen_consumption_rate)

    return (;
        heat_balance=heat_balance_val,
        leaf_temperature,
        surface_temperature,
        energy_balance,
        mass_balance,
        evaporation_out,
        solar_out,
        longwave_gain_out,
        longwave_loss_out,
        convection_out,
    )
end
