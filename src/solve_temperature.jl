# Passive thermoregulation solver (ectotherms, torpid endotherms, leaves)

"""
    solve_temperature(organism, environment; temperature_bracket=(270.0u"K", 370.0u"K"))

Find the steady-state temperature of an organism (animal or leaf) by root-finding on
`heat_balance`. Dispatches on `evaluation_strategy(organism)`:

- `SingleBody` (Naked animals, leaves): root-finds then returns the full `heat_balance` output.
- `MultiSided` (insulated animals): root-finds using a weighted per-side residual, then assembles
  a rich `(; thermoregulation, morphology, energy_flows, mass_flows)` output via
  `_assemble_multisided_output`.

Returns air temperature as the fallback `core_temperature` if root-finding fails.
"""
function solve_temperature(organism, environment; temperature_bracket=(270.0u"K", 370.0u"K"))
    solve_temperature(evaluation_strategy(organism), organism, environment; temperature_bracket)
end

function solve_temperature(::SingleBody, organism, environment; temperature_bracket=(270.0u"K", 370.0u"K"))
    lo = ustrip(u"K", temperature_bracket[1])
    hi = ustrip(u"K", temperature_bracket[2])
    fallback_temperature = environment.environment_vars.air_temperature
    try
        solved_temperature = zbrent(
            T -> ustrip(u"W", heat_balance(T * u"K", organism, environment).heat_balance),
            lo, hi, 1e-3,
        )
        heat_balance(solved_temperature * u"K", organism, environment)
    catch
        heat_balance(fallback_temperature, organism, environment)
    end
end

function solve_temperature(::MultiSided, organism, environment; temperature_bracket=(270.0u"K", 370.0u"K"))
    environment_pars = stripparams(environment.environment_pars)
    environment_vars = environment.environment_vars
    resp_pars  = respiration_pars(organism)
    metab_pars = metabolism_pars(organism)

    # Initial guesses: offset from setpoint (or midpoint of bracket)
    reference_temperature = metab_pars.core_temperature
    initial_skin_temperature = reference_temperature - 5u"K"
    initial_insulation_temperature  = reference_temperature - 3u"K"

    lo = ustrip(u"K", temperature_bracket[1])
    hi = ustrip(u"K", temperature_bracket[2])

    function residual(temperature)
        core_temperature = temperature * u"K"
        packed = _pack_sides(organism, environment, core_temperature, initial_skin_temperature, initial_insulation_temperature)
        dmult = packed.sky_factor_ref + packed.vegetation_factor_ref
        vmult = 1 - dmult
        net_metabolic = packed.temps_out[1].flows.net_metabolic * dmult +
                        packed.temps_out[2].flows.net_metabolic * vmult
        skin_mean = (packed.temps_out[1].skin_temperature + packed.temps_out[2].skin_temperature) * 0.5
        lung_temperature = (core_temperature + skin_mean) * 0.5
        metabolic_heat_flow = metabolic_rate(metab_pars.model, organism.body.shape.mass, core_temperature)
        resp_atmos = AtmosphericConditions(environment_vars)
        respiration_heat_flow = respiration(
            MetabolicRates(; metabolic=metabolic_heat_flow, sum=metabolic_heat_flow, minimum=metabolic_heat_flow),
            resp_pars, resp_atmos, organism.body.shape.mass, lung_temperature,
            environment_vars.air_temperature;
            gas_fractions=environment_pars.gas_fractions, O2conversion=Kleiber1961(),
        ).respiration_heat_flow
        ustrip(u"W", metabolic_heat_flow - respiration_heat_flow - net_metabolic)
    end

    equilibrium_temperature = try
        zbrent(residual, lo, hi, 1e-3) * u"K"
    catch
        environment_vars.air_temperature
    end

    # Re-run _pack_sides at converged temperature to assemble full output
    packed = _pack_sides(organism, environment, equilibrium_temperature, initial_skin_temperature, initial_insulation_temperature)
    dmult = packed.sky_factor_ref + packed.vegetation_factor_ref
    vmult = 1 - dmult
    skin_mean = packed.temps_out[1].skin_temperature * dmult + packed.temps_out[2].skin_temperature * vmult
    lung_temperature = (equilibrium_temperature + skin_mean) * 0.5
    metabolic_heat_flow = metabolic_rate(metab_pars.model, organism.body.shape.mass, equilibrium_temperature)
    resp_atmos = AtmosphericConditions(environment_vars)
    respiration_out = respiration(
        MetabolicRates(; metabolic=metabolic_heat_flow, sum=metabolic_heat_flow, minimum=metabolic_heat_flow),
        resp_pars, resp_atmos, organism.body.shape.mass, lung_temperature,
        environment_vars.air_temperature;
        gas_fractions=environment_pars.gas_fractions, O2conversion=Kleiber1961(),
    )
    metabolic_heat_flow = respiration_out.metabolic_heat_flow

    _assemble_multisided_output(organism, environment, equilibrium_temperature, metabolic_heat_flow, respiration_out, packed)
end
