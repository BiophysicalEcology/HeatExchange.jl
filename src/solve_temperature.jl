# Passive thermoregulation solver (ectotherms, torpid endotherms, leaves)

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
    solve_temperature(organism, environment; T_bracket=(270.0u"K", 370.0u"K"))

Find the steady-state temperature of an organism (animal or leaf) by root-finding on
`heat_balance`. Dispatches on `evaluation_strategy(organism)`:

- `SingleBody` (Naked animals, leaves): root-finds then returns the full `heat_balance` output.
- `MultiSided` (insulated animals): root-finds using a weighted per-side residual, then assembles
  a rich `(; thermoregulation, morphology, energy_flows, mass_flows)` output via
  `_assemble_multisided_output`.

Returns air temperature as the fallback `core_temperature` if root-finding fails.
"""
function solve_temperature(organism, environment; T_bracket=(270.0u"K", 370.0u"K"))
    solve_temperature(evaluation_strategy(organism), organism, environment; T_bracket)
end

function solve_temperature(::SingleBody, organism, environment; T_bracket=(270.0u"K", 370.0u"K"))
    lo = ustrip(u"K", T_bracket[1])
    hi = ustrip(u"K", T_bracket[2])
    T_fallback = environment.environment_vars.air_temperature
    try
        T_sol = zbrent(
            T -> ustrip(u"W", heat_balance(T * u"K", organism, environment).heat_balance),
            lo, hi, 1e-3,
        )
        heat_balance(T_sol * u"K", organism, environment)
    catch
        heat_balance(T_fallback, organism, environment)
    end
end

function solve_temperature(::MultiSided, organism, environment; T_bracket=(270.0u"K", 370.0u"K"))
    environment_pars = stripparams(environment.environment_pars)
    environment_vars = environment.environment_vars
    resp_pars  = respiration_pars(organism)
    metab_pars = metabolism_pars(organism)

    # Initial guesses: offset from setpoint (or midpoint of bracket)
    T_ref = metab_pars.core_temperature
    T_skin_init = T_ref - 5u"K"
    T_ins_init  = T_ref - 3u"K"

    lo = ustrip(u"K", T_bracket[1])
    hi = ustrip(u"K", T_bracket[2])

    function residual(T_val)
        T_core = T_val * u"K"
        packed = _pack_sides(organism, environment, T_core, T_skin_init, T_ins_init)
        dmult = packed.sky_factor_ref + packed.vegetation_factor_ref
        vmult = 1 - dmult
        net_generated = packed.temps_out[1].flows.net_generated * dmult +
                        packed.temps_out[2].flows.net_generated * vmult
        skin_mean = (packed.temps_out[1].skin_temperature + packed.temps_out[2].skin_temperature) * 0.5
        lung_temperature = (T_core + skin_mean) * 0.5
        Q_gen = metabolic_rate(metab_pars.model, organism.body.shape.mass, T_core)
        resp_atmos = AtmosphericConditions(environment_vars)
        Q_resp = respiration(
            MetabolicRates(; metabolic=Q_gen, sum=Q_gen, minimum=Q_gen),
            resp_pars, resp_atmos, organism.body.shape.mass, lung_temperature,
            environment_vars.air_temperature;
            gas_fractions=environment_pars.gas_fractions, O2conversion=Kleiber1961(),
        ).respiration_heat_flow
        ustrip(u"W", Q_gen - Q_resp - net_generated)
    end

    T_eq = try
        zbrent(residual, lo, hi, 1e-3) * u"K"
    catch
        environment_vars.air_temperature
    end

    # Re-run _pack_sides at converged temperature to assemble full output
    packed = _pack_sides(organism, environment, T_eq, T_skin_init, T_ins_init)
    dmult = packed.sky_factor_ref + packed.vegetation_factor_ref
    vmult = 1 - dmult
    skin_mean = packed.temps_out[1].skin_temperature * dmult + packed.temps_out[2].skin_temperature * vmult
    lung_temperature = (T_eq + skin_mean) * 0.5
    Q_gen = metabolic_rate(metab_pars.model, organism.body.shape.mass, T_eq)
    resp_atmos = AtmosphericConditions(environment_vars)
    respiration_out = respiration(
        MetabolicRates(; metabolic=Q_gen, sum=Q_gen, minimum=Q_gen),
        resp_pars, resp_atmos, organism.body.shape.mass, lung_temperature,
        environment_vars.air_temperature;
        gas_fractions=environment_pars.gas_fractions, O2conversion=Kleiber1961(),
    )
    generated_heat_flow = respiration_out.generated_heat_flow

    _assemble_multisided_output(organism, environment, T_eq, generated_heat_flow, respiration_out, packed)
end
