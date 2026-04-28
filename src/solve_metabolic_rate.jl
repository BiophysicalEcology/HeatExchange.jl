# Active thermoregulation solver (endotherms, thermogenic plants)
# Finds the metabolic_heat_flow that closes the heat budget at a fixed core temperature.

"""
    SolveMetabolicRateOptions

Options controlling the endotherm metabolic rate solver.

# Fields
- `respire`: Whether to include respiration in the heat balance (default: true)
- `temperature_error_tolerance`: Convergence tolerance for simultaneous temperature solution (default: 1e-3 K)
- `resp_tolerance`: Relative convergence tolerance for respiration root finding (default: 1e-5)
"""
Base.@kwdef struct SolveMetabolicRateOptions{RE,ST,BT} <: AbstractModelParameters
    respire::RE = Param(true)
    temperature_error_tolerance::ST = Param(1e-3u"K")
    resp_tolerance::BT = Param(1e-5)
end

"""
    solve_metabolic_rate(o::Organism, e, skin_temperature, insulation_temperature)

Solve for the metabolic rate that balances heat production with heat loss for an endotherm.

Finds `metabolic_heat_flow` at `metab_pars.core_temperature` by:
1. Solving skin and insulation surface temperatures iteratively per body side (via `_pack_sides`).
2. Root-finding (zbrent) on the respiration heat balance.

# Arguments
- `o::Organism`: Organism with body geometry and traits
- `e`: Environment containing `environment_pars` and `environment_vars`
- `skin_temperature`: Initial guess for skin temperature
- `insulation_temperature`: Initial guess for insulation surface temperature

# Returns
NamedTuple with:
- `thermoregulation`: Temperature and conductivity outputs
- `morphology`: Body geometry and areas
- `energy_flows`: Heat flux components (solar, longwave, convection, etc.)
- `mass_flows`: Water and gas exchange rates
"""
function solve_metabolic_rate(o::Organism, e, skin_temperature, insulation_temperature)
    environment_pars = stripparams(e.environment_pars)
    environment_vars = e.environment_vars
    opts = options(o)
    resp_pars = respiration_pars(o)
    metab_pars = metabolism_pars(o)

    packed = _pack_sides(o, e, metab_pars.core_temperature, skin_temperature, insulation_temperature)
    (; temps_out, side_bodies, sky_factor_ref, vegetation_factor_ref) = packed

    max_skin_temperature = max(temps_out[1].skin_temperature, temps_out[2].skin_temperature)

    gen_d = temps_out[1].flows.net_metabolic
    gen_v = temps_out[2].flows.net_metabolic
    dmult = sky_factor_ref + vegetation_factor_ref
    vmult = 1 - dmult
    flow_sum = gen_d * dmult + gen_v * vmult

    skin_temperature = (temps_out[1].skin_temperature + temps_out[2].skin_temperature) * 0.5
    insulation_temperature = (temps_out[1].insulation_temperature + temps_out[2].insulation_temperature) * 0.5
    lung_temperature = (metab_pars.core_temperature + skin_temperature) * 0.5

    if opts.respire
        minimum_flow = metab_pars.metabolic_heat_flow
        flux_m1 = metab_pars.metabolic_heat_flow * (-2.0)
        flux_m2 = metab_pars.metabolic_heat_flow * 10.0
        if max_skin_temperature >= metab_pars.core_temperature
            flux_m2 = metab_pars.metabolic_heat_flow * 1.01
        end
        resp_atmos = AtmosphericConditions(environment_vars)
        f = x -> ustrip(u"W", respiration(
            MetabolicRates(; metabolic=x * u"W", sum=flow_sum, minimum=minimum_flow),
            resp_pars,
            resp_atmos,
            side_bodies[2].shape.mass,
            lung_temperature,
            environment_vars.air_temperature;
            gas_fractions=environment_pars.gas_fractions,
            O2conversion=Kleiber1961(),
        ).balance)

        metabolic_heat_flow = zbrent(
            f,
            ustrip(u"W", flux_m1),
            ustrip(u"W", flux_m2),
            opts.resp_tolerance * ustrip(u"W", metab_pars.metabolic_heat_flow),
        ) * u"W"

        respiration_out = respiration(
            MetabolicRates(; metabolic=metabolic_heat_flow, sum=flow_sum, minimum=minimum_flow),
            resp_pars,
            resp_atmos,
            side_bodies[2].shape.mass,
            lung_temperature,
            environment_vars.air_temperature;
            gas_fractions=environment_pars.gas_fractions,
            O2conversion=Kleiber1961(),
        )
        metabolic_heat_flow = respiration_out.metabolic_heat_flow
    else
        metabolic_heat_flow = flow_sum
        minimum_flow = metab_pars.metabolic_heat_flow
    end

    resp_out_for_assembly = opts.respire ? respiration_out : nothing
    return _assemble_multisided_output(o, e, metab_pars.core_temperature, metabolic_heat_flow, resp_out_for_assembly, packed)
end
