"""
    solve_coupled_metabolic_rate(body::CompositeBody, parts; neighbours = nothing, core_temperature,
        skin_temperature, insulation_temperature, temperature_tolerance, respire, respiration_pars,
        lung_mass, air_temperature, atmos, gas_fractions, metabolic_heat_flow_setpoint,
        resp_tolerance, extra_net_metabolic = 0 W, smoothing = HardBound())

The metabolic heat flow of `body` with every part's core held at `core_temperature`.

`parts` is keyed like `body.parts`, each entry holding the keywords of
[`solve_part_surface`](@ref) other than `body`, `covered_area` and the temperatures, which
come from `body` and the keywords here. `neighbours` is keyed like `body.parts`, each entry the
fractions of that part's view taken by every part, as the `neighbours` of
`silhouette_factors`; with them the parts exchange longwave radiation with each other.

`extra_net_metabolic` is heat the parts lose by other routes, such as conduction to other
compartments.

Returns `(; metabolic_heat_flow, parts, net_metabolic_total, skin_temperature,
insulation_temperature, lung_temperature, respiration_out)`, where `parts` holds each part's
`solve_part_surface` result and the temperatures are means over the parts.
"""
function solve_coupled_metabolic_rate(body::CompositeBody, parts::NamedTuple;
    neighbours=nothing,
    core_temperature,
    skin_temperature,
    insulation_temperature,
    temperature_tolerance,
    respire::Bool,
    respiration_pars,
    lung_mass,
    air_temperature,
    atmos,
    gas_fractions,
    metabolic_heat_flow_setpoint,
    resp_tolerance,
    extra_net_metabolic=zero(metabolic_heat_flow_setpoint),
    smoothing::SmoothingStrategy=HardBound(),
)
    covered = BiophysicalGeometry.covered_areas(body.parts, body.joins)
    setups = map((part, setup, covered_area) -> merge(setup, (; body = part, covered_area)), body.parts, parts, covered)
    solved = _solve_parts(setups, neighbours, skin_temperature, insulation_temperature, temperature_tolerance, smoothing)

    net_metabolic_total = sum(part -> part.net_metabolic, solved) + extra_net_metabolic
    skin_mean = sum(part -> part.skin_temperature, solved) / length(solved)
    insulation_mean = sum(part -> part.insulation_temperature, solved) / length(solved)
    max_skin_temperature = maximum(part -> part.skin_temperature, solved)
    lung_temperature = (core_temperature + skin_mean) * 0.5

    if respire
        minimum_flow = metabolic_heat_flow_setpoint
        flux_lower = metabolic_heat_flow_setpoint * (-2.0)
        flux_upper = max_skin_temperature >= core_temperature ?
            metabolic_heat_flow_setpoint * 1.01 : metabolic_heat_flow_setpoint * 10.0

        respiration_balance(metabolic_heat_flow) = respiration(
            MetabolicRates(; metabolic=metabolic_heat_flow * u"W", sum=net_metabolic_total, minimum=minimum_flow),
            respiration_pars, atmos, lung_mass, lung_temperature, air_temperature;
            gas_fractions, O2conversion=Kleiber1961(), smoothing,
        ).balance

        metabolic_heat_flow = zbrent(
            metabolic_heat_flow -> ustrip(u"W", respiration_balance(metabolic_heat_flow)),
            ustrip(u"W", flux_lower),
            ustrip(u"W", flux_upper),
            resp_tolerance * ustrip(u"W", metabolic_heat_flow_setpoint),
        ) * u"W"

        respiration_out = respiration(
            MetabolicRates(; metabolic=metabolic_heat_flow, sum=net_metabolic_total, minimum=minimum_flow),
            respiration_pars, atmos, lung_mass, lung_temperature, air_temperature;
            gas_fractions, O2conversion=Kleiber1961(), smoothing,
        )
        metabolic_heat_flow = respiration_out.metabolic_heat_flow
    else
        metabolic_heat_flow = net_metabolic_total
        respiration_out = nothing
    end

    return (;
        metabolic_heat_flow,
        parts = solved,
        net_metabolic_total,
        skin_temperature = skin_mean,
        insulation_temperature = insulation_mean,
        lung_temperature,
        respiration_out,
    )
end

_solve_parts(setups, ::Nothing, skin_temperature, insulation_temperature, temperature_tolerance, smoothing) =
    map(setup -> solve_part_surface(; setup..., skin_temperature, insulation_temperature, temperature_tolerance, smoothing), setups)

# Each part's surface sees its neighbours' surfaces, so the parts are solved again with
# their neighbours at the last surface temperatures until those stop changing.
function _solve_parts(setups, neighbours::NamedTuple, skin_temperature, insulation_temperature, temperature_tolerance, smoothing)
    temperatures = map(_ -> insulation_temperature, setups)
    solved = _solve_beside(setups, neighbours, temperatures, skin_temperature, insulation_temperature, temperature_tolerance, smoothing)
    for _ in 1:49
        temperatures = map(part -> part.insulation_temperature, solved)
        solved = _solve_beside(setups, neighbours, temperatures, skin_temperature, insulation_temperature, temperature_tolerance, smoothing)
        maximum(map((part, t) -> abs(part.insulation_temperature - t), solved, temperatures)) < temperature_tolerance && break
    end
    return solved
end

function _solve_beside(setups, neighbours, temperatures, skin_temperature, insulation_temperature, temperature_tolerance, smoothing)
    map(setups, neighbours) do setup, fractions
        beside = values(map((fraction, temperature) -> (; fraction, temperature), fractions, temperatures))
        environment_vars = merge(setup.environment_vars, (; neighbours = beside))
        solve_part_surface(; setup..., environment_vars, skin_temperature, insulation_temperature, temperature_tolerance, smoothing)
    end
end
