# Coupled multi-part metabolic-rate solve for a single regulated compartment.
#
# This is the part-based reformulation of `solve_metabolic_rate`: the regulated
# compartment's core is held at the setpoint, each part's surface is solved at that
# core via `solve_part_surface`, and the metabolic heat flow that closes the
# whole-organism respiration balance is found by the same zbrent root-find used
# today. The difference from the dorsal/ventral path is that the internal heat that
# must be produced is the *plain sum* of genuine per-part `net_metabolic`, rather
# than the view-factor-weighted mean of two evaluations of one body.
#
# Floating compartments (legs, head — cores determined by `solve_core_temperatures`
# rather than pinned to the setpoint) are layered on top of this for general
# multi-compartment organisms; the single-`SharedCore`-compartment case here is the
# dorsal/ventral equivalence gate.

"""
    solve_coupled_metabolic_rate(; part_surface_setups, core_temperature,
        skin_temperature, insulation_temperature, temperature_tolerance,
        respire, respiration_pars, lung_mass, air_temperature, atmos, gas_fractions,
        metabolic_heat_flow_setpoint, resp_tolerance, smoothing = HardBound())

Solve for the metabolic heat flow of a single regulated compartment made of one or
more parts, all sharing the setpoint `core_temperature`.

`part_surface_setups` is a tuple of NamedTuples, one per part, each holding the
keyword inputs `solve_part_surface` needs *except* the temperatures and tolerance
(`body`, `insulation_pars`, `traits`, `environment_vars`, `conduction_fraction`,
`conductance_coefficient`, `ventral_fraction`, `longwave_depth_fraction`).

Returns a NamedTuple:
- `metabolic_heat_flow` — the closed metabolic heat rate (W)
- `parts` — per-part `solve_part_surface` results (skin, insulation, net_metabolic, …)
- `net_metabolic_total` — Σ per-part `net_metabolic` (the internal heat that must be produced)
- `skin_temperature`, `insulation_temperature` — part-mean surface temperatures
- `lung_temperature`
- `respiration_out` — the closing `respiration(...)` result, or `nothing` when `respire` is false
"""
function solve_coupled_metabolic_rate(;
    part_surface_setups::Tuple,
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
    smoothing::SmoothingStrategy=HardBound(),
)
    # 1. Per-part surface solve at the shared setpoint core.
    parts = map(part_surface_setups) do setup
        solve_part_surface(;
            setup...,
            skin_temperature,
            insulation_temperature,
            temperature_tolerance,
            smoothing,
        )
    end

    # 2. Internal heat that must be produced = plain sum of per-part core→skin flow.
    net_metabolic_total = sum(part -> part.net_metabolic, parts)

    # Part-mean surface temperatures and lung temperature.
    number_of_parts = length(parts)
    skin_mean = sum(part -> part.skin_temperature, parts) / number_of_parts
    insulation_mean = sum(part -> part.insulation_temperature, parts) / number_of_parts
    max_skin_temperature = maximum(part -> part.skin_temperature, parts)
    lung_temperature = (core_temperature + skin_mean) * 0.5

    # 3. Close the whole-organism respiration balance for metabolic_heat_flow.
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
        parts,
        net_metabolic_total,
        skin_temperature = skin_mean,
        insulation_temperature = insulation_mean,
        lung_temperature,
        respiration_out,
    )
end
