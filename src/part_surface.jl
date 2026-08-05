# Per-part surface solve — the side-agnostic replacement for one side of `_pack_sides`.
#
# Given a part's geometry, insulation, physiology, environment (view factors +
# boundaries), and an *imposed* compartment core temperature, converge the part's
# skin and insulation-surface temperatures against the same validated surface
# physics `solve_temperatures` runs today. Insulation and radiation exposure enter
# as independent per-part inputs (they are only coincidentally aligned in the
# dorsal/ventral special case): the part's own `insulation_pars` sets insulation,
# the part's own packed `environment_vars` (view factors, solar) sets radiation
# exposure, and `conduction_fraction`/`conductance_coefficient` set ground contact.
#
# `side = Dorsal()` is used uniformly: in the surface solve the side label only
# selects which stored insulation a part uses (via `_side_value`), and a single
# part carries one insulation (represented with equal dorsal/ventral fibres), so
# the label is immaterial. Ground contact is governed entirely by
# `conduction_fraction` + `conductance_coefficient`, not the side label.

"""
    solve_part_surface(; body, insulation_pars, traits, environment_vars,
                       conduction_fraction, conductance_coefficient,
                       ventral_fraction, longwave_depth_fraction,
                       skin_temperature, insulation_temperature,
                       temperature_tolerance, smoothing = HardBound())

Converge one part's skin and insulation-surface temperatures at the imposed core
temperature `traits.core_temperature`, then report the heat conducted from core to
skin and the effective core→skin conductance the compartment solve needs.

Reuses `solve_temperatures` — no physics is duplicated. Returns a NamedTuple:

- `skin_temperature`, `insulation_temperature` — converged surface temperatures
- `net_metabolic` — heat conducted core→skin through flesh/fat (W)
- `flesh_conductance` — effective `net_metabolic / (core − skin)` (W/K), the part's
  contribution to its compartment's diagonal flesh conductance
- `flows` — the full `HeatFlows` for this part
- `insulation_conductivity`, `tolerance`, `success`, `ntry` — solver diagnostics
"""
function solve_part_surface(;
    body::AbstractBody,
    insulation_pars::InsulationParameters,
    traits::NamedTuple,
    environment_vars::NamedTuple,
    conduction_fraction,
    conductance_coefficient,
    ventral_fraction,
    longwave_depth_fraction,
    skin_temperature,
    insulation_temperature,
    temperature_tolerance,
    covered_area=zero(BiophysicalGeometry.total_area(body)),
    characteristic_dim=characteristic_dimension(VolumeCubeRoot(), body),
    smoothing::SmoothingStrategy=HardBound(),
)
    core_temperature = traits.core_temperature
    insulation_temperature_mean = insulation_temperature * 0.7 + skin_temperature * 0.3
    insulation = insulation_properties(insulation_pars, insulation_temperature_mean, ventral_fraction; smoothing)
    geometry_vars = GeometryVariables(;
        side = Dorsal(),
        conductance_coefficient,
        ventral_fraction,
        conduction_fraction,
        longwave_depth_fraction,
    )
    # A joined part exposes only its uncovered surface: the flat face(s) that mate
    # with neighbouring parts (the SharedCore/conductive join patch) are internal,
    # so the convective/radiative/evaporative area is total − covered.
    geometry = (;
        total_area         = BiophysicalGeometry.total_area(body) - covered_area,
        area_evaporation   = evaporation_area(body) - covered_area,
        characteristic_dim,
    )
    result = solve_temperatures(;
        body,
        insulation_pars,
        insulation,
        geometry_vars,
        environment_vars,
        traits,
        temperature_tolerance,
        skin_temperature,
        insulation_temperature,
        geometry,
        smoothing,
    )
    net_metabolic = result.flows.net_metabolic
    flesh_conductance = _effective_flesh_conductance(net_metabolic, core_temperature, result.skin_temperature)
    return (;
        skin_temperature = result.skin_temperature,
        insulation_temperature = result.insulation_temperature,
        net_metabolic,
        flesh_conductance,
        flows = result.flows,
        insulation_conductivity = result.insulation_conductivity,
        tolerance = result.tolerance,
        success = result.success,
        ntry = result.ntry,
    )
end

# =============================================================================
# Non-iterative residual twin of `solve_part_surface` (multi-part NLP primitive).
#
# `solve_part_surface` *root-finds* a part's skin and insulation temperatures at
# an imposed core. Its twin here takes those temperatures (plus core, metabolic,
# and the per-part effectors) as *explicit inputs* and returns the heat-balance
# residuals — the form an NLP/IPOPT solver drives, where skin and insulation
# become decision variables instead of solver outputs. It rebuilds exactly the
# same geometry/insulation/geometry_vars `solve_part_surface` does (so the two
# agree residual-for-residual at a converged point) and calls the shared
# `solve_part_heat_balance` primitive — zero physics duplication.
#
# Two residuals per part are returned:
#   `surface_balance` = residual_energy_balance − residual_internal_conduction
#       — metabolic and respiration algebraically cancel, leaving pure surface
#       physics (solar − surface_losses + net_metabolic_internal = 0). This is
#       the same cancellation the dorsal/ventral `MultiSided` path uses.
#   `residual_skin_temperature` — skin_temperature minus the value implied by
#       the insulation-side heat balance.
# `net_metabolic_heat_internal` (flesh-conducted heat, independent of metabolic /
# respiration) is returned for the whole-organism balance the caller assembles
# once (Σ per-part internal heat = metabolic − respiration).
# =============================================================================

"""
    part_surface_residuals(setup, core_temperature, skin_temperature,
                           insulation_temperature, metabolic_heat_flow;
                           k_flesh, pant, skin_wetness, resp_pars,
                           smoothing = HardBound())

Evaluate one part's surface heat-balance residuals at explicit temperatures. The
`setup` is the same per-part NamedTuple `solve_part_surface` consumes (`body`,
`insulation_pars`, `traits`, `environment_vars`, `conduction_fraction`,
`conductance_coefficient`, `ventral_fraction`, `longwave_depth_fraction`,
`covered_area`, `characteristic_dim`).

Returns `(; surface_balance, residual_skin_temperature, net_metabolic_heat_internal, balance)`,
where `balance` is the full `solve_part_heat_balance` result.
"""
function part_surface_residuals(
    setup::NamedTuple,
    core_temperature,
    skin_temperature,
    insulation_temperature,
    metabolic_heat_flow;
    k_flesh,
    pant,
    skin_wetness,
    resp_pars,
    smoothing::SmoothingStrategy=HardBound(),
)
    (; body, insulation_pars, traits, environment_vars, conduction_fraction,
       conductance_coefficient, ventral_fraction, longwave_depth_fraction,
       covered_area, characteristic_dim) = setup

    insulation_temperature_mean = insulation_temperature * 0.7 + skin_temperature * 0.3
    insulation = insulation_properties(insulation_pars, insulation_temperature_mean, ventral_fraction; smoothing)
    geometry_vars = GeometryVariables(;
        side = Dorsal(),
        conductance_coefficient,
        ventral_fraction,
        conduction_fraction,
        longwave_depth_fraction,
    )
    geometry = (;
        total_area         = BiophysicalGeometry.total_area(body) - covered_area,
        area_evaporation   = evaporation_area(body) - covered_area,
        characteristic_dim,
    )
    balance = solve_part_heat_balance(
        core_temperature, skin_temperature, insulation_temperature, metabolic_heat_flow;
        body,
        geometry,
        insulation_pars,
        insulation,
        geometry_vars,
        environment_vars,
        traits,
        resp_pars,
        k_flesh,
        pant,
        skin_wetness,
        smoothing,
    )
    surface_balance = balance.residual_energy_balance - balance.residual_internal_conduction
    return (;
        surface_balance,
        residual_skin_temperature = balance.residual_skin_temperature,
        net_metabolic_heat_internal = balance.net_metabolic_heat_internal,
        balance,
    )
end

# Effective core→skin conductance G_flesh = net_metabolic / (core − skin). The
# core-skin gradient is strictly positive for a heat-generating endotherm; guard a
# vanishing gradient so the compartment diagonal never sees Inf/NaN. The floor is
# well below any physical resolution, so it never perturbs a real solve.
function _effective_flesh_conductance(net_metabolic, core_temperature, skin_temperature)
    gradient = core_temperature - skin_temperature
    floor_gradient = 1e-6u"K"
    safe_gradient = abs(gradient) < floor_gradient ? oftype(gradient, floor_gradient) : gradient
    return net_metabolic / safe_gradient
end
