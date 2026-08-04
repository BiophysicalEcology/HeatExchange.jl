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
