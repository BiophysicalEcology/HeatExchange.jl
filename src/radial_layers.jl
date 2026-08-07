# =============================================================================
# Radial conduction as an ordered stack of concentric layers.
#
# Design + rationale: docs/radial_layer_model.md. The current core→skin physics
# (`net_metabolic_heat`) is a radial conduction network already solved in closed
# form — a uniformly heat-generating flesh core in series with passive conductive
# shells (fat, then fur). This module makes that network explicit: an ordered list of
# layers, each contributing a shape-dispatched thermal resistance, whose series sum
# reproduces `net_metabolic_heat` exactly (gated in test/radial_layers.jl). Adding a
# layer — a second flesh shell for a large animal, extra fur, clothing — becomes an
# entry in the list, not a new bespoke formula.
#
# Phase 8 foundation: the conduction chain only (core→skin, reproduce-gated). The
# outer boundary node (convection/radiation/evaporation), the ground-contact branch,
# and the discretised radiative source are the next phases — see the design doc's
# "target shape" and "open design decisions" sections. Nothing here is wired into the
# production heat balance yet; it is the proven foundation the redesign builds on.
# =============================================================================

abstract type AbstractRadialLayer end

"""
    GeneratingCore(conductivity)

The innermost region (flesh): a solid body with uniform volumetric heat generation.
Its centre→surface resistance is the analytic Poisson result (`r²/4kV` for a cylinder,
`r²/6kV` for a sphere), which is what makes it *not* a plain resistor — the metabolic
heat is produced throughout its volume, so the temperature profile is parabolic.
"""
struct GeneratingCore{K} <: AbstractRadialLayer
    conductivity::K
end

"""
    ConductiveShell(conductivity, r_inner, r_outer)

A passive concentric shell (fat, fur, clothing, snow) carrying heat by conduction
between its inner and outer radii.
"""
struct ConductiveShell{K,R} <: AbstractRadialLayer
    conductivity::K
    r_inner::R
    r_outer::R
end

# --- Per-layer thermal resistance (K/W), shape-dispatched --------------------
#
# These reproduce the per-term resistances inside `net_metabolic_heat` exactly. They
# are the standard radial-conduction results; `flesh_volume` supplies the length/area
# normalisation (`flesh_radius²/flesh_volume = 1/(π·length)` for a cylinder,
# `flesh_radius³/flesh_volume = 3/(4π)` for a sphere), so the same expression works for
# a shell at any radii, not only the flesh/fat boundary.

# Cylinder / slab (grouped exactly as `net_metabolic_heat` groups them).
_layer_resistance(l::GeneratingCore, ::Union{AbstractCylindrical,AbstractSlab}, body) =
    flesh_radius(body)^2 / (4 * l.conductivity * flesh_volume(body))
_layer_resistance(l::ConductiveShell, ::Union{AbstractCylindrical,AbstractSlab}, body) =
    flesh_radius(body)^2 / (2 * l.conductivity * flesh_volume(body)) * log(l.r_outer / l.r_inner)

# Sphere.
_layer_resistance(l::GeneratingCore, ::AbstractSpherical, body) =
    flesh_radius(body)^2 / (6 * l.conductivity * flesh_volume(body))
_layer_resistance(l::ConductiveShell, ::AbstractSpherical, body) =
    flesh_radius(body)^3 / (3 * l.conductivity * flesh_volume(body)) *
    ((l.r_outer - l.r_inner) / (l.r_inner * l.r_outer))

# Ellipsoid is genuinely anisotropic (its "shells" are parametrised by three semi-axes,
# not a single radius), so it is not a clean instance of the r_inner/r_outer shell yet.
# Deferred — see the design doc; `net_metabolic_heat` still handles it directly.

"""
    stack_resistance(stack, body) -> resistance

Total series thermal resistance (K/W) of an ordered tuple of radial layers, summing
each layer's shape-dispatched resistance. Layers in series simply add.
"""
stack_resistance(stack::Tuple, body) =
    sum(l -> _layer_resistance(l, shape(body), body), stack)

"""
    core_to_skin_stack(body, conductivities) -> Tuple

The current model's core→skin radial stack: a generating flesh core in series with a
passive fat shell (`flesh_radius → skin_radius`). The fur shell (`skin_radius →
insulation_radius`) is the next segment out, handled at the surface node.
"""
core_to_skin_stack(body, c::ThermalConductivities) = (
    GeneratingCore(c.flesh),
    ConductiveShell(c.fat, flesh_radius(body), skin_radius(body)),
)

"""
    radial_net_metabolic_heat(body, conductivities, core_temperature, skin_temperature)

Net metabolic heat conducted core→skin (W), computed as `(core − skin) / R` where `R`
is the series resistance of the explicit radial stack. Reproduces `net_metabolic_heat`
to machine precision for cylinders/slabs and spheres (gated), demonstrating that the
current closed form *is* this radial network collapsed.
"""
radial_net_metabolic_heat(body, c::ThermalConductivities, core_temperature, skin_temperature) =
    (core_temperature - skin_temperature) / stack_resistance(core_to_skin_stack(body, c), body)
