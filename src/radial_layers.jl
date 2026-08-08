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
# This is the production core→skin conduction path: `net_metabolic_heat` (net_metabolic_heat.jl)
# is a thin wrapper over `radial_net_metabolic_heat` here. The outer boundary is the separate
# per-part energy balance in `solve_part_heat_balance`; the ground-contact branch and the
# discretised radiative source remain design decisions — see docs/radial_layer_model.md.
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
between its inner and outer radii. Inner and outer radii are typed independently: for
some geometries they arrive in different (but dimension-compatible) unit representations.
"""
struct ConductiveShell{K,Ri,Ro} <: AbstractRadialLayer
    conductivity::K
    r_inner::Ri
    r_outer::Ro
end

# --- Per-layer thermal resistance (K/W), shape-dispatched --------------------
#
# These reproduce the per-term resistances inside `net_metabolic_heat` exactly. They
# are the standard radial-conduction results; `flesh_volume` supplies the length/area
# normalisation (`flesh_radius²/flesh_volume = 1/(π·length)` for a cylinder,
# `flesh_radius³/flesh_volume = 3/(4π)` for a sphere), so the same expression works for
# a shell at any radii, not only the flesh/fat boundary.

# The `smoothing` kwarg is used only by the ellipsoid (its safe_min semi-axis guard); the
# cylinder/sphere resistances accept and ignore it so `stack_resistance` can pass it uniformly.

# Cylinder / slab (grouped exactly as `net_metabolic_heat` groups them).
_layer_resistance(l::GeneratingCore, ::Union{AbstractCylindrical,AbstractSlab}, body; smoothing::SmoothingStrategy = HardBound()) =
    flesh_radius(body)^2 / (4 * l.conductivity * flesh_volume(body))
_layer_resistance(l::ConductiveShell, ::Union{AbstractCylindrical,AbstractSlab}, body; smoothing::SmoothingStrategy = HardBound()) =
    flesh_radius(body)^2 / (2 * l.conductivity * flesh_volume(body)) * log(l.r_outer / l.r_inner)

# Sphere.
_layer_resistance(l::GeneratingCore, ::AbstractSpherical, body; smoothing::SmoothingStrategy = HardBound()) =
    flesh_radius(body)^2 / (6 * l.conductivity * flesh_volume(body))
_layer_resistance(l::ConductiveShell, ::AbstractSpherical, body; smoothing::SmoothingStrategy = HardBound()) =
    flesh_radius(body)^3 / (3 * l.conductivity * flesh_volume(body)) *
    ((l.r_outer - l.r_inner) / (l.r_inner * l.r_outer))

# Ellipsoid. Concentric ellipsoidal shells have no exact separable conduction solution, so
# this uses the established equivalent-sphere approximation: an equivalent radius
# `r_eq = sqrt(3·ssqg)` (with `ssqg` the core's semi-axis combination), and the b-semi-minor
# axis as the radial coordinate. It telescopes across shells exactly like the sphere and
# reduces to the sphere when the body is spherical. `r_inner`/`r_outer` on the shell are the
# b-semi-minor axes of its inner/outer surfaces.
@inline function _ellipsoid_ssqg(body, smoothing::SmoothingStrategy)
    len = body.geometry.length
    a, b, c, fat = len.a_semi_major_skin, len.b_semi_minor_skin, len.c_semi_minor_skin, len.fat
    a2 = safe_min(smoothing, (a - fat)^2, a^2; scale = oneunit(a^2))
    b2 = safe_min(smoothing, (b - fat)^2, b^2; scale = oneunit(b^2))
    c2 = safe_min(smoothing, (c - fat)^2, c^2; scale = oneunit(c^2))
    return (a2 * b2 * c2) / (a2 * b2 + a2 * c2 + b2 * c2)
end

_layer_resistance(l::GeneratingCore, ::AbstractEllipsoidal, body; smoothing::SmoothingStrategy = HardBound()) =
    _ellipsoid_ssqg(body, smoothing) / (2 * l.conductivity * flesh_volume(body))
_layer_resistance(l::ConductiveShell, ::AbstractEllipsoidal, body; smoothing::SmoothingStrategy = HardBound()) =
    sqrt(3 * _ellipsoid_ssqg(body, smoothing))^3 / (3 * l.conductivity * flesh_volume(body)) *
    ((l.r_outer - l.r_inner) / (l.r_inner * l.r_outer))

"""
    stack_resistance(stack, body) -> resistance

Total series thermal resistance (K/W) of an ordered tuple of radial layers, summing
each layer's shape-dispatched resistance. Layers in series simply add.
"""
stack_resistance(stack::Tuple, body; smoothing::SmoothingStrategy = HardBound()) =
    sum(l -> _layer_resistance(l, shape(body), body; smoothing), stack)

"""
    core_to_skin_stack(body, conductivities) -> Tuple

The current model's core→skin radial stack: a generating flesh core in series with a
passive fat shell (`flesh_radius → skin_radius`). The fur shell (`skin_radius →
insulation_radius`) is the next segment out, handled at the surface node.
"""
core_to_skin_stack(body, c::ThermalConductivities; smoothing::SmoothingStrategy = HardBound()) =
    _core_to_skin_stack(shape(body), body, c, smoothing)

# Cylinder/sphere/slab use the actual flesh/skin radii as the shell's radial coordinates.
_core_to_skin_stack(::Union{AbstractCylindrical,AbstractSlab,AbstractSpherical}, body, c, ::SmoothingStrategy) = (
    GeneratingCore(c.flesh),
    ConductiveShell(c.fat, flesh_radius(body), skin_radius(body)),
)
# Ellipsoid uses the b-semi-minor axis of the flesh/skin surfaces as the shell coordinates
# (see `_layer_resistance` above). A zero-thickness fat shell contributes zero resistance.
function _core_to_skin_stack(::AbstractEllipsoidal, body, c, smoothing::SmoothingStrategy)
    len = body.geometry.length
    b_skin  = len.b_semi_minor_skin
    b_flesh = safe_min(smoothing, b_skin, b_skin - len.fat; scale = oneunit(b_skin))
    return (GeneratingCore(c.flesh), ConductiveShell(c.fat, b_flesh, b_skin))
end

"""
    radial_net_metabolic_heat(body, conductivities, core_temperature, skin_temperature)

Net metabolic heat conducted core→skin (W), computed as `(core − skin) / R` where `R`
is the series resistance of the explicit radial stack. Reproduces `net_metabolic_heat`
to machine precision for cylinders/slabs and spheres (gated), demonstrating that the
current closed form *is* this radial network collapsed.
"""
radial_net_metabolic_heat(body, c::ThermalConductivities, core_temperature, skin_temperature;
                          smoothing::SmoothingStrategy = HardBound()) =
    (core_temperature - skin_temperature) /
    stack_resistance(core_to_skin_stack(body, c; smoothing), body; smoothing)
