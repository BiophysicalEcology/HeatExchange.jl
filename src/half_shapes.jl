# ── Half shapes → parent-family thermal correlations ─────────────────────
#
# `Half{S}` is not a subtype of the family supertype `S <: AbstractShape`
# dispatches on (it wraps, not inherits), so it would miss the family methods
# below. In every one of these functions `shape` is a pure dispatch tag — all
# geometry comes from `body` (which already carries the halved dimensions). So a
# Half simply dispatches as its parent. The one half-sensitive quantity,
# `_shell_angle_fraction`, reads `body.shape` and so still sees the Half.

nusselt_free(h::Half, args...) = nusselt_free(h.parent, args...)
nusselt_forced(h::Half, args...) = nusselt_forced(h.parent, args...)
surface_and_lung_temperature(h::Half, args...) = surface_and_lung_temperature(h.parent, args...)
_core_to_skin_stack(h::Half, args...) = _core_to_skin_stack(h.parent, args...)
_layer_resistance(l, h::Half, args...; kw...) = _layer_resistance(l, h.parent, args...; kw...)
radiant_temperature(h::Half, args...; kw...) = radiant_temperature(h.parent, args...; kw...)
insulation_radiant_temperature(h::Half, args...; kw...) = insulation_radiant_temperature(h.parent, args...; kw...)
compressed_radiant_temperature(h::Half, args...; kw...) = compressed_radiant_temperature(h.parent, args...; kw...)
mean_skin_temperature(h::Half, args...; kw...) = mean_skin_temperature(h.parent, args...; kw...)
