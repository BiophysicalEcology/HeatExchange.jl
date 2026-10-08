# A half shape takes its parent's correlations: in these the shape only selects the
# method, and the geometry comes from the body.

nusselt_free(h::Half, args...) = nusselt_free(h.parent, args...)
nusselt_forced(h::Half, args...) = nusselt_forced(h.parent, args...)
surface_and_lung_temperature(h::Half, args...) = surface_and_lung_temperature(h.parent, args...)
_core_to_skin_stack(h::Half, args...) = _core_to_skin_stack(h.parent, args...)
_layer_resistance(l, h::Half, args...; kw...) = _layer_resistance(l, h.parent, args...; kw...)
radiant_temperature(h::Half, args...; kw...) = radiant_temperature(h.parent, args...; kw...)
mean_skin_temperature(h::Half, args...; kw...) = mean_skin_temperature(h.parent, args...; kw...)
