# Gradients, resistances and flows

The heat budget of an organism, set out in [Solving a heat balance](heat_balance.md), is a sum of flows of heat.
Almost every one of them has one form. Something flows because there is a difference in a *potential*
between two places, and the size of the flow is that difference divided by a *resistance*, or multiplied by a
conductance:

```math
\text{flow} = \frac{\text{potential difference}}{\text{resistance}}
```

| Flow | Potential | Resistance set by |
|:--|:--|:--|
| conduction through flesh, fat and fur | temperature | the thickness and thermal conductivity of each layer, see [Layers as a radial graph](radial_layers.md) |
| conduction to the substrate | temperature | the area of contact and the conductivity of the substrate |
| conduction between the parts of a body | temperature | the area of the join and the distance to it, see [Bodies of many parts](multipart.md) |
| convection | temperature | the boundary layer of air: size, shape and wind speed |
| longwave radiation | the fourth power of temperature | emissivity and view factor. In its linear form, a temperature difference over a radiative resistance |
| evaporation from skin, fur and leaves | vapour density, set at a surface by its temperature and water potential | the boundary layer, the wet fraction of the skin, the stomata |
| respiration | vapour density and temperature of the air breathed in and out | the volume of air breathed, which follows the metabolic rate |

Solar radiation and metabolic heat are the exceptions. They are sources: heat that arrives at a node of the
system whatever its temperature. Everything else is a flow down a gradient, and is zero when the gradient is zero.
An organism at the temperature of its surroundings, with no sun and no metabolism, exchanges nothing.

## The budget as a network

Seen this way the heat budget is a network: nodes that each have a potential, joined by resistances, with
sources at some of the nodes. The steady state is the set of potentials at which what flows into each node equals
what flows out. That is the picture drawn in [Layers as a radial graph](radial_layers.md) for one part of a body,
and it is the reason that a residual is the natural thing to compute: a residual is the net flow into a node.

## In the other packages

The same form runs through the other packages. In
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) heat flows down gradients of temperature
through the soil and the air above it, and water flows down gradients of water potential through the soil, each
against a resistance set by the properties of the soil and the turbulence of the air. The potentials at the
surface of the organism in this package, the temperature and vapour density of the air and the temperatures of
the ground and sky, are the output of those flows.

## What an organism can change

An organism has two ways to change a flow, and they are the two terms of the equation.

It can change a **resistance**, by changing one of its own parameters. Raising the fur thickens a layer. Sending
blood to the skin raises the conductivity of the flesh. Curling up reduces the area through which every flow
passes. Wetting the skin lowers the resistance to evaporation. Changing colour changes how much of a source is
absorbed.

It can change a **gradient**, by changing where it is. Moving into shade, down a burrow or up a bush replaces
the potentials around it with others. Letting its own temperature rise changes the gradient from the other end.

These are the responses listed in [Environments and the ecosystem](ecosystem.md), and each is a change to the
inputs of this package followed by another solution.

## Gradients of information

What decides which change is made is also a gradient, of a different kind. An organism responds to the
difference between its present state and a target state: a body temperature above the preferred one, or a
metabolic rate below the minimum that it can produce. That difference drives its behaviour as a temperature
difference drives a flow of heat, and the response continues until the difference is gone or the means are
exhausted. But it is a gradient of *information*. Nothing flows down it. It is sensed, and acted upon. The
physical gradients decide what will happen to an organism in a given state and place, and are the subject of this
package and of Microclimate.jl. The gradient between state and target decides what the organism does about it,
and is the subject of
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), see
[Gradients and control](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/gradients) and [Behaviour as control](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/control) in its documentation.
