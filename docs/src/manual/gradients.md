# Gradients, resistances and flows

The budget solved in [Solving a heat balance](heat_balance.md) is made of exchanges of heat, sources of heat
and, away from steady state, heat stored. With it go exchanges of mass: water and the gases of respiration.
Many of these exchanges share one form. Something flows because a potential differs between two places, and the
flow is that difference divided by a *resistance*, or multiplied by a *conductance*:

```math
\text{flow} = \frac{\text{potential difference}}{\text{resistance}} = \text{conductance} \times \text{potential difference}
```

Resistances here are those of a whole surface or layer, so area is part of the conductance: conduction is
``Q = (kA/L)\,\Delta T`` and convection ``Q = hA\,\Delta T``.

## Flows of heat

| Exchange | Potential | Conductance set by | Fit to the form |
|:--|:--|:--|:--|
| conduction through flesh, fat and fur | temperature | area, thickness and conductivity of each layer, see [Layers as a radial graph](radial_layers.md) | exact |
| conduction to the substrate | temperature | contact area, conductivity of the substrate | exact |
| conduction between parts of a body | temperature | area of the join and distance to it, see [Bodies of many parts](multipart.md) | exact |
| convection | temperature | area and boundary layer: size, shape, wind speed | exact, though the coefficient depends on the temperatures in free convection |
| longwave radiation | ``T^4`` | area, emissivity, view factor | exact in ``T^4`` |
| evaporation | vapour density | wet area, boundary layer, stomata | with the notes below |
| respiration | temperature and vapour density of the air breathed | volume of air breathed | partly: it is transport by moving air |

**Radiation.** Net exchange with one surface is ``Q = \varepsilon \sigma A F\,(T_s^4 - T_{env}^4)``. As a
difference of temperature the conductance contains the temperatures and is not constant. It is the *net*
exchange that vanishes at equal temperatures.

**Evaporation.** Water moves down its chemical potential, the water potential. Vapour density is the convenient
variable through the air, set at a surface by its temperature and water potential. Wetting more skin raises the
conductance, through the area. A drier or saltier surface lowers the potential.

**Respiration.** Air is carried in and out at a rate the organism sets. That is advection, not a passive
resistance. Blood carrying heat within the body is the same, and is represented by an effective conductivity
of the flesh.

## Flows of mass

The package is named for heat, but every solution also returns `mass_flows`: water evaporated, air breathed,
and the oxygen, carbon dioxide, nitrogen and water vapour entering and leaving, see [`MolarFluxes`](@ref).

| Flow | Potential | How it is computed here |
|:--|:--|:--|
| water vapour from skin, fur and leaves | vapour density | down the gradient, through the same boundary layer as heat |
| water vapour in the breath | vapour density of inhaled and exhaled air | by advection, exhaled at a set humidity and temperature |
| oxygen | partial pressure | from demand: what the metabolic rate requires |
| carbon dioxide | partial pressure | from the oxygen consumed and the respiratory quotient |
| air | pressure | from the oxygen required, its fraction in air and the fraction extracted, times any panting |

Only the first row is solved as a gradient. The partial-pressure gradients that move the gases are not resolved:
the oxygen an animal needs is assumed to arrive. That holds where oxygen is not limiting. In a burrow, at
altitude or in water it may not, and the air can be changed through `gas_fractions` in
[`EnvironmentalPars`](@ref).

Heat and mass are coupled at three places, each a fixed ratio between a flow of one and a flow of the other:

- **Latent heat.** Each gram of water evaporated carries about 2.4 kJ.
- **Heat of oxidation.** Each litre of oxygen consumed releases about 20 kJ, and carbon dioxide follows through
  the respiratory quotient, see [`O2_to_Joules`](@ref).
- **Ventilation.** The air that brings oxygen in takes heat and water out, so metabolic rate and respiratory
  loss are solved together. Panting moves more air than the oxygen demand requires.

So one metabolic rate implies a set of flows: heat, oxygen in, carbon dioxide out, air, and water in the breath.

The metabolic rate is itself, for now, an empirical function of mass and temperature, see
[Metabolism](metabolism.md), so the gas flows inherit its allometry. Dynamic Energy Budget theory 
(Kooijman 2010) replaces it with the chemistry: assimilation, maintenance, growth and reproduction, each with its 
own balance of mass. That couples the flows of food, water, oxygen, carbon dioxide and nitrogenous waste to each
other, and gives the heat and the entropy that the organism produces from first principles. It is to enter
through [AnimalMapper.jl](https://github.com/BiophysicalEcology/AnimalMapper.jl), by way of
[DEBtool_J.jl](https://github.com/add-my-pet/DEBtool_J.jl).

## Sources and storage

Solar radiation and metabolic heat are *sources*: no difference of temperature between organism and
surroundings drives them. Yet, they are not fixed; solar radiation depends on orientation, silhouette, 
absorptivity and shade, metabolism depends on activity, state and body temperature.

Heat *stored* is zero at steady state. Away from it, the residual of the budget is the rate of storage, see
[Solving a heat balance](heat_balance.md#Steady-state-and-storage).

If the air and every surface are at the organism's temperature, the air is saturated, and there is no sun and no
metabolic heat, every net exchange is zero. Remove one condition and something flows: wet skin in dry air at
body temperature still evaporates, and cools.

## The budget as a network

The budget is, therefore, a network: nodes with potentials, joined by conductances, with sources at some nodes. Steady
state is the set of potentials at which what flows into each node equals what flows out, see
[Layers as a radial graph](radial_layers.md). A residual is the net flow into a node, which is why it is the
natural thing to compute. Here, it is the heat gained less the heat lost: positive means the organism would warm.

Bond graphs (Paynter 1961) name the parts of such a network. Each *bond* carries an *effort*, the potential,
and a *flow*:

| Element | What it does | Here |
|:--|:--|:--|
| **0-junction** | one effort, flows sum to zero | a node: core, skin, fur surface. Its sum of flows is the residual |
| **1-junction** | one flow, efforts add | layers in series, see [`stack_resistance`](@ref) |
| **resistor** | relates flow to the effort across it | conduction, convection, net longwave, the boundary layer for water vapour: [`ConductiveShell`](@ref), [`ConductiveCoupling`](@ref) |
| **source of effort** | holds a potential for whatever flows | the surroundings: [`EnvironmentalVars`](@ref) |
| **source of flow** | supplies a flow, whatever the potential | absorbed solar radiation, metabolic heat: [`GeneratingCore`](@ref) |
| **capacitor** | stores what flows in, and its effort rises | heat capacity. Not yet in the package |
| **transformer** | turns one kind of flow into another at a fixed ratio | evaporation (water to heat), metabolism (oxygen to heat) |

The qualifications above are then differences of element, not exceptions. Sun and metabolism are sources of
flow, the surroundings sources of effort, and evaporation and metabolism the transformers that join the
networks of heat, water and gas. A transient budget is the same network with a capacitor at each node.

The correspondence is loose in three ways. Temperature with heat flow is a *pseudo* bond graph: in a true one
the flow is of entropy, so that effort times flow is power. The resistors are not linear. And a layer that
makes heat has its source spread through its resistance, see [Layers as a radial graph](radial_layers.md).

The package uses the vocabulary and does not build its equations from a bond graph.
[BondGraphs.jl](https://github.com/jedforrest/BondGraphs.jl) (Forrest et al. 2023) does, and comes from the
thermodynamics of biochemical networks (Gawthrop and Crampin 2014), the reactions that make metabolic heat.

## In the other packages

In [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) heat flows down gradients of
temperature, and water down gradients of water potential, through soil and air. The results are the boundary
conditions here: the temperature and vapour density of the air, and the temperatures of the ground and of the
sky, the last an effective radiative temperature.

## What an organism can change

An organism can act on each kind of term, see [Environments and the ecosystem](ecosystem.md).

- **A conductance**, by changing a parameter of its own. Raising fur thickens a layer. Blood sent to the skin is
  represented as more conductive flesh. Wetting skin enlarges the area that evaporates. Posture changes the areas
  exposed to air, sun, sky and ground, each differently.
- **A potential**, at either side of the organism/environment divide. Moving into shade, down a burrow or up a bush 
  replaces the environmental potentials. Letting its own body temperature change alters the potential on the 
  organism side.
- **A source, a sink or storage.** Colour and orientation change the solar radiation absorbed. Shivering and
  activity raise metabolic heat. Panting raises the air flow that carries heat and water away. A body temperature
  that varies transiently through the day (lagging behind the steady state) stores heat and releases it later.

Each is a change to the inputs of this package, followed by another solution.

## Gradients of information

The deciding factor for actions an organism may take, i.e. behaviors, is also a difference: between the state 
an organism senses and the state it expects or prefers. For example, a body temperature above the preferred 
one, a threshold level of dehydration, or, for an endotherm, a steady state heat budget implying a metabolic rate 
lower than it can produce. In control theory this is an error signal. Nothing flows down it. It is sensed, and acted 
upon.

It is a gradient in more than name. Under the free energy principle (Friston 2010) a living thing must stay
within the small set of states in which it remains what it is. A state outside that set is improbable for it,
and so surprising, and organisms act to keep surprise low. What they minimise is a bound on surprise, the
*free energy*: at its simplest a sum of squared prediction errors, each weighted by its precision (*not* physical 
energy flow of the organism). Action descends the gradient of free energy as heat descends the gradient of 
temperature.

Minimising that prediction error is *homeostasis* (Cannon 1932): the expected value is the set point, the
error the departure from it, the regulatory responses the actions that reduce it. A response made ahead of the
departure is *allostasis* (Sterling 2012). Thus, homeostasis is feedback, allostasis is feed-forward.

The two kinds of gradient run opposite ways. Physics takes heat and water down gradients of potential, towards
equilibrium with the surroundings. Homeostasis takes the organism down a gradient of prediction error, away
from that equilibrium, by rearranging the physical gradients. The first is the subject of this package and of
Microclimate.jl, the second of
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), see
[Gradients and control](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/gradients).

Two caveats should be borne in mind here. The error is not driven to zero: a response costs water, energy, time 
and safety, and stops when a further step is not worth it. And the free energy principle is a mathematical framing
from a physics perspective, not a mechanism that the controllers of BiophysicalBehaviour.jl implement.
