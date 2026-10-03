# Layers as a radial graph

Heat produced at the core of an animal is conducted outwards through shells: flesh, fat, skin, fur. NicheMapR holds
this as a formula for each shape, with one layer of each kind. In HeatExchange.jl it is a chain of nodes joined by thermal
resistances, so that the number and kind of layers is data, not code. This page describes that chain, how far it
has been taken, and the plan for future versions of the package.

```@setup layers
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## The heat balance as a network

The heat balance of one part of a body is a network of temperature *nodes* joined by thermal resistances, with
heat entering or leaving at each node:

```@example layers
radial_network_diagram() # hide
```

The closed-form equations of the endotherm model are solutions to this network. In bond-graph terms the nodes
are 0-junctions, the layers are resistors in series, and the heat entering at a node is a source of flow, see
[Gradients, resistances and flows](gradients.md#The-budget-as-a-network).

**Nodes**

| Node | In the code | How it is found |
|:--|:--|:--|
| core | `core_temperature` | given, or solved for |
| boundary of flesh and fat | none | eliminated: the two resistances are added |
| skin | `skin_temperature` | solved for |
| outer surface of the fur | `insulation_temperature` | solved for |
| depth in the fur at which longwave radiation is exchanged | `radiant_temperature` | derived. The outer surface by default |
| compressed fur against the ground | `compressed_insulation_temperature` | derived, where the body touches the ground |

**Resistances**

| Between | Conductance, for a cylinder | Where |
|:--|:--|:--|
| core and boundary, through flesh that generates heat | ``4 k_{flesh} V / R^2`` | [`GeneratingCore`](@ref) |
| boundary and skin, through fat | ``2 k_{fat} V / (R^2 \ln(R_s / R))`` | [`ConductiveShell`](@ref) |
| skin and fur surface | ``2 \pi k_{fur} L / \ln(R_{fa} / R_s)`` | the surface solve, see [Insulation](insulation.md) |
| skin and substrate, through compressed fur | as for fur, with the compressed depth and conductivity, over the fraction in contact | the surface solve |
| fur surface and environment | convection, and longwave radiation in its linear form | [`solve_part_heat_balance`](@ref) |

**Heat entering and leaving**

| Node | Heat |
|:--|:--|
| core | metabolic heat generated, less the heat lost in breathing |
| flesh | the generation is spread through its volume, which gives the core its form of resistance |
| skin | less evaporation from the skin |
| fur surface | plus sunlight absorbed, less evaporation from wet fur, less convection and radiation |

From the core to the skin the network is a simple chain, the subject of the rest of this page. Beyond the skin
it is not, for three reasons given [at the end](#Where-the-chain-is-not-a-chain).

## Resistances in series

The heat that flows between two temperatures is their difference divided by a thermal resistance, in K/W, and
resistances in series add. For an animal with its core at ``T_c`` and its skin at ``T_s``:

```math
Q_{gen,net} = \frac{T_c - T_s}{R_{flesh} + R_{fat}}
```

The two resistances are of different kinds.

Fat is a passive shell. For a cylinder of length ``L``, a shell of conductivity ``k`` between radii ``r_{in}``
and ``r_{out}`` has

```math
R_{shell} = \frac{\ln(r_{out} / r_{in})}{2 \pi k L}
```

Flesh is not a plain resistor, because heat is produced throughout it. With heat generated evenly through a
cylinder of radius ``R`` and volume ``V``, the temperature falls as a parabola from the centre, and the
resistance from centre to surface is (Porter and Kearney 2009)

```math
R_{core} = \frac{R^2}{4 k V}
```

For a sphere the 4 becomes a 6.

## Layers

These are the two types of [`AbstractRadialLayer`](@ref):

| Type | Layer | Resistance |
|:--|:--|:--|
| [`GeneratingCore`](@ref)`(conductivity)` | flesh, making heat evenly through its volume | from the centre to its surface |
| [`ConductiveShell`](@ref)`(conductivity, r_inner, r_outer)` | fat, and any other passive layer | from its inner to its outer radius |

A *stack* is a tuple of layers from the centre outwards. [`stack_resistance`](@ref) adds their resistances, each
by a method for the family of the shape. [`core_to_skin_stack`](@ref) builds the stack of a body from its radii:

```@example layers
using HeatExchange, BiophysicalGeometry, Unitful
import HeatExchange: ThermalConductivities

fur = FibrousLayer(20.0u"mm", 30.0u"μm", 3000.0u"cm^-2")
fat = FatLayer(0.2, 901.0u"kg/m^3")
body = Body(Cylinder(10.0u"kg", 1000.0u"kg/m^3", 3.0), CompositeInsulation(fur, fat))
conductivities = ThermalConductivities(0.9u"W/m/K", 0.23u"W/m/K", nothing)   # flesh, fat, fur

stack = core_to_skin_stack(body, conductivities)
map(layer -> nameof(typeof(layer)), stack)
```

```@example layers
uconvert(u"K/W", stack_resistance(stack, body))
```

The body, with its fur and fat cut away:

```@example layers
shape_gallery("" => body; size = (420, 300)) # hide
```

The heat conducted from core to skin is the temperature difference over that resistance,
[`radial_net_metabolic_heat`](@ref). The solvers use it through [`net_metabolic_heat`](@ref):

```@example layers
core_temperature, skin_temperature = u"K"(37.0u"°C"), u"K"(32.0u"°C")
radial_net_metabolic_heat(body, conductivities, core_temperature, skin_temperature),
net_metabolic_heat(; body, conductivities, core_temperature, skin_temperature)
```

The stack of one core and one shell gives the same numbers as the closed forms of NicheMapR, which it replaced,
to rounding error for every shape. The tests of the package hold it to that.

## The temperature through the layers

With a resistance for each layer, the temperature at each boundary follows from the heat flowing through it.
Here the fat, though much the thinner layer, accounts for about two fifths of the difference between core and
skin:

```@example layers
flesh, fat_shell = stack
flow = radial_net_metabolic_heat(body, conductivities, core_temperature, skin_temperature)
boundary_temperature = core_temperature - flow * stack_resistance((flesh,), body)

surface_temperature = u"K"(15.0u"°C")   # of the fur, for the figure
radial_stack_diagram(["flesh", "fat", "fur"], # hide
    [flesh_radius(body), skin_radius(body), insulation_radius(body)], # hide
    [core_temperature, boundary_temperature, skin_temperature, surface_temperature]) # hide
```

Within the flesh the temperature follows a parabola, and within each shell a logarithm. The straight lines join
the values at the boundaries.

## Shapes

The resistance of a layer depends on the family of the shape:

| Family | Core | Shell |
|:--|:--|:--|
| cylinders and plates | ``R^2 / 4kV`` | ``\dfrac{R^2}{2kV} \ln \dfrac{r_{out}}{r_{in}}`` |
| spheres | ``R^2 / 6kV`` | ``\dfrac{R^3}{3kV} \dfrac{r_{out} - r_{in}}{r_{in} \, r_{out}}`` |
| ellipsoids | ``S^2 / 2kV`` | as a sphere of equivalent radius ``\sqrt{3 S^2}``, in the minor semi-axis |

```@example layers
layers = CompositeInsulation(fur, fat) # hide
layer_sections("Cylinder" => Body(Cylinder(10.0u"kg", 1000.0u"kg/m^3", 3.0), layers), # hide
    "Sphere" => Body(BiophysicalGeometry.Sphere(10.0u"kg", 1000.0u"kg/m^3"), layers), # hide
    "Ellipsoid" => Body(Ellipsoid(10.0u"kg", 1000.0u"kg/m^3", 3.0, 3.0), layers)) # hide
```

where ``R`` is the radius of the flesh, ``V`` its volume, and ``S^2 = a^2 b^2 c^2 / (a^2 b^2 + a^2 c^2 + b^2 c^2)``
for the semi-axes ``a``, ``b``, ``c`` of the flesh (Porter and Kearney 2009). Concentric ellipsoidal shells have
no exact solution of this kind, and the equivalent sphere is an approximation, exact for a sphere. The half
shapes of BiophysicalGeometry.jl use the methods of the whole shape.

## What a list of layers allows

A change in the structure of the animal is a change to the list:

- **No fat** is a shell of zero thickness, with zero resistance.
- **More layers** are more entries: blubber under muscle, clothing, snow on the back of an animal.
- **A large animal** can have its flesh divided into shells, each making heat, for a warm core inside cooler
  outer tissue.

```@example layers
clothed = (stack..., ConductiveShell(0.04u"W/m/K", skin_radius(body), skin_radius(body) + 5.0u"mm"))
uconvert(u"K/W", stack_resistance(clothed, body))
```

## Where this stands

The chain from core to skin is in use: every solver conducts heat through it. The chain from the skin through the fur to
the environment is not yet a list of layers. It is the surface solve of [`solve_part_heat_balance`](@ref), with
one layer of fur, see [Insulation](insulation.md).

| | State |
|:--|:--|
| One layer of flesh and one of fat, all shapes | in use, and reproduces NicheMapR |
| Any number of passive shells between core and skin | available through [`stack_resistance`](@ref); the solvers build the stack of flesh and fat |
| One layer of fur, by the surface solve | in use |
| Several shells of heat-generating flesh | planned. It needs the body to describe a boundary within the flesh |
| Several layers of fur or clothing | planned |
| Contact with the ground as a branch from the skin node | present in the surface solve, to become a branch of the graph |
| Longwave radiation exchanged at a depth within the fur | planned as a source at a node within the fur, in place of `longwave_depth_fraction` |
| Heat capacity at each node, for transients | planned, see [Solving a heat balance](heat_balance.md#Steady-state-and-storage) |

## Where the chain is not a chain

Three features of the heat balance beyond the skin are not resistances in series. Each was measured in the
current code before deciding what to do about it.

**The skin temperature is written as the mean of two estimates.** One comes from the core side: core
temperature less heat flow times the resistance of flesh and fat. The other from the environment side: fur
surface temperature plus heat flow times the resistance of the fur. `residual_skin_temperature` is the
difference between the skin temperature and their mean. In a network the skin node is fixed by one condition,
that heat arriving equals heat leaving. The two are the same: at a solution the estimates agree to about
10⁻¹³ K, for bodies from 0.1 to 100 kg, with and without ground contact and wet skin.

**Contact with the ground is a branch.** Heat leaves the skin by two paths side by side: through the free fur
to the air, and through the compressed fur to the substrate. It is not a small term. For a 1 kg cylinder under
a sky at −10 °C on a substrate at 27 °C:

| Fraction of the surface on the ground | Share of the heat exchanged, other than by evaporation, that passes through the ground |
|:--|:--|
| 0.1 | 8 % |
| 0.3 | 26 % |
| 0.5 | 45 % |

with the compressed fur a little under 1 K warmer than the free fur surface. So the structure must be a graph,
with a node for the substrate reached from the skin, not a list.

**Radiation can act within the fur.** `longwave_depth_fraction` places the exchange of longwave radiation at a
depth in the coat, and the closed form solves conduction and radiation at that depth together. Every example and
test uses a fraction of 1, the outer surface, where this coupling does nothing. For other values the closed
form divides by ``\ln(R_{fa} / R_{rad})``, which goes to zero as the depth approaches the surface, and for thin
fur the solution runs away. So radiation at depth is off by default and not usable when on. With the fur
divided into shells, it becomes heat entering at a node inside the fur, with nothing to divide by. The same
holds for sunlight, which in a deep or sparse coat is absorbed below the surface.

## Where this is going

The target is one small graph for each part of a body:

| | |
|:--|:--|
| nodes | core, shells of flesh, skin, shells of fur or clothing, surface |
| edges | a conductance for each shell, by the family of the shape, with a radiative part for fur |
| sources | metabolic heat in the flesh, sunlight and longwave radiation at nodes in the fur, evaporation at the skin and the surface |
| boundary | convection and radiation at the outer node, the one nonlinear node |
| branch | from the skin through compressed fur to the substrate |

Bare skin is then a graph with no shells of fur, a large animal one with several shells of flesh, and clothing
or snow more shells outside. The solution is that of a chain of linear conduction with one nonlinear boundary.
For an optimiser, each node temperature is a variable with one residual, which generalises the two variables
each part has now, see [Differentiability and the NLP interface](autodiff.md).

The same structure, nodes joined by conductances with a source at each node, is applied sideways between the
parts of a body in [Bodies of many parts](multipart.md).
