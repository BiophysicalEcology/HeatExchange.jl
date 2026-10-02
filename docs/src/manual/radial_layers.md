# Layers as a radial graph

Heat made in the core of an animal is conducted outwards through a series of shells: flesh, then fat, then skin,
then fur. NicheMapR holds this as a formula for each shape, with one layer of each kind. Here it is held as what it
is, a chain of nodes joined by thermal resistances, so that the number and kind of layers is data and not code.
This page describes that chain, how far it has been taken, and where it is going.

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

The closed-form equations of the endotherm model are this network already solved. Writing it out shows what is
in it.

**Nodes**

| Node | In the code | How it is found |
|:--|:--|:--|
| core | `core_temperature` | given, or solved for |
| boundary of flesh and fat | none | eliminated: the two resistances are added |
| skin | `skin_temperature` | solved for |
| outer surface of the fur | `insulation_temperature` | solved for |
| depth in the fur at which longwave radiation is exchanged | `radiant_temperature` | derived. It is the outer surface by default |
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
| flesh | the generation is spread through its volume, which is what gives the core its form of resistance |
| skin | less evaporation from the skin |
| fur surface | plus sunlight absorbed, less evaporation from wet fur, less convection and radiation |

The part of this network from the core to the skin is a simple chain, and is the subject of the rest of this
page. The part beyond the skin is not, for the three reasons given [at the end](#Where-the-chain-is-not-a-chain).

## Resistances in series

By analogy with Ohm's law, the heat that flows between two temperatures is their difference divided by a thermal
resistance, in K/W, and resistances in series add. Temperature is the potential and heat is what flows, see
[Gradients, resistances and flows](gradients.md). For an animal with its core at ``T_c`` and its skin at ``T_s``:

```math
Q_{gen,net} = \frac{T_c - T_s}{R_{flesh} + R_{fat}}
```

The two resistances are of different kinds.

Fat is a passive shell. For a cylinder of length ``L``, the resistance of a shell of conductivity ``k`` between
radii ``r_{in}`` and ``r_{out}`` is

```math
R_{shell} = \frac{\ln(r_{out} / r_{in})}{2 \pi k L}
```

Flesh is not a plain resistor, because the heat is produced throughout it. With heat generated evenly through a
cylinder of radius ``R`` and volume ``V``, the temperature falls as a parabola from the centre, and the resistance
from the centre to the surface is (Porter and Kearney 2009)

```math
R_{core} = \frac{R^2}{4 k V}
```

For a sphere the 4 becomes a 6.

## Layers

These two kinds are the two types of [`AbstractRadialLayer`](@ref):

| Type | Layer | Resistance |
|:--|:--|:--|
| [`GeneratingCore`](@ref)`(conductivity)` | flesh, making heat evenly through its volume | from the centre to its surface |
| [`ConductiveShell`](@ref)`(conductivity, r_inner, r_outer)` | fat, and any other passive layer | from its inner to its outer radius |

A *stack* is a tuple of layers from the centre outwards, and [`stack_resistance`](@ref) adds their resistances,
each found by a method for the family of the shape. [`core_to_skin_stack`](@ref) builds the stack of a body from
its radii:

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

The heat conducted from the core to the skin is then the temperature difference over that resistance,
[`radial_net_metabolic_heat`](@ref), and this is the function that the solvers use, through
[`net_metabolic_heat`](@ref):

```@example layers
core_temperature, skin_temperature = u"K"(37.0u"°C"), u"K"(32.0u"°C")
radial_net_metabolic_heat(body, conductivities, core_temperature, skin_temperature),
net_metabolic_heat(; body, conductivities, core_temperature, skin_temperature)
```

The stack of one core and one shell gives the same numbers as the closed forms of NicheMapR, which it replaced, to
within rounding error for every shape. The tests of the package hold it to that.

## The temperature through the layers

With a resistance for each layer, the temperature at each boundary follows from the heat flowing through it. Here
the fat, though much the thinner layer, accounts for about two fifths of the difference between core and skin:

```@example layers
flesh, fat_shell = stack
flow = radial_net_metabolic_heat(body, conductivities, core_temperature, skin_temperature)
boundary_temperature = core_temperature - flow * stack_resistance((flesh,), body)

surface_temperature = u"K"(15.0u"°C")   # of the fur, for the figure
radial_stack_diagram(["flesh", "fat", "fur"], # hide
    [flesh_radius(body), skin_radius(body), insulation_radius(body)], # hide
    [core_temperature, boundary_temperature, skin_temperature, surface_temperature]) # hide
```

Within the flesh the temperature in fact follows a parabola, and within each shell a logarithm. The straight
lines join the values at the boundaries.

## Shapes

The resistance of a layer depends on the family of the shape of the body:

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
for the semi-axes ``a``, ``b``, ``c`` of the flesh (Porter and Kearney 2009). Concentric ellipsoidal shells have no
exact solution of this kind, and the equivalent sphere is an approximation that is exact for a sphere. The half
shapes of BiophysicalGeometry.jl use the methods of the whole shape.

## What a list of layers allows

Because the layers are a list, a change in the structure of the animal is a change to the list:

- **No fat** is a shell of zero thickness, with zero resistance.
- **More layers** are more entries. A second shell of a different conductivity can stand for a layer of blubber
  under a layer of muscle, or for clothing, or for snow on the back of an animal.
- **A large animal** can have its flesh divided into shells, each making heat, to represent a warm core inside
  cooler outer tissue.

```@example layers
clothed = (stack..., ConductiveShell(0.04u"W/m/K", skin_radius(body), skin_radius(body) + 5.0u"mm"))
uconvert(u"K/W", stack_resistance(clothed, body))
```

## Where this stands

The chain from the core to the skin is in use: every solver in the package conducts heat through it. The rest of
the path, from the skin through the fur to the environment, is not yet a list of layers. It is the surface solve of
[`solve_part_heat_balance`](@ref), which has one layer of fur, see [Insulation](insulation.md). So, in this
version:

| | State |
|:--|:--|
| One layer of flesh and one of fat, all shapes | in use, and reproduces NicheMapR |
| Any number of passive shells between the core and the skin | available through [`stack_resistance`](@ref); the solvers build the stack of flesh and fat |
| One layer of fur, by the surface solve | in use |
| Several shells of heat-generating flesh | planned. It needs the body to describe a boundary within the flesh |
| Several layers of fur or clothing | planned |
| Contact with the ground as a branch from the skin node | present in the surface solve, to become a branch of the graph |
| Longwave radiation exchanged at a depth within the fur | planned as a source at a node within the fur, in place of `longwave_depth_fraction` |

## Where the chain is not a chain

Three features of the heat balance beyond the skin are not resistances in series. Each was measured in the
current code before deciding what to do about it.

**The skin temperature is written as the mean of two estimates.** One estimate comes from the core side, the
core temperature less the heat flow times the resistance of flesh and fat. The other comes from the
environment side, the fur surface temperature plus the heat flow times the resistance of the fur. The
`residual_skin_temperature` of [`solve_part_heat_balance`](@ref) is the difference between the skin temperature
and their mean. In a network the skin node is fixed by one condition, that the heat arriving equals the heat
leaving. The two are the same thing: at a solution the two estimates agree to about 10⁻¹³ K, for bodies from
0.1 to 100 kg, with and without contact with the ground and with wet skin. The mean is an unusual way of writing
continuity of heat flow at the skin, and a node balance gives the same numbers.

**Contact with the ground is a branch.** Heat leaves the skin by two paths side by side, through the free fur to
the air and through the compressed fur to the substrate, as in the diagram above. It is not a small term. For a
1 kg cylinder under a sky at −10 °C on a substrate at 27 °C:

| Fraction of the surface on the ground | Share of the heat exchanged, other than by evaporation, that passes through the ground |
|:--|:--|
| 0.1 | 8 % |
| 0.3 | 26 % |
| 0.5 | 45 % |

with the compressed fur a little under 1 K warmer than the free fur surface. So the structure must be a graph,
with a node for the substrate reached from the skin, and not a list.

**Radiation can act within the fur.** `longwave_depth_fraction` places the exchange of longwave radiation at a
depth in the coat, and the closed form then solves conduction and radiation at that depth together. Every
example and test uses a fraction of 1, the outer surface, where this coupling does nothing. For other values
the closed form divides by ``\ln(R_{fa} / R_{rad})``, which goes to zero as the depth approaches the surface, and
for thin fur the solution runs away. So radiation at depth is off by default and not usable when on. With the fur
divided into shells, radiation at a depth becomes heat entering at a node inside the fur, with nothing to divide
by. The same holds for sunlight, which in a deep or sparse coat is absorbed below the surface.

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
For an optimiser, each node temperature is a variable with one residual, which generalises the two variables,
skin and fur surface, that each part has now.

The same idea, nodes joined by conductances with a source of heat at each node, is applied sideways between the
parts of a body in [Bodies of many parts](multipart.md). The radial direction and the lateral one are two uses of
one structure.
