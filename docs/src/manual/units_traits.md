# Units, dimensions and functional traits

Every quantity in this package carries units. That is a convenience, and it is also a statement about the kind of
model this is. This page is about what the units reveal: which equations are physical and which are descriptive,
where the two meet in the code, and how that distinction defines what a functional trait is.

```@setup units
using Main.FigureHelpers
using CairoMakie
```

## Units as a check on the physics

A physical equation is dimensionally consistent: the two sides have the same dimensions, whatever units are used.
Convective heat loss is a heat transfer coefficient times an area times a temperature difference, and with
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) the product is a power without anyone saying so:

```@example units
using HeatExchange, BiophysicalGeometry, Unitful

heat_transfer_coefficient, area, temperature_difference = 4.0u"W/m^2/K", 780.0u"cm^2", 10.0u"K"
uconvert(u"W", heat_transfer_coefficient * area * temperature_difference)
```

The area was given in cm² and the coefficient per m². An equation with a term missing, or with a radius where a
diameter was meant, gives a result in the wrong dimension and fails at once. The whole heat budget is checked in
this way each time it runs.

The parameters of such an equation have dimensions that can be interpreted. A thermal conductivity, in
W m⁻¹ K⁻¹, is the heat that crosses a metre of the material for each degree across it. It is a property of the
material, it can be measured on a sample by itself, and it means the same thing in a mouse and in an elephant.

## Descriptive equations have no such dimensions

Not every equation in the package is of this kind. The basal metabolic rate of a mammal is computed as
(Kleiber 1947)

```math
Q_{basal} = 3.39 \, M^{0.75}
```

in W, with the mass ``M`` in kg. The mass raised to the power 0.75 has units that mean nothing:

```@example units
(65.0u"kg")^0.75
```

and for the equation to give watts the coefficient 3.39 would need units of W kg⁻⁰·⁷⁵. There is no physical
quantity with those dimensions. The coefficient is not a property of anything. It is the intercept of a line
fitted to measurements of many species, and its value depends on the units in which they were plotted. In the
code, the only way to evaluate such an equation is to take the units off the mass, compute with the bare number,
and put the units of the answer back on:

```julia
function metabolic_rate(::Kleiber, mass, core_temperature=nothing)
    mass_kg = ustrip(u"kg", mass)
    m_rate_kcal_day = 70.0 * mass_kg^0.75
    return (m_rate_kcal_day * 4.185 / (24 * 3.6))u"W"
end
```

**Wherever `ustrip` is needed to evaluate an equation, the model has left parameters with interpretable
dimensions, tied to an explicit physical mechanism, and is using a description of data in its place.** The need
to strip the units is the sign.

The distinction is not between good equations and bad ones. Kleiber's equation is one of the best supported
relations in biology. It is between an equation that says *how* something comes about and one that says *what*
is observed. A descriptive equation rests on the data that it was fitted to. It holds within the range of those data and for
the kinds of organism in them, and to use it outside that range is to extrapolate. A physical equation does not
rest on a set of data in that way. Its parameters are properties of the organism and its environment, so applying
it to a new size, shape or environment is not extrapolation: nothing is being carried beyond the observations
that define it. What must hold is that the mechanism it describes still operates. This is the reason for building
a model from physical and chemical principles where they are available, and for marking the places where they
are not (Kearney and Porter 2020).

## Where units are stripped in this package

There are two reasons for `ustrip` in the source, and only one of them marks a descriptive equation.

**For arithmetic.** A root finder and a linear solver work on plain numbers. The residual of the heat budget is
stripped to W before it is given to [`zbrent`](@ref), and the conductance matrix of a body of several parts is
stripped to W/K before it is solved. The units are put back on the answer, and nothing about the meaning of the
equations has changed. See [Differentiability and the NLP interface](autodiff.md).

**For description.** These are the places where an empirical relation stands in for a mechanism:

| Where | Relation | In place of |
|:--|:--|:--|
| [`metabolic_rate`](@ref) with [`Kleiber`](@ref), [`McKechnieWolf`](@ref), [`AndrewsPough2`](@ref) | power laws of mass, and ``10^{0.038 \, T_b}`` with ``T_b`` in °C | the chemical transformations that produce the heat |
| [`metabolic_rate`](@ref) with [`PlantDarkRespiration`](@ref) | proportional to mass, with an Arrhenius term | the respiration of a leaf |
| the basal rate in [`ellipsoid_endotherm`](@ref) | Kleiber's equation with a ``Q_{10}`` | as above |
| the conductivity of the body in [`ellipsoid_endotherm`](@ref) | a linear function of the radius in m | the movement of heat within the body by tissue and blood |
| the cap on respiratory water loss in [`respiration`](@ref) | a maximum per kg of body mass (Welch 1980) | the behaviour of an animal that would move away |
| `DesertIguana` and `LeopardFrog` of BiophysicalGeometry.jl | surface area as a power of mass | the shape of the animal |

The last can be seen in the units of the result, which no length has:

```@example units
Body(DesertIguana(40.0u"g", 1000.0u"kg/m^3"), Naked()).geometry.length
```

Everything else in the heat budget, radiation, conduction through flesh, fat and fur, evaporation and the molar
balance of the lungs, is computed without removing a unit.

### Dimensionless numbers

Convection sits between the two. The heat transfer coefficient comes from relations such as
``Nu = 0.35 \, Re^{0.6}``, which were fitted to measurements in wind tunnels and are descriptive. But the Nusselt
and Reynolds numbers are dimensionless, so the coefficient and the exponent are pure numbers, and no units are
stripped to evaluate them:

```@example units
air_density, wind_speed, dimension, viscosity = 1.2u"kg/m^3", 2.0u"m/s", 5.0u"cm", 1.8e-5u"Pa*s"
reynolds_number = uconvert(NoUnits, air_density * wind_speed * dimension / viscosity)
```

This is what dimensional analysis buys. The physics decides which combinations of variables matter, and the
measurements are fitted only as a relation between those combinations. One relation then serves for any fluid,
speed and size at which the flow is similar, which is why a correlation from a wind tunnel can be used for a
lizard in air and a fish in water. The fitted part is confined to the shape. See
[Convection and conduction](convection_conduction.md).

## `allometric` as the marker

[BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl) collects the descriptive
relations of the ecosystem behind one function, `allometric`:

```@example units
import BiologicalScaling: allometric, power_law, reference, BasalMetabolicRate, EutherianMammal

allometric(BasalMetabolicRate(), EutherianMammal(), 65.0u"kg")
```

The relation behind it is an object that says outright what a descriptive equation is: a coefficient and an
exponent, the unit that the input must be in, the unit that the output is in, and the source of the data:

```@example units
relation = power_law(BasalMetabolicRate(), EutherianMammal())
(; relation.coefficient, relation.exponent, relation.input_unit, relation.output_unit)
```

```@example units
reference(relation)
```

The removal of units, which is unavoidable, happens once, in one place, with the units that the relation was
fitted in recorded beside it.

The name is meant to be read as a flag. A call to `allometric` in a model marks a point at which a process is
being described and not computed from its parts. The metabolic rate equations of this package are of that kind,
and are moving to BiologicalScaling.jl so that they are reached through `allometric`, see
[Metabolism](metabolism.md). The empirical surface areas now in BiophysicalGeometry.jl are moving for the same
reason: a `Cylinder` gives an area from geometry, and a `DesertIguana` gives one from a regression, and the two
should not look alike.

## Processes and sub-processes

One way to define *mechanistic* against *phenomenological* is by whether a process is captured as a set of
sub-processes. By that definition the terms are relative, and apply at each level of a model.

The heat budget is mechanistic at its top level. Body temperature is not a function fitted to air temperature. It
follows from solar and longwave radiation, convection, conduction, evaporation and metabolism, each computed
from the properties of the organism and its surroundings, and summed, see
[Solving a heat balance](heat_balance.md).

Each of those terms can be asked the same question.

- **Conduction through fur** can be a conductivity given as a number, as in the
  [ellipsoid model](../tutorials/ellipsoid.md), or it can be computed from the depth of the coat and the
  diameter, length, density and conductivity of its fibres, with the radiation between them, as in the full
  model, see [Insulation](insulation.md). The second is the first with its sub-processes made explicit.
- **Conduction through flesh** uses a conductivity that is varied to stand for the blood flow to the skin. The
  blood is not modelled. A model of perfusion between the parts of a body would take this down a level, see
  [Bodies of many parts](multipart.md).
- **Evaporation from the skin** uses the fraction of the skin that acts as a free water surface. Behind that
  number are the resistance of the skin to water, sweat glands, and the wetting of fur, none of which is
  explicit.
- **Metabolic heat**, when it is not what is being solved for, is a power law of mass and a function of
  temperature. This is the least mechanistic term in the budget.

For metabolism the sub-processes are known. Dynamic Energy Budget theory (Kooijman 2010) gives the metabolic
rate of an organism as the sum of the energy dissipated in its chemical transformations: somatic maintenance,
maturity maintenance, maturation, and the overheads of assimilation, growth and reproduction. Each has its own
parameters, with interpretable dimensions, and the metabolic rate at any size, temperature and food level follows
from them. So does the rate of oxygen consumption, of carbon dioxide and water production, and of nitrogenous
waste, which a power law does not give. The allometric scaling of metabolic rate with mass is then a result and
not an input. This is how the ectotherm model of NicheMapR computes metabolism when its DEB model is switched on
(Kearney and Porter 2020), and it is where the `allometric` call for metabolic rate in this package points: at a
process that can be replaced by its parts.

When the metabolic rate is the unknown, as for an endotherm in the cold, the question does not arise.
[`solve_metabolic_rate`](@ref) finds the heat that must be produced from the heat budget alone, and the
allometric equation is used only for the minimum below which the animal cannot go.

## Functional traits

The same reasoning says what a functional trait is. Kearney et al. (2021) argued that the link between a trait
and the performance of an organism, its survival, development, growth and reproduction, should be made through
a model of the organism as a thermodynamic system. A trait is *functional* if it has a defined role in such a
model, and *descriptive* if it does not.

### State variables are not traits

Such a model is a dynamical system, and its terms are of four kinds (Kearney et al. 2021):

| Term | What it is | Here |
|:--|:--|:--|
| **state variable** | a quantity that defines the state of the organism, and changes through time | the temperatures of the core, the skin and the fur surface |
| **process** | a rate at which the state is changed, or at which energy or matter crosses the boundary of the organism | each heat flow of the budget; the metabolic rate; the rates of water loss from the skin and the lungs; the rate of oxygen consumption |
| **environmental variable** | an aspect of the environment that acts on the state | air, sky and ground temperature, wind speed, humidity, solar radiation, see [`EnvironmentalVars`](@ref) |
| **parameter** | a term, usually constant, that sets how the state responds to the environment | the properties of the organism in its parameter structs, see [Parameters](parameters.md) |

What this package *computes* is the state and the processes. Body temperature is the state variable that
[`solve_temperature`](@ref) solves for, the one at which the heat stored is zero. The metabolic rate that
[`solve_metabolic_rate`](@ref) returns, and the water loss that comes with either solution, are processes. None
of these is a trait. A body temperature of 31 °C belongs to a lizard *in a place at a time*, as the metabolic rate
of 112 W belongs to an animal *in air at 0 °C*. In another environment the same organism, with the same traits,
has another body temperature and another metabolic rate. A table of body temperatures or of field metabolic rates
is a table of states and processes under conditions that are mostly not recorded, and it cannot be given to a
model as a property of the species.

What can be a trait is a *particular value* of a state variable or a process: one at which something happens.
The body temperature at which an animal dies, or loses coordination, or stops foraging and seeks shade, is a
property of the organism, and it does not change with where the animal is. So is the lowest metabolic rate
that it can have, and the highest rate at which it can sweat. These are threshold traits. They are where the
state that this package computes is turned into a consequence for the organism: an hour in which it can be
active or cannot, a response that it must make, a death.

A threshold is also one end of a gradient. The difference between the present state and a threshold or target
state is what an organism responds to, a gradient of information where the heat budget has gradients of physical
potential, see [Solving a heat balance](gradients.md).

The heat budget itself has few thresholds. It will return a body temperature of 60 °C, or a metabolic rate
below the minimum, without comment. The thresholds belong to the model of what the organism does about its
state, which is [BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl).
[States, thresholds and traits](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/states_traits) in the documentation of that package describes them: the body
temperatures that bound activity and trigger each behaviour of an ectotherm, and the limits of core temperature,
flesh conductivity, panting and skin wetness within which an endotherm responds.

### Four kinds of functional trait

Kearney et al. (2021) distinguished four ways in which a trait can have a role in a model, and each has a place
in this package:

| Class | Role in the model | In HeatExchange.jl |
|:--|:--|:--|
| **parameter** | a parameter of the model | the fields of the parameter structs: solar reflectance of the fur, emissivity, flesh and fat conductivity, the fraction of the skin that is wet, oxygen extraction efficiency, target core temperature, see [Parameters](parameters.md) |
| **threshold** | a value of a state variable or process at which the organism fails or changes its behaviour | the minimum metabolic rate, below which a solution means that the animal cannot lose its heat, and the target core temperature. The rest are in BiophysicalBehaviour.jl |
| **model** | decides the structure of the model | whether the body is `Naked()` or has a `FibrousLayer`; the family of its shape; [`AnimalEvaporationParameters`](@ref) or [`LeafEvaporationParameters`](@ref); whether it is solved for temperature or for metabolic rate, which is ectothermy and endothermy; the number of its parts and how they are coupled |
| **estimation** | an observation from which a parameter is estimated, given the conditions of the measurement | fibre diameter, length, density and depth, from which the conductivity of the fur is computed; mass, density and proportions, from which areas are computed; oxygen consumption, from which metabolic heat is computed |

Two consequences of this view are visible in how the package is used.

**Some familiar measures are not traits.** The lower and upper critical temperatures of an endotherm are widely
tabulated as measures of its thermal tolerance. In this scheme they are neither parameters nor thresholds of the
organism. They are air temperatures, outside the organism, and they hold only for the wind, humidity, radiation
and posture of the experiment in which they were measured (Kearney et al. 2021). Here they are *outputs*. The
tutorial [An endotherm: metabolic rate](../tutorials/endotherm.md) finds the lower critical temperature of an
animal from its traits, and shows it move by many degrees as the animal changes posture. What the theory says to
measure instead is what the heat budget needs: shape and posture, size, areas, the properties of the fur, solar
reflectance, the wet fraction of the skin, the target core temperature and the minimum metabolic rate. A
metabolic chamber experiment still gives one of these directly, the basal metabolic rate, and it tests the
model, if the wind speed, humidity and posture are reported with it.

**A descriptive trait can become functional with the right metadata.** The colour of an animal says something
about its solar heat load, but the heat budget needs the reflectance over the whole solar spectrum. The depth of
a coat is of no use without the diameter and density of its fibres. The model states which measurements belong
together and what must be recorded with them. The units and bounds that each parameter carries, see
[Parameters](parameters.md), are a small part of that metadata, and the dimensions of a parameter are the first
check that a value from the literature is the quantity that the model means.

The structure of the package follows the classification. The parameter and estimation traits are the data of an
[`Organism`](@ref). The model traits are its types, on which the methods are chosen. The environmental variables
are a separate object, [`EnvironmentalVars`](@ref), because they are not properties of the organism. And the
state variables and processes, body temperature, metabolic rate and water loss, are what the solvers return.
