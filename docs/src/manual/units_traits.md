# Units, dimensions and functional traits

Every quantity in this package carries units. This page is about what the units reveal: which equations are
physical and which descriptive, where the two meet in the code, and what that says about functional traits.

```@setup units
using Main.FigureHelpers
using CairoMakie
```

## Units as a check on the physics

A physical equation is dimensionally consistent whatever units are used. Convective heat loss is a heat
transfer coefficient times an area times a temperature difference, and with
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) the product is a power without anyone saying so:

```@example units
using HeatExchange, BiophysicalGeometry, Unitful

heat_transfer_coefficient, area, temperature_difference = 4.0u"W/m^2/K", 780.0u"cm^2", 10.0u"K"
uconvert(u"W", heat_transfer_coefficient * area * temperature_difference)
```

The area was in cm² and the coefficient per m². An equation with a term missing, or a radius where a diameter
was meant, has the wrong dimension and fails at once. The whole heat budget is checked this way each time it
runs.

The parameters of such an equation can be interpreted. A thermal conductivity, in W m⁻¹ K⁻¹, is the heat that
crosses a metre of material for each degree across it. It is a property of the material, can be measured on a
sample, and means the same in a mouse and an elephant.

## Descriptive equations have no such dimensions

The basal metabolic rate of a mammal is computed as (Kleiber 1947)

```math
Q_{basal} = 3.39 \, M^{0.75}
```

in W, with ``M`` in kg. Mass to the power 0.75 has units that mean nothing:

```@example units
(65.0u"kg")^0.75
```

and the coefficient would need units of W kg⁻⁰·⁷⁵. No physical quantity has those dimensions. The coefficient
is the intercept of a line fitted to many species, and its value depends on the units they were plotted in. The
only way to evaluate it is to strip the units, compute, and put units back:

```julia
function metabolic_rate(::Kleiber, mass, core_temperature=nothing)
    mass_kg = ustrip(u"kg", mass)
    m_rate_kcal_day = 70.0 * mass_kg^0.75
    return (m_rate_kcal_day * 4.185 / (24 * 3.6))u"W"
end
```

**Wherever `ustrip` is needed to evaluate an equation, the model has left parameters with interpretable
dimensions, tied to a physical mechanism, and is using a description of data in its place.**

The distinction is not between good and bad equations. Kleiber's is among the best supported relations in
biology. It is between an equation that says *how* something comes about and one that says *what* is observed.

- A **descriptive** equation rests on the data it was fitted to. It holds within their range and for the kinds
  of organism in them. To use it outside that range is to extrapolate.
- A **physical** equation does not rest on a set of data in that way. Its parameters are properties of the
  organism and environment, so applying it to a new size, shape or environment is not extrapolation. What must
  hold is that the mechanism still operates.

This is the reason to build a model from physical and chemical principles where they are available, and to mark
the places where they are not (Kearney and Porter 2020).

## Where units are stripped in this package

There are two reasons for `ustrip` in the source. Only one marks a descriptive equation.

**For arithmetic.** Root finders and linear solvers work on plain numbers. The residual is stripped to W before
[`zbrent`](@ref), and the conductance matrix of a body of several parts to W/K. Units go back on the answer and
nothing about the equations has changed. See [Differentiability and the NLP interface](autodiff.md).

**For description.** Here an empirical relation stands in for a mechanism:

| Where | Relation | In place of |
|:--|:--|:--|
| [`metabolic_rate`](@ref) with [`Kleiber`](@ref), [`McKechnieWolf`](@ref), [`AndrewsPough2`](@ref) | power laws of mass, and ``10^{0.038 \, T_b}`` with ``T_b`` in °C | the chemical transformations that produce the heat |
| [`metabolic_rate`](@ref) with [`PlantDarkRespiration`](@ref) | proportional to mass, with an Arrhenius term | the respiration of a leaf |
| the basal rate in [`ellipsoid_endotherm`](@ref) | Kleiber's equation with a ``Q_{10}`` | as above |
| the conductivity of the body in [`ellipsoid_endotherm`](@ref) | a linear function of the radius in m | heat moved within the body by tissue and blood |
| the cap on respiratory water loss in [`respiration`](@ref) | a maximum per kg of body mass (Welch 1980) | the behaviour of an animal that would move away |
| `DesertIguana` and `LeopardFrog` of BiophysicalGeometry.jl | surface area as a power of mass | the shape of the animal |

The last shows in the units of the result, which no length has:

```@example units
Body(DesertIguana(40.0u"g", 1000.0u"kg/m^3"), Naked()).geometry.length
```

Everything else, radiation, conduction through flesh, fat and fur, evaporation and the molar balance of the
lungs, is computed without removing a unit.

### Dimensionless numbers

Convection sits between the two. The heat transfer coefficient comes from relations such as
``Nu = 0.35 \, Re^{0.6}``, fitted to wind-tunnel measurements and so descriptive. But the Nusselt and Reynolds
numbers are dimensionless, so coefficient and exponent are pure numbers and no units are stripped:

```@example units
air_density, wind_speed, dimension, viscosity = 1.2u"kg/m^3", 2.0u"m/s", 5.0u"cm", 1.8e-5u"Pa*s"
reynolds_number = uconvert(NoUnits, air_density * wind_speed * dimension / viscosity)
```

That is what dimensional analysis buys. Physics decides which combinations of variables matter, and the data
are fitted only as a relation between them. One relation then serves any fluid, speed and size at which the
flow is similar: a lizard in air, a fish in water. The fitted part is confined to the shape. See
[Convection and conduction](convection_conduction.md).

## `allometric` as the marker

[BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl) collects the descriptive
relations of the ecosystem behind one function, `allometric`:

```@example units
import BiologicalScaling: allometric, power_law, reference, BasalMetabolicRate, EutherianMammal

allometric(BasalMetabolicRate(), EutherianMammal(), 65.0u"kg")
```

The relation behind it says outright what a descriptive equation is: a coefficient, an exponent, the unit the
input must be in, the unit of the output, and the source of the data:

```@example units
relation = power_law(BasalMetabolicRate(), EutherianMammal())
(; relation.coefficient, relation.exponent, relation.input_unit, relation.output_unit)
```

```@example units
reference(relation)
```

The removal of units happens once, in one place, with the fitted units recorded beside it.

The name is a flag. A call to `allometric` marks a point where a process is described and not computed from its
parts. The metabolic rate equations of this package are of that kind and are moving to BiologicalScaling.jl,
see [Metabolism](metabolism.md). So are the empirical surface areas now in BiophysicalGeometry.jl: a `Cylinder`
gives an area from geometry and a `DesertIguana` from a regression, and the two should not look alike.

## Processes and sub-processes

A process is *mechanistic* if it is captured as a set of sub-processes, and *phenomenological* if not. The
terms are relative, and apply at each level of a model.

The heat budget is mechanistic at its top level: body temperature is not fitted to air temperature but follows
from radiation, convection, conduction, evaporation and metabolism, see
[Solving a heat balance](heat_balance.md). Each of those can be asked the same question.

- **Conduction through fur** is a given conductivity in the [ellipsoid model](../tutorials/ellipsoid.md), and
  is computed from the depth of the coat and the properties of its fibres in the full model, see
  [Insulation](insulation.md).
- **Conduction through flesh** uses a conductivity that stands for blood flow to the skin. The blood is not
  modelled, see [Bodies of many parts](multipart.md).
- **Evaporation from the skin** uses the fraction of skin that acts as a free water surface. The resistance of
  the skin, the sweat glands and the wetting of fur behind that number are not explicit.
- **Metabolic heat**, when it is not what is solved for, is a power law of mass and a function of temperature.
  It is the least mechanistic term in the budget.

For metabolism the sub-processes are known. Dynamic Energy Budget theory (Kooijman 2010) gives the metabolic
rate as the energy dissipated in somatic and maturity maintenance, maturation, and the overheads of
assimilation, growth and reproduction. Each has parameters with interpretable dimensions. Oxygen consumption,
carbon dioxide, water and nitrogenous waste follow too, which a power law does not give, and the allometric
scaling of metabolic rate becomes a result, not an input. This is how NicheMapR computes metabolism with its DEB
model switched on (Kearney and Porter 2020), and it is where the `allometric` call here points. The replacement
is to come through [AnimalMapper.jl](https://github.com/BiophysicalEcology/AnimalMapper.jl), by way of
[DEBtool_J.jl](https://github.com/add-my-pet/DEBtool_J.jl), see
[Flows of mass](gradients.md#Flows-of-mass).

When metabolic rate is the unknown, as for an endotherm in the cold, the question does not arise:
[`solve_metabolic_rate`](@ref) finds the heat required from the heat budget alone, and the allometric equation
gives only the minimum.

## Functional traits

Kearney et al. (2021) argued that the link between a trait and the performance of an organism should be made
through a model of the organism as a thermodynamic system. A trait is *functional* if it has a defined role in
such a model, and *descriptive* if not.

### State variables are not traits

Such a model is a dynamical system, with terms of four kinds (Kearney et al. 2021):

| Term | What it is | Here |
|:--|:--|:--|
| **state variable** | defines the state of the organism, and changes through time | the temperatures of the core, the skin and the fur surface |
| **process** | a rate at which the state changes, or at which energy or matter crosses the boundary | each heat flow; the metabolic rate; water loss; oxygen consumption |
| **environmental variable** | an aspect of the environment that acts on the state | [`EnvironmentalVars`](@ref) |
| **parameter** | sets how the state responds to the environment | the parameter structs, see [Parameters](parameters.md) |

The package *computes* the state and the processes. Body temperature, from [`solve_temperature`](@ref), is a
state. Metabolic rate, from [`solve_metabolic_rate`](@ref), and water loss are processes. None is a trait. A
body temperature of 31 °C belongs to a lizard *in a place at a time*, and a metabolic rate of 112 W to an animal
*in air at 0 °C*. A table of body temperatures or field metabolic rates is a table of states and processes
under conditions mostly not recorded, and cannot be given to a model as a property of a species.

What can be a trait is a *particular value* of a state or process, one at which something happens: the body
temperature at which an animal dies, loses coordination, or stops foraging and seeks shade; the lowest
metabolic rate it can have; the highest rate at which it can sweat. These are threshold traits. They turn a
computed state into a consequence: an hour of activity, a response, a death.

A threshold is also one end of a gradient: the difference between the present state and a threshold is what an
organism responds to, see [Gradients of information](gradients.md#Gradients-of-information).

The heat budget has few thresholds. It returns a body temperature of 60 °C without comment. The thresholds
belong to [BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), see
[States, thresholds and traits](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/states_traits)
in its documentation.

### Four kinds of functional trait

Kearney et al. (2021) distinguished four roles a trait can have in a model:

| Class | Role in the model | In HeatExchange.jl |
|:--|:--|:--|
| **parameter** | a parameter of the model | the fields of the parameter structs: reflectance, emissivity, conductivities, wet fraction of the skin, oxygen extraction, target core temperature, see [Parameters](parameters.md) |
| **threshold** | a value of a state or process at which the organism fails or changes its behaviour | the minimum metabolic rate and the target core temperature. The rest are in BiophysicalBehaviour.jl |
| **model** | decides the structure of the model | `Naked()` or a `FibrousLayer`; the family of the shape; [`AnimalEvaporationParameters`](@ref) or [`LeafEvaporationParameters`](@ref); solved for temperature or for metabolic rate; the number of parts and their coupling |
| **estimation** | an observation from which a parameter is estimated | fibre diameter, length, density and depth, for fur conductivity; mass, density and proportions, for areas; oxygen consumption, for metabolic heat |

Two consequences show in how the package is used.

**Some familiar measures are not traits.** The lower and upper critical temperatures of an endotherm are air
temperatures, outside the organism, and hold only for the wind, humidity, radiation and posture of the
experiment (Kearney et al. 2021). Here they are *outputs*:
[An endotherm: metabolic rate](../tutorials/endotherm.md) finds the lower critical temperature from the traits,
and shows it move by many degrees with posture. What to measure instead is what the heat budget needs: shape
and posture, size, areas, the fur, reflectance, the wet fraction of the skin, the target core temperature and
the minimum metabolic rate. A metabolic chamber still gives the last directly, and tests the model, if wind,
humidity and posture are reported.

**A descriptive trait can become functional with the right metadata.** Colour says something about solar heat
load, but the heat budget needs reflectance over the solar spectrum. Coat depth is of no use without fibre
diameter and density. The model says which measurements belong together. The units and bounds each parameter
carries, see [Parameters](parameters.md), are part of that metadata, and the dimensions of a parameter are the
first check that a value from the literature is the quantity the model means.

The structure of the package follows the classification:

- parameter and estimation traits are the data of an [`Organism`](@ref);
- model traits are its types, on which methods are chosen;
- environmental variables are a separate object, [`EnvironmentalVars`](@ref);
- state variables and processes are what the solvers return.
