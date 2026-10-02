# Get started

HeatExchange.jl solves the heat budget of an organism in a given environment, with every quantity a
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantity.

```julia
using Pkg
Pkg.add(url = "https://github.com/BiophysicalEcology/HeatExchange.jl")
```

```@setup get_started
using Main.FigureHelpers
using CairoMakie
```

## An organism

An [`Organism`](@ref) is a body and a set of traits. The body is a shape from
[BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl), sized from a mass and a
density. The traits are the surface properties and physiology that the heat budget needs. Here is a 40 g lizard
with bare skin and the default traits of an ectotherm:

```@example get_started
using HeatExchange, BiophysicalGeometry, Unitful

shape = DesertIguana(40.0u"g", 1000.0u"kg/m^3")
body = Body(shape, Naked())
lizard = Organism(body, example_ectotherm_heat_exchange_traits(; shape_pars = shape))
nothing # hide
```

## An environment

An environment has two parts: `environment_vars`, the conditions that change from hour to hour, and
`environment_pars`, the properties of the site. Here the lizard is in the sun on a cool morning with a light
wind:

```@example get_started
environment_vars = example_environment_vars(;
    air_temperature = u"K"(20.0u"°C"),
    wind_speed = 1.0u"m/s",
    global_radiation = 400.0u"W/m^2",
    zenith_angle = 30.0u"°",
)
environment_pars = example_environment_pars()
environment = (; environment_pars, environment_vars)
nothing # hide
```

`example_environment_vars` sets the sky, ground and substrate to the air temperature unless told otherwise. See
[Parameters](manual/parameters.md) for every field, and
[Environments and the ecosystem](manual/ecosystem.md) for where real values come from.

## Solve for body temperature

A lizard makes little metabolic heat and does not adjust it to defend a core temperature, so its body
temperature is whatever balances the heat it gains and loses. [`solve_temperature`](@ref) finds it:

```@example get_started
out = solve_temperature(lizard, environment)
u"°C"(out.core_temperature) # default output is in Kelvin
```

The lizard is above air temperature. The terms of its heat budget show why:

```@example get_started
flow_table(out.energy_balance) # hide
```

```@example get_started
b = out.energy_balance # hide
budget_bars(["Solar" => (b.solar_flow, FLOW_COLOURS.solar), "Longwave in" => (b.longwave_flow_in, FLOW_COLOURS.longwave), # hide
             "Metabolism" => (b.metabolic_heat_flow, FLOW_COLOURS.metabolism)], # hide
            ["Longwave out" => (b.longwave_flow_out, FLOW_COLOURS.longwave), "Convection" => (b.convection_heat_flow, FLOW_COLOURS.convection), # hide
             "Conduction" => (b.conduction_flow, FLOW_COLOURS.conduction), "Evaporation" => (b.evaporation_heat_flow, FLOW_COLOURS.evaporation), # hide
             "Respiration" => (b.respiration_heat_flow, FLOW_COLOURS.respiration)]) # hide
```

Sunlight and longwave radiation from the sky and ground come in, and leave as longwave radiation, by convection
to the air, and by conduction to the ground, which is at air temperature here and touches a tenth of the lizard.
Metabolism and evaporation are a small part of this animal's budget. The water it loses is in
`out.mass_balance`. See [Solving a heat balance](manual/heat_balance.md).

## Solve for metabolic rate

A mammal or bird holds its core temperature steady, and changes the heat it produces. Here is a 65 kg animal
with 2 mm of fur, the default of [`example_heat_exchange_traits`](@ref), in still air at 0 °C:

```@example get_started
shape = Ellipsoid(65.0u"kg", 1000.0u"kg/m^3", 1.1, 1.1)
fibres = FibreProperties(; diameter = 30.0u"μm", length = 23.9u"mm", density = 3000.0u"cm^-2", depth = 2.0u"mm",
                           reflectance = 0.2, conductivity = 0.209u"W/m/K")
fur = FibrousLayer(fibres.depth, fibres.diameter, fibres.density)
fat = FatLayer(0.0, 901.0u"kg/m^3")
body = Body(shape, CompositeInsulation(fur, fat))

traits = example_heat_exchange_traits(;
    shape_pars = shape,
    insulation_pars = InsulationParameters(; dorsal = fibres, ventral = fibres, depth_compressed = fibres.depth,
                                             longwave_depth_fraction = 1.0),
)
mammal = Organism(body, traits)

cold = (; environment_pars = example_environment_pars(),
          environment_vars = example_environment_vars(; air_temperature = u"K"(0.0u"°C")))
nothing # hide
```

The fur appears twice. The body needs its depth, for the outer area and radius. The traits need the properties
of its fibres too, for how well it conducts heat, see [Insulation](manual/insulation.md).

[`solve_metabolic_rate`](@ref) finds the metabolic rate that holds the core at the temperature in
`metabolism_pars`, 37 °C here. It takes first guesses of the skin and fur surface temperatures:

```@example get_started
result = solve_metabolic_rate(mammal, cold, u"K"(34.0u"°C"), u"K"(0.0u"°C"))
result.energy_flows.metabolic_heat_flow
```

The basal metabolic rate of an animal of this mass is 77.6 W (Kleiber 1947), so it must produce about 1.4 times
basal to stay warm. The skin and the outer surface of the fur are found at the same time:

```@example get_started
u"°C"(result.thermoregulation.skin_temperature), u"°C"(result.thermoregulation.insulation_temperature)
```

The result has four groups of output, `thermoregulation`, `morphology`, `energy_flows` and `mass_flows`, see
[Temperature or metabolic rate](manual/solvers.md).

## Any combination

The two solvers are not tied to kinds of animal. The furred body can be solved for its temperature, that of an
animal that has stopped thermoregulating, with its metabolic rate at basal:

```@example get_started
torpid = solve_temperature(mammal, cold)
u"°C"(torpid.thermoregulation.core_temperature)
```

## Where next

- [Solving a heat balance](manual/heat_balance.md) explains the budget and how it is solved, and
  [Gradients, resistances and flows](manual/gradients.md) what its terms have in common.
- The tutorials work through [an ectotherm](tutorials/ectotherm.md), [an endotherm](tutorials/endotherm.md) and
  [a leaf](tutorials/leaf.md).
- [Bodies of many parts](manual/multipart.md) and the tutorials on [two halves](tutorials/two_parts.md) and
  [a human](tutorials/human.md) cover bodies with more than one part.
- [For NicheMapR users](manual/nichemapr.md) maps the models of NicheMapR onto this package.
- For what an animal does about its heat budget, see the documentation of
  [BiophysicalBehaviour.jl](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/).
