# Introduction

An organism is an open thermodynamic system. Heat, water and the other inputs and outputs of metabolism cross
its surface, and by the first law these flows must balance or be stored (Porter and Gates 1969, Kearney and
Porter 2020). HeatExchange.jl computes the flows of heat, and the flows of water and gas that go with them, and
finds the state at which they balance.

```@setup intro
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## What is computed

Heat comes in from metabolism and from solar and longwave radiation. It leaves by emitted longwave radiation and
by evaporation, and is lost or gained by convection and conduction:

```@example intro
heat_budget_diagram() # hide
```

Each term depends on the temperature of the surface, which depends on the temperature of the core and on how
well heat is conducted between them. The package answers one of two questions:

- Given the metabolic heat production, **what body temperature** balances the budget? The question for an
  ectotherm, a leaf, or an endotherm that is not regulating. See [`solve_temperature`](@ref).
- Given the core temperature, **what metabolic rate** balances the budget? The question for a mammal or bird.
  The answer is its energy cost, with the water it must evaporate. See [`solve_metabolic_rate`](@ref).

Both are steady states: no heat is stored. See [Solving a heat balance](heat_balance.md) for how they are found,
[Temperature or metabolic rate](solvers.md) for how to call them, and
[Gradients, resistances and flows](gradients.md) for what the terms have in common.

## Organisms

An [`Organism`](@ref) has two parts:

```@example intro
using HeatExchange, BiophysicalGeometry, Unitful

shape = Cylinder(40.0u"g", 1000.0u"kg/m^3", 6.0)
body = Body(shape, Naked())
organism = Organism(body, example_ectotherm_heat_exchange_traits(; shape_pars = shape))
nothing # hide
```

```@example intro
fur = FibrousLayer(20.0u"mm", 30.0u"μm", 3000u"cm^-2") # hide
layers = CompositeInsulation(fur, FatLayer(0.2, 901.0u"kg/m^3")) # hide
shape_gallery("Cylinder" => Body(Cylinder(2.0u"kg", 1000.0u"kg/m^3", 2.0), layers), # hide
    "Sphere" => Body(BiophysicalGeometry.Sphere(2.0u"kg", 1000.0u"kg/m^3"), layers), # hide
    "Ellipsoid" => Body(Ellipsoid(2.0u"kg", 1000.0u"kg/m^3", 2.0, 2.0), layers)) # hide
```

Bodies of three shapes, each of 2 kg, with part of the fur and fat cut away to show the flesh.

The **body** is a shape with its layers of fat and fur, from
[BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl). It supplies every area,
length and volume:

| From the body | Used for |
|:--|:--|
| `total_area`, `skin_area`, `evaporation_area` | radiation, convection and evaporation |
| `silhouette` | sunlight absorbed from the direct beam |
| `flesh_radius`, `skin_radius`, `insulation_radius`, `flesh_volume` | conduction from the core through flesh, fat and fur |
| the family of the shape (cylinder, sphere, ellipsoid, plate) | the equations for convection and for conduction |

The **traits** are a [`HeatExchangeTraits`](@ref), one parameter struct for each process:

| Field | Type | Process |
|:--|:--|:--|
| `insulation_pars` | [`InsulationParameters`](@ref) | [fur and feathers](insulation.md) |
| `conduction_pars_external` | [`ExternalConductionParameters`](@ref) | [contact with the ground](convection_conduction.md) |
| `conduction_pars_internal` | [`InternalConductionParameters`](@ref) | [flesh and fat](radial_layers.md) |
| `radiation_pars` | [`RadiationParameters`](@ref) | [solar and longwave radiation](radiation.md) |
| `convection_pars` | [`ConvectionParameters`](@ref) | [convection](convection_conduction.md) |
| `evaporation_pars` | [`AnimalEvaporationParameters`](@ref) or [`LeafEvaporationParameters`](@ref) | [evaporation](evaporation_respiration.md) |
| `hydraulic_pars` | [`HydraulicParameters`](@ref) | [water potential](evaporation_respiration.md) |
| `respiration_pars` | [`RespirationParameters`](@ref) | [respiration](evaporation_respiration.md) |
| `metabolism_pars` | [`MetabolismParameters`](@ref) | [metabolism](metabolism.md) |
| `options` | [`SolveMetabolicRateOptions`](@ref) | solver tolerances |

See [Parameters](parameters.md).

Nothing in an organism says whether it is an ectotherm or an endotherm, an animal or a plant. Those are
combinations of traits and of the question asked. A leaf is a thin plate with
[`LeafEvaporationParameters`](@ref), solved for temperature, see [A leaf](../tutorials/leaf.md). An endotherm
is any body solved for metabolic rate.

## Environments

An environment is a NamedTuple of the conditions around the organism, [`EnvironmentalVars`](@ref), and the
properties of the site, [`EnvironmentalPars`](@ref):

```@example intro
environment = (;
    environment_pars = example_environment_pars(),
    environment_vars = EnvironmentalVars(;
        air_temperature = u"K"(20.0u"°C"),
        sky_temperature = u"K"(-5.0u"°C"),
        ground_temperature = u"K"(30.0u"°C"),
        substrate_temperature = u"K"(30.0u"°C"),
        relative_humidity = 0.2,
        wind_speed = 1.0u"m/s",
        atmospheric_pressure = 101325.0u"Pa",
        zenith_angle = 20.0u"°",
        substrate_conductivity = 0.5u"W/m/K",
        global_radiation = 1000.0u"W/m^2",
        diffuse_fraction = 0.1,
        shade = 0.0,
    ),
)
u"°C"(solve_temperature(organism, environment).core_temperature)
```

These are the conditions where the organism is: air temperature and wind speed at its height, the ground under
it, the sky above it. They come from a microclimate model, see
[Environments and the ecosystem](ecosystem.md).

## Units

All inputs and outputs are [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities, in any units of
the right dimension. Fractions (relative humidity, shade, skin wetness) are numbers from 0 to 1. Temperatures
are absolute: convert with `u"K"(20.0u"°C")`, and show a result with `u"°C"(temperature)`. See
[Units, dimensions and functional traits](units_traits.md) for what the units say about the equations.

## Origins

The package is a re-design, in Julia, of the heat budget code of
[NicheMapR](https://github.com/mrke/NicheMapR). The ectotherm model (Kearney and Porter 2020) goes back to
Porter, Mitchell, Beckman and DeWitt (1973). The endotherm model (Kearney et al. 2021) brought together
distributed heat generation in the flesh (Porter and Kearney 2009), heat transfer through porous fur (Conley and
Porter 1986) and the simultaneous solution for skin and fur temperatures under solar radiation (Mathewson and
Porter 2013). The package reproduces their results, and the tutorials show the comparisons.

Kearney et al. (2021) found that the endotherm model usually needed structural changes to suit a species, so
presented it as modules to be combined, and sketched its extension to several body parts. This package takes
that further:

- the two models are one set of heat-flow functions and two solvers;
- layers of tissue and insulation are a list, see [Layers as a radial graph](radial_layers.md);
- a body can have any number of parts, see [Bodies of many parts](multipart.md);
- thermoregulation is the domain of
  [BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), by rules or by
  optimisation.

See [For NicheMapR users](nichemapr.md).
