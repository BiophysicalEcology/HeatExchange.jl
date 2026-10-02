# Introduction

An organism can be seen as an open thermodynamic system. Its boundary is its outer surface, and across that
boundary flow heat, water and the other inputs and outputs of metabolism. By the first law of thermodynamics these
flows must balance, or the difference is stored in the organism (Porter and Gates 1969, Kearney and Porter 2020).
HeatExchange.jl computes the flows of heat, and the evaporation of water that is part of them, and finds the state
of the organism at which they balance.

```@setup intro
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## What is computed

The heat budget of an organism has inputs from metabolism and from solar and longwave radiation. Heat is lost by
the longwave radiation that the surface emits and by evaporation from the skin and the lungs, and is lost or gained
by convection to the air and conduction to the ground:

```@example intro
heat_budget_diagram() # hide
```

Each term depends on the temperature of the surface of the organism, and that depends on the temperature of its
core and on how well heat is conducted between the two. The package answers one of two questions:

- For a given rate of metabolic heat production, **what body temperature** balances the budget? This is the
  question for an ectotherm, a leaf, or an endotherm that is not regulating. See [`solve_temperature`](@ref).
- For a given core temperature, **what metabolic rate** balances the budget? This is the question for a mammal or
  a bird, and the answer is its energy cost, with the water it must evaporate. See [`solve_metabolic_rate`](@ref).

Both are steady-state solutions: no heat is stored. [Solving a heat balance](heat_balance.md) explains how they are
found, and [Temperature or metabolic rate](solvers.md) how to call them.

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

The **traits** are a [`HeatExchangeTraits`](@ref), a set of parameter structs, one for each process:

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
[`LeafEvaporationParameters`](@ref), solved for temperature. An endotherm is any body solved for metabolic rate.

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

These are the conditions at the place where the organism is: the air temperature and wind speed at its height,
the temperature of the ground under it, and the radiant temperature of the sky above it. They are the output of a
microclimate model, see [Environments and the ecosystem](ecosystem.md).

## Units

All inputs and outputs are [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities, in any units of
the right dimension. Fractions, such as relative humidity, shade and skin wetness, are numbers between 0 and 1.
Temperatures are absolute: a temperature in °C is converted with `u"K"(20.0u"°C")`, and a result is shown in °C
with `u"°C"(temperature)`. What the units say about the equations of the model, and about its parameters as
functional traits, is the subject of [Units, dimensions and functional traits](units_traits.md).

## Origins

The package is a re-design, in Julia, of the heat budget code of
[NicheMapR](https://github.com/mrke/NicheMapR): the ectotherm model (Kearney and Porter 2020), which goes back to
the program of Porter, Mitchell, Beckman and DeWitt (1973), and the endotherm model (Kearney et al. 2021), which
brought together distributed heat generation in the flesh (Porter and Kearney 2009), heat transfer through porous
fur (Conley and Porter 1986) and the simultaneous solution for skin and fur temperatures under solar radiation
(Mathewson and Porter 2013). It reproduces their results, and the tutorials show the comparisons.

Kearney et al. (2021) found that the endotherm model, unlike the ectotherm model, usually needed structural
changes to suit a species, and so presented it as modules to be combined, and sketched how it could be extended to
several body parts. This package takes that further. The two models are one set of heat-flow functions and two
solvers. Layers of tissue and insulation are a list. A body can have any number of parts. And thermoregulation is
is now the domain of [BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), 
which can apply it by rules or by optimisation. See [For NicheMapR users](nichemapr.md).
