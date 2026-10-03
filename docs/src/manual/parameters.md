# Parameters

The traits of an organism and the properties of its environment are held in structs, one for each process. This
page lists them with their defaults, and the functions that build ready-made sets.

```@setup parameters
using Main.FigureHelpers
using CairoMakie
using HeatExchange, BiophysicalGeometry, Unitful
```

## Traits

A [`HeatExchangeTraits`](@ref) holds one of each struct below. Build it by position or with an `example_`
function:

```@example parameters
using HeatExchange, BiophysicalGeometry, Unitful

shape = Ellipsoid(65.0u"kg", 1000.0u"kg/m^3", 1.1, 1.1)
traits = example_heat_exchange_traits(; shape_pars = shape,
    radiation_pars = example_radiation_pars(; body_absorptivity_dorsal = 0.9),
    evaporation_pars = example_evaporation_pars(; skin_wetness = 0.01))
organism = Organism(Body(shape, Naked()), traits)
radiation_pars(organism).body_absorptivity_dorsal
```

Each group is read back by a function named for its field: [`insulation_pars`](@ref),
[`conduction_pars_external`](@ref), [`conduction_pars_internal`](@ref), [`radiation_pars`](@ref),
[`convection_pars`](@ref), [`evaporation_pars`](@ref), [`hydraulic_pars`](@ref), [`respiration_pars`](@ref),
[`metabolism_pars`](@ref) and [`options`](@ref).

The defaults below are those of the structs themselves.

::: tabs

== Radiation

[`RadiationParameters`](@ref), see [Radiation](radiation.md).

```@example parameters
parameter_table(RadiationParameters()) # hide
```

`silhouette_area`, `total_area` and `conduction_area` are not used by the solvers, which take these from the body.

== Conduction

[`ExternalConductionParameters`](@ref) and [`InternalConductionParameters`](@ref), see
[Convection and conduction](convection_conduction.md) and [Layers as a radial graph](radial_layers.md).

```@example parameters
parameter_table(ExternalConductionParameters()) # hide
```

```@example parameters
parameter_table(InternalConductionParameters()) # hide
```

== Convection

[`ConvectionParameters`](@ref), see [Convection and conduction](convection_conduction.md).

```@example parameters
parameter_table(ConvectionParameters()) # hide
```

== Evaporation

[`AnimalEvaporationParameters`](@ref), or for a leaf [`LeafEvaporationParameters`](@ref), and
[`HydraulicParameters`](@ref), see [Evaporation and respiration](evaporation_respiration.md).

```@example parameters
parameter_table(AnimalEvaporationParameters()) # hide
```

```@example parameters
parameter_table(LeafEvaporationParameters()) # hide
```

```@example parameters
parameter_table(HydraulicParameters()) # hide
```

== Respiration

[`RespirationParameters`](@ref), see [Evaporation and respiration](evaporation_respiration.md).

```@example parameters
parameter_table(RespirationParameters()) # hide
```

== Metabolism

[`MetabolismParameters`](@ref), see [Metabolism](metabolism.md).

```@example parameters
parameter_table(MetabolismParameters()) # hide
```

== Insulation

[`InsulationParameters`](@ref) and its [`FibreProperties`](@ref), see [Insulation](insulation.md).

```@example parameters
parameter_table("Dorsal" => InsulationParameters().dorsal, "Ventral" => InsulationParameters().ventral) # hide
```

== Options

[`SolveMetabolicRateOptions`](@ref): whether respiration is included, and the tolerances of the surface solve and
of the respiration balance.

```@example parameters
parameter_table(SolveMetabolicRateOptions()) # hide
```

:::

## Environment

[`EnvironmentalPars`](@ref) holds the properties of the site:

```@example parameters
parameter_table(EnvironmentalPars()) # hide
```

[`EnvironmentalVars`](@ref) holds the conditions. It has no defaults, except that `reference_air_temperature`,
`bush_temperature` and `vegetation_temperature` are the air temperature unless given:

| Field | Meaning |
|:--|:--|
| `air_temperature` | air temperature at the height of the organism |
| `reference_air_temperature` | air temperature at the reference height, taken as that of vegetation overhead |
| `sky_temperature` | radiant temperature of the sky |
| `ground_temperature` | temperature of the ground the organism faces |
| `substrate_temperature` | temperature of the surface the organism touches |
| `bush_temperature`, `vegetation_temperature` | temperature of vegetation beside and above the organism |
| `relative_humidity` | at the height of the organism, 0 to 1 |
| `wind_speed` | at the height of the organism |
| `atmospheric_pressure` | atmospheric pressure |
| `zenith_angle` | angle of the sun from overhead |
| `global_radiation` | solar radiation on a horizontal surface |
| `diffuse_fraction` | fraction of the solar radiation that is diffuse, 0 to 1 |
| `shade` | fraction of the sky above the organism shaded by vegetation, 0 to 1 |
| `substrate_conductivity` | thermal conductivity of the substrate |

[`EnvironmentalVarsVec`](@ref) has the same fields as vectors, for a series of hours. See
[Environments and the ecosystem](ecosystem.md) for where the values come from.

## Example sets

Two sets of functions return parameters ready for use, each taking keywords for the values to change.

- `example_…`: the defaults of the NicheMapR endotherm model `endoR`. A 65 kg ellipsoid with 2 mm of fur and a
  core at 37 °C, in still, dry air with sky and ground at air temperature and no sun (Kearney et al. 2021).
- `example_ectotherm_…`: the defaults of the NicheMapR ectotherm model. A 40 g lizard with a tenth of its
  surface on the ground and skin that is 0.1 % wet (Kearney and Porter 2020).

| Function | Returns |
|:--|:--|
| [`example_heat_exchange_traits`](@ref), [`example_ectotherm_heat_exchange_traits`](@ref) | a whole [`HeatExchangeTraits`](@ref) |
| [`example_shape_pars`](@ref), [`example_ellipsoid_shape_pars`](@ref) | a shape |
| [`example_insulation_pars`](@ref) | [`InsulationParameters`](@ref) |
| [`example_conduction_pars_external`](@ref), [`example_ectotherm_conduction_pars_external`](@ref) | [`ExternalConductionParameters`](@ref) |
| [`example_conduction_pars_internal`](@ref), [`example_ectotherm_conduction_pars_internal`](@ref) | [`InternalConductionParameters`](@ref) |
| [`example_radiation_pars`](@ref), [`example_ectotherm_radiation_pars`](@ref) | [`RadiationParameters`](@ref) |
| [`example_convection_pars`](@ref) | [`ConvectionParameters`](@ref) |
| [`example_evaporation_pars`](@ref), [`example_ectotherm_evaporation_pars`](@ref), [`example_leaf_evaporation_pars`](@ref) | evaporation parameters |
| [`example_hydraulic_pars`](@ref), [`example_ectotherm_hydraulic_pars`](@ref) | [`HydraulicParameters`](@ref) |
| [`example_respiration_pars`](@ref), [`example_ectotherm_respiration_pars`](@ref) | [`RespirationParameters`](@ref) |
| [`example_metabolism_pars`](@ref), [`example_ectotherm_metabolism_pars`](@ref) | [`MetabolismParameters`](@ref) |
| [`example_metabolic_rate_options`](@ref) | [`SolveMetabolicRateOptions`](@ref) |
| [`example_environment_vars`](@ref), [`example_environment_pars`](@ref) | [`EnvironmentalVars`](@ref), [`EnvironmentalPars`](@ref) |

```@example parameters
parameter_table("Endotherm" => example_radiation_pars(), "Ectotherm" => example_ectotherm_radiation_pars()) # hide
```

The ranges over which an animal can change these, and the thresholds at which it does, are the parameters of
BiophysicalBehaviour.jl, see
[Parameters](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/parameters) in its
documentation.

## Bounds, and fitting

Each parameter is a functional trait in the sense of Kearney et al. (2021), see
[Units, dimensions and functional traits](units_traits.md).

The default of each field is a `Param` of [ModelParameters.jl](https://github.com/rafaqz/ModelParameters.jl),
carrying a value with its units and, for many fields, bounds:

```@example parameters
RadiationParameters().sky_view_factor
```

A model built from structs of `Param`s can be shown as a table, and its parameters read and set as a vector,
which is what a sensitivity analysis or an optimiser needs. The bounds are those within which a trait makes
physical sense. The solvers remove the wrappers before computing, so plain numbers and quantities work as well,
and the `example_` functions return them.
