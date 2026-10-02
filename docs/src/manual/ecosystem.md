# Environments and the ecosystem

HeatExchange.jl solves the heat budget of one organism in one set of conditions. It does not know where the
conditions came from, what the organism does about them, or what the result means for its growth and survival.
Those are the work of the other packages of the [BiophysicalEcology](https://github.com/BiophysicalEcology)
ecosystem, which together make a mechanistic niche model (Kearney and Porter 2009, 2020):

```@raw html
<div style="text-align:center; font-family: var(--vp-font-family-mono); font-size: 0.85em; line-height: 2.0;">
climate and terrain<br>
↓ <i>Microclimate.jl, MicroclimateMapper.jl, SolarRadiation.jl</i><br>
<b>conditions where the organism is</b><br>
↓ <i>HeatExchange.jl</i>, with bodies from <i>BiophysicalGeometry.jl</i><br>
<b>body temperature, metabolic rate, water loss</b><br>
↕ <i>BiophysicalBehaviour.jl</i> chooses the place, the posture and the physiological state<br>
<b>activity, energy and water budgets</b><br>
↓ growth, development and reproduction
</div>
```

## Where the environment comes from

The variables in [`EnvironmentalVars`](@ref) are those at the organism, not those reported by a weather station.
Air temperature and wind speed change steeply in the first metre above the ground. The ground under a lizard may
be 30 °C hotter than the air at 2 m, and the clear sky above it 20 °C colder. A microclimate model computes these
from weather, terrain and soil, and each variable has a source:

| [`EnvironmentalVars`](@ref) | What it is | Source |
|:--|:--|:--|
| `air_temperature`, `wind_speed`, `relative_humidity` | at the height of the organism | the profiles of [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl), at the chosen height |
| `reference_air_temperature` | at the reference height, about 2 m, taken as the temperature of the vegetation that casts shade | Microclimate.jl |
| `ground_temperature`, `substrate_temperature` | the surface that the organism sees below it, and the one it touches | the soil surface temperature of Microclimate.jl, or the soil temperature at the depth of a burrow |
| `sky_temperature` | the radiant temperature of the sky | Microclimate.jl |
| `global_radiation`, `diffuse_fraction`, `zenith_angle` | sunlight on a horizontal surface, the part of it that is scattered, and the angle of the sun from overhead | [SolarRadiation.jl](https://github.com/BiophysicalEcology/SolarRadiation.jl), by itself or through Microclimate.jl |
| `shade` | the fraction of the sky above the organism that is vegetation | chosen, or found by BiophysicalBehaviour.jl |
| `substrate_conductivity` | of the soil that the organism lies on | the soil thermal properties of Microclimate.jl |
| `atmospheric_pressure` | | from elevation, with [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl) |
| `bush_temperature`, `vegetation_temperature` | of vegetation beside and above the organism | air temperature at its height and at the reference height |

A microclimate model run gives these for every hour of a day or a year, for full sun and for deep shade, above the
ground and at depths in the soil. One heat budget is then solved for each hour and each place. In outline, for a
lizard on the surface in the sun, with the names of the output of Microclimate.jl:

```julia
using Microclimate, HeatExchange, Unitful

micro = solve(example_microclimate_problem())

body_temperature = map(eachindex(micro.sky_temperature)) do hour
    environment_vars = EnvironmentalVars(;
        air_temperature = micro.profile.air_temperature[hour, 1],      # at the lowest height
        wind_speed = micro.profile.wind_speed[hour, 1],
        relative_humidity = micro.profile.relative_humidity[hour, 1],
        sky_temperature = micro.sky_temperature[hour],
        ground_temperature = micro.soil_temperature[hour, 1],          # the soil surface
        substrate_temperature = micro.soil_temperature[hour, 1],
        global_radiation = micro.global_radiation[hour],
        diffuse_fraction = micro.diffuse_fraction[hour],
        zenith_angle, atmospheric_pressure, substrate_conductivity, shade = 0.0,
    )
    solve_temperature(lizard, (; environment_pars, environment_vars)).core_temperature
end
```

This is a sketch and is not run here. The documentation of Microclimate.jl gives the full description of its
output.

### Over space

[MicroclimateMapper.jl](https://github.com/BiophysicalEcology/MicroclimateMapper.jl) runs Microclimate.jl over
rasters and sets of points. It fetches the terrain, weather, soil and land-cover data for each place from gridded
datasets, through [RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl) and
[Rasters.jl](https://github.com/rafaqz/Rasters.jl), and returns the same microclimate output for every cell or
point. A map of body temperature, of hours of possible activity or of the energy and water costs of an endotherm
is then a heat budget solved for each cell and each hour.

## What the organism does about it

An animal is not fixed in one place and one state. A lizard moves between sun and shade, changes its posture, and
goes underground. A mammal changes the depth of its fur, its posture and the blood flow to its skin, and when
those are not enough it pants or sweats (Kearney et al. 2021).
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl) models these, and each is
a change to the organism or to the environment that is given to this package, followed by another solve. Each
changes either a resistance to a flow of heat or the gradient that drives it, see
[Solving a heat balance](gradients.md):

| Response | What changes |
|:--|:--|
| seek shade, climb, go underground | [`EnvironmentalVars`](@ref) |
| orient to the sun | `solar_orientation` in [`RadiationParameters`](@ref) |
| curl up or stretch out | the axis ratio of the shape, and so the body |
| raise or flatten the fur | the depth of the fibrous layer |
| dilate or constrict blood vessels in the skin | `flesh_conductivity` |
| pant | `pant` in [`RespirationParameters`](@ref) |
| sweat, lick | `skin_wetness` in [`AnimalEvaporationParameters`](@ref) |
| let the core temperature rise | `core_temperature` in [`MetabolismParameters`](@ref) |

BiophysicalBehaviour.jl also assembles the heat budgets of bodies with several parts, and decides what an animal
does while its body temperature is changing. The transient heat budget itself, for an animal too large to be at
steady state, belongs to this package, see [Solving a heat balance](heat_balance.md#Steady-state-and-storage). See
[Differentiability and the NLP interface](autodiff.md) for how this package is written to be driven in that way.

## The packages

| Package | Link with HeatExchange.jl |
|:--|:--|
| [BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl) | Provides the body: its shape, layers, areas, radii, silhouette, and for bodies of several parts the joins and the views between parts |
| [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl) | Provides the properties of dry and humid air and of water used in convection, evaporation and respiration |
| [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) | Provides the conditions at the organism, above and below the ground, in sun and shade |
| [MicroclimateMapper.jl](https://github.com/BiophysicalEcology/MicroclimateMapper.jl) | Provides those conditions over rasters and sets of points, from gridded climate, terrain and soil data |
| [SolarRadiation.jl](https://github.com/BiophysicalEcology/SolarRadiation.jl) | Provides the direct and diffuse sunlight and the position of the sun |
| [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl) | Provides allometric equations for metabolic rate, surface area and body proportions. The metabolic rate equations now in this package are moving there, see [Metabolism](metabolism.md) |
| [ThermalPhysiology.jl](https://github.com/BiophysicalEcology/ThermalPhysiology.jl) | Uses the body temperatures found here in thermal performance curves and models of thermal death |
| [BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl) | Calls this package to find the response of an animal to its environment by behaviour and physiology, at steady state and through time |

These are to be brought together in [NicheMapper.jl](https://github.com/BiophysicalEcology/NicheMapper.jl) (in
development).
