# The ellipsoid model

Before the full endotherm model there was a simpler one: an animal as a furred ellipsoid in an environment
where air, ground and sky are at one temperature and there is no sun (Porter and Kearney 2009). The heat budget
then has a closed-form solution, with no iteration, and shows how size, shape and fur set the range of
temperatures in which an endotherm can live. It is `ellipsoid_endo` in NicheMapR and
[`ellipsoid_endotherm`](@ref) here.

```@setup ellipsoid
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## The model

Heat made evenly through an ellipsoid of flesh passes through four resistances: the flesh, the fur, and then
convection and radiation side by side at the outer surface:

```math
Q_{gen} = \frac{T_c - T_a}{R_{flesh} + R_{fur} + \dfrac{R_{conv} \, R_{rad}}{R_{conv} + R_{rad}}}
```

This is the network of [Gradients, resistances and flows](../manual/gradients.md#The-budget-as-a-network) at
its simplest: resistors in series and in parallel, one source of flow and one source of effort.

The fur has a fixed conductivity. The shape is a prolate ellipsoid whose long axis is `posture` times its short
axes: a posture of 1 is a ball, and a larger one an animal stretched out.

```@example ellipsoid
using BiophysicalGeometry, Unitful # hide
shape_gallery(("posture $posture" => Body(Ellipsoid(0.5u"kg", 1000.0u"kg/m^3", posture, posture), # hide
    CompositeInsulation(FibrousLayer(5.0u"mm", 30.0u"μm", 3000u"cm^-2"), FatLayer(0.0, 901.0u"kg/m^3"))) for posture in (1.1, 2.0, 4.5))...) # hide
```

A 500 g animal with 5 mm of fur in three postures, with part of the fur cut away.

```@example ellipsoid
using HeatExchange, Unitful

ellipsoid(air_temperature; mass = 0.5u"kg", posture = 4.5, insulation_depth = 5.0u"mm") = ellipsoid_endotherm(
    u"K"(air_temperature * u"°C"), 0.1u"m/s", 0.5, 101325.0u"Pa";     # air temperature, wind speed, humidity, pressure
    mass, posture, insulation_depth,
    density = 1000.0u"kg/m^3",
    insulation_conductivity = 0.04u"W/m/K",
    emissivity = 0.95,
    core_temperature = u"K"(37.0u"°C"),
    q10 = 3.0,
    oxygen_fraction = 0.2094,
    oxygen_extraction_efficiency = 0.2,
    stress_factor = 0.6,
)

out = ellipsoid(5.0)
keys(out)
```

```@example ellipsoid
out.required_metabolic_heat_production, out.basal_metabolic_rate_fraction, out.skin_temperature
```

A 500 g animal with 5 mm of fur, stretched out at 5 °C, needs three and a half times its basal metabolic rate,
which is estimated from its mass by the equation of Kleiber (1947) unless `minimum_metabolic_rate` is given.

## The thermoneutral zone

Because the heat required is the temperature difference over a fixed resistance, the air temperature at which
it equals the basal rate can be written down directly. That is the *lower critical temperature*:

```@example ellipsoid
out.lower_critical_air_temperature, out.upper_critical_air_temperature
```

Above it the animal makes more heat than it loses and must lose the rest by evaporating water. The *upper
critical temperature* is defined here as the air temperature at which the heat to be lost reaches a fraction,
`stress_factor`, of the basal rate. Between the two is the thermoneutral zone.

```@example ellipsoid
air_temperatures = 5.0:1.0:45.0
results = ellipsoid.(air_temperatures)

fig = Figure(size = (760, 320))
ax1 = Axis(fig[1, 1]; xlabel = "Air temperature (°C)", ylabel = "Metabolic rate (W)")
lines!(ax1, air_temperatures, [ustrip(u"W", r.final_metabolic_heat_production) for r in results]; linewidth = 2)
vlines!(ax1, [ustrip(u"°C", out.lower_critical_air_temperature), ustrip(u"°C", out.upper_critical_air_temperature)];
        color = :black, linestyle = :dash)
ax2 = Axis(fig[1, 2]; xlabel = "Air temperature (°C)", ylabel = "Water loss (g/h)")
lines!(ax2, air_temperatures, [ustrip(u"g/hr", r.total_water_loss_rate) for r in results]; linewidth = 2, label = "total")
lines!(ax2, air_temperatures, [ustrip(u"g/hr", r.respiratory_water_loss_rate) for r in results]; linewidth = 2, label = "respiratory")
axislegend(ax2; position = :lt)
fig
```

`final_metabolic_heat_production` is the heat required, or the basal rate where that is higher. The water loss
is a rough estimate: the heat that must be lost above the lower critical temperature, as water evaporated.

## Size, shape and fur

The critical temperatures depend on the three things the model has: size, posture and fur.

```@example ellipsoid
masses = 10 .^ range(-2, 3; length = 40)   # kg
fig = Figure(size = (760, 320))
ax1 = Axis(fig[1, 1]; xlabel = "Mass (kg)", ylabel = "Lower critical temperature (°C)", xscale = log10)
for depth in (2.0, 10.0, 30.0)
    lines!(ax1, masses, [ustrip(u"°C", ellipsoid(20.0; mass = m * u"kg", posture = 2.0, insulation_depth = depth * u"mm").lower_critical_air_temperature)
                         for m in masses]; linewidth = 2, label = "$depth mm of fur")
end
axislegend(ax1; position = :lb)
ax2 = Axis(fig[1, 2]; xlabel = "Mass (kg)", ylabel = "Lower critical temperature (°C)", xscale = log10)
for posture in (1.1, 2.0, 4.5)
    lines!(ax2, masses, [ustrip(u"°C", ellipsoid(20.0; mass = m * u"kg", posture, insulation_depth = 10.0u"mm").lower_critical_air_temperature)
                         for m in masses]; linewidth = 2, label = "posture $posture")
end
axislegend(ax2; position = :lb)
fig
```

A large animal has less surface for its mass and holds its heat. With the same depth of fur, its lower
critical temperature is far below that of a small one, and a small animal can only match it with fur it could
not carry. Curling up lowers the critical temperature by a few degrees at any size. These are the patterns from
which Porter and Kearney (2009) drew general conclusions about the thermal niches of endotherms.

## Compared with NicheMapR

The output of `ellipsoid_endo` for the same animal over a range of air temperatures is kept with the tests:

```@example ellipsoid
using DelimitedFiles

path(file) = joinpath(pkgdir(HeatExchange), "test", "data", file)
table, header = readdlm(path("ellipsoid_output.csv"), ','; header = true)
column(name) = Float64.(table[:, findfirst(==(name), vec(header))])
reference = ellipsoid.(column("AirTemp"))

maximum(abs.([ustrip(u"W", r.required_metabolic_heat_production) for r in reference] .- column("Qgen"))),
maximum(abs.([ustrip(u"°C", r.lower_critical_air_temperature) for r in reference] .- column("LCT")))
```

The largest differences in the heat required, in W, and in the lower critical temperature, in °C, are within
rounding.

## The ellipsoid model and the full model

The ellipsoid model leaves out what makes the full model hard: fur whose conductivity depends on its
temperature, radiation within the fur, sunlight, a sky and ground at different temperatures, a coat that
differs between back and belly, contact with the ground, and evaporation from the skin as part of the balance.
With those absent or constant, the full model of [the endotherm tutorial](endotherm.md) reduces to nearly the
same thing. Its fixed conductivities are descriptive where the full model is mechanistic, see
[Units, dimensions and functional traits](../manual/units_traits.md#Processes-and-sub-processes).

The ellipsoid model is useful as a first estimate, for broad comparisons across many species, and for teaching.
[`solve_metabolic_rate`](@ref) is for anything in a real environment.
