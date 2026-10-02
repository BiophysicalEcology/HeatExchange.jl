# An ectotherm: body temperature

The body temperature of an ectotherm is set by the heat that it exchanges with its surroundings. This tutorial
computes the steady-state body temperature and water loss of a lizard, takes its heat budget apart, compares it
with the ectotherm model of NicheMapR (Kearney and Porter 2020), and then asks how the answer changes with the
sun, the wind, the size of the animal and the way it holds itself.

```@setup ectotherm
using Main.FigureHelpers
using CairoMakie
```

## The lizard

The animal is the default of the NicheMapR function `ectoR_devel`: a 40 g lizard with the shape of a desert iguana
(Porter et al. 1973), whose surface area follows a relation measured for that species:

```@example ectotherm
using HeatExchange, BiophysicalGeometry, Unitful

shape = DesertIguana(40.0u"g", 1000.0u"kg/m^3")
body = Body(shape, Naked())
uconvert(u"cm^2", total_area(body)), uconvert(u"cm^2", silhouette(shape, NormalToSun())),
uconvert(u"cm^2", silhouette(shape, ParallelToSun()))
```

`DesertIguana` is an empirical relation of area to mass, not a geometric shape, and it and `LeopardFrog` are to
move from BiophysicalGeometry.jl to BiologicalScaling.jl. The last two are the area that intercepts the direct beam of the sun when the lizard is side-on to it and when it
points at it. Its traits are those of [`example_ectotherm_heat_exchange_traits`](@ref), with the lizard side-on to
the sun and its eyes closed:

```@example ectotherm
lizard_traits(; shape = shape, orientation = NormalToSun(), skin_wetness = 0.001) = example_ectotherm_heat_exchange_traits(;
    shape_pars = shape,
    radiation_pars = example_ectotherm_radiation_pars(; solar_orientation = orientation),
    evaporation_pars = example_ectotherm_evaporation_pars(; skin_wetness, eye_fraction = 0.0),
)
lizard = Organism(body, lizard_traits())
nothing # hide
```

| Trait | Value | In |
|:--|:--|:--|
| solar absorptivity | 0.85 | [`RadiationParameters`](@ref) |
| emissivity | 0.95 | |
| view of the sky and of the ground | 0.4 each | |
| fraction of the surface on the ground | 0.1 | [`ExternalConductionParameters`](@ref) |
| thermal conductivity of flesh | 0.5 W m⁻¹ K⁻¹ | [`InternalConductionParameters`](@ref) |
| fraction of the skin that is wet | 0.001 | [`AnimalEvaporationParameters`](@ref) |
| water potential of the body | −707 J/kg | [`HydraulicParameters`](@ref) |
| oxygen extraction efficiency, respiratory quotient | 0.2, 0.8 | [`RespirationParameters`](@ref) |
| metabolic rate | Andrews and Pough (1985), standard | [`MetabolismParameters`](@ref) |

## The environment

The lizard is on warm ground under a clear sky in strong sun, with a light breeze and very dry air, in light
shade:

```@example ectotherm
import HeatExchange: GasFractions

environment_vars = EnvironmentalVars(;
    air_temperature = u"K"(20.0u"°C"),
    sky_temperature = u"K"(-5.0u"°C"),
    ground_temperature = u"K"(30.0u"°C"),
    substrate_temperature = u"K"(30.0u"°C"),
    relative_humidity = 0.05,
    wind_speed = 1.0u"m/s",
    atmospheric_pressure = 101325.0u"Pa",
    zenith_angle = 20.0u"°",
    substrate_conductivity = 0.5u"W/m/K",
    global_radiation = 1000.0u"W/m^2",
    diffuse_fraction = 0.1,
    shade = 0.1,
)
environment_pars = example_environment_pars(; gas_fractions = GasFractions(0.2095, 0.0003, 0.7902))
environment = (; environment_pars, environment_vars)
nothing # hide
```

## The simplest operation

```@example ectotherm
out = solve_temperature(lizard, environment)
u"°C"(out.core_temperature), u"°C"(out.surface_temperature), u"°C"(out.lung_temperature)
```

The lizard settles 11 °C above the air. Its core is a twentieth of a degree warmer than its skin, as its
metabolic heat is small. The terms of the heat budget, `out.energy_balance`, are the `enbal` table of NicheMapR:

```@example ectotherm
flow_table(out.energy_balance) # hide
```

```@example ectotherm
b = out.energy_balance # hide
budget_bars(["Solar" => (b.solar_flow, FLOW_COLOURS.solar), "Longwave in" => (b.longwave_flow_in, FLOW_COLOURS.longwave), # hide
             "Metabolism" => (b.metabolic_heat_flow, FLOW_COLOURS.metabolism)], # hide
            ["Longwave out" => (b.longwave_flow_out, FLOW_COLOURS.longwave), "Convection" => (b.convection_heat_flow, FLOW_COLOURS.convection), # hide
             "Conduction" => (b.conduction_flow, FLOW_COLOURS.conduction), "Evaporation" => (b.evaporation_heat_flow, FLOW_COLOURS.evaporation), # hide
             "Respiration" => (b.respiration_heat_flow, FLOW_COLOURS.respiration)]) # hide
```

Radiation dominates. The lizard absorbs sunlight and longwave radiation, and loses heat as longwave radiation and
by convection to the cooler air. Conduction to the ground is small, as the ground is close to the temperature of the lizard. Metabolism is a few hundredths of a watt. The heat budget of a small ectotherm in the sun is, to a
close approximation, a balance between radiation and convection.

The water lost and the oxygen used are in `out.mass_balance`, the `masbal` table of NicheMapR:

```@example ectotherm
m = out.mass_balance
uconvert(u"ml/hr", m.oxygen_flow), uconvert(u"mg/hr", m.respiration_mass), uconvert(u"mg/hr", m.cutaneous_mass)
```

## Compared with NicheMapR

The same animal and environment were given to `ectoR_devel` of NicheMapR, and its output is kept with the tests
of this package:

```@example ectotherm
using DelimitedFiles

table = readdlm(joinpath(pkgdir(HeatExchange), "test", "data", "ectoR_output.csv"), ','; skipstart = 1)
saved = Dict(String(table[i, 1]) => Float64(table[i, 2]) for i in axes(table, 1))
enbal = ("QSOL", "QIRIN", "QMET", "QRESP", "QEVAP", "QIROUT", "QCONV", "QCOND")   # the order of the enbal table
masbal = ("O2_ml", "H2OResp_g", "H2OCut_g")
nichemapr = merge(saved, Dict(name => saved["enbal$i"] for (i, name) in enumerate(enbal)),
                  Dict(name => saved["masbal$i"] for (i, name) in enumerate(masbal)))

markdown_table(["Quantity", "HeatExchange.jl", "NicheMapR"], [ # hide
    ("Core temperature", celsius(out.core_temperature), nichemapr["TC"] * u"°C"), # hide
    ("Skin temperature", celsius(out.surface_temperature), nichemapr["TSKIN"] * u"°C"), # hide
    ("Solar radiation absorbed", b.solar_flow, nichemapr["QSOL"] * u"W"), # hide
    ("Longwave radiation absorbed", b.longwave_flow_in, nichemapr["QIRIN"] * u"W"), # hide
    ("Longwave radiation emitted", b.longwave_flow_out, nichemapr["QIROUT"] * u"W"), # hide
    ("Convection", b.convection_heat_flow, nichemapr["QCONV"] * u"W"), # hide
    ("Conduction", uconvert(u"W", b.conduction_flow), nichemapr["QCOND"] * u"W"), # hide
    ("Evaporation", b.evaporation_heat_flow, nichemapr["QEVAP"] * u"W"), # hide
    ("Metabolism", b.metabolic_heat_flow, nichemapr["QMET"] * u"W"), # hide
    ("Respiration", b.respiration_heat_flow, nichemapr["QRESP"] * u"W"), # hide
    ("Oxygen consumption", uconvert(u"ml/hr", m.oxygen_flow), nichemapr["O2_ml"] * u"ml/hr"), # hide
    ("Cutaneous water loss", uconvert(u"mg/hr", m.cutaneous_mass), nichemapr["H2OCut_g"] * 1000 * u"mg/hr"), # hide
]) # hide
```

## Sun and shade

Shade is the main thing that a lizard can change. Here the same lizard is solved from full sun to full shade.
In this simple case only the sunlight is shaded, and the temperatures of the air and ground are left as they were:

```@example ectotherm
with(; kw...) = (; environment_pars, environment_vars = EnvironmentalVars(; (; (name => getfield(environment_vars, name)
                    for name in fieldnames(typeof(environment_vars)))..., kw...)...))
shades = 0.0:0.05:1.0
orientations = ("Side-on to the sun" => NormalToSun(), "Intermediate" => Intermediate(), "Pointing at the sun" => ParallelToSun())

fig, ax = figure_axis("Shade (%)", "Body temperature (°C)")
for (label, orientation) in orientations
    animal = Organism(body, lizard_traits(; orientation))
    temperatures = [ustrip(u"°C", solve_temperature(animal, with(shade = s)).core_temperature) for s in shades]
    lines!(ax, 100 .* shades, temperatures; linewidth = 2, label)
end
hlines!(ax, [20.0]; color = :black, linestyle = :dash)
axislegend(ax; position = :rt)
fig
```

The dashed line is the air temperature. In full shade the lizard is close to it, a little below because it
faces a cold sky. In the sun it can change its temperature by several degrees by turning. These are the
choices that the thermoregulation routines of NicheMapR make in a fixed order: change posture, seek shade, climb,
and go underground (Kearney and Porter 2020). They are made here by
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), which calls
[`solve_temperature`](@ref) as this loop does.

## Wind and size

Convection ties the body to the air temperature, and how strongly depends on the wind and on the size of the
animal:

```@example ectotherm
wind_speeds = 10 .^ range(-1, 1; length = 30)
masses = (4.0u"g", 40.0u"g", 400.0u"g", 4000.0u"g")

fig, ax = figure_axis("Wind speed (m/s)", "Body temperature (°C)"; xscale = log10)
for mass in masses
    sized = DesertIguana(mass, 1000.0u"kg/m^3")
    animal = Organism(Body(sized, Naked()), lizard_traits(; shape = sized))
    temperatures = [ustrip(u"°C", solve_temperature(animal, with(wind_speed = v * u"m/s")).core_temperature) for v in wind_speeds]
    lines!(ax, wind_speeds, temperatures; linewidth = 2, label = string(mass))
end
axislegend(ax; position = :rt)
fig
```

In still air a large lizard in the sun runs far hotter than a small one, and wind cools them all towards the air.
A steady state is a fair description of a 40 g lizard, which comes to a new temperature in minutes. A 4 kg lizard
takes much longer, and a transient heat budget is then needed, see
[Solving a heat balance](../manual/heat_balance.md#Steady-state-and-storage).

## Wet skin

A frog, with a skin that is wet all over, loses heat by evaporation that a lizard does not:

```@example ectotherm
wetness = 10 .^ range(-3, 0; length = 30)
results = [solve_temperature(Organism(body, lizard_traits(; skin_wetness = w)), environment) for w in wetness]

fig = Figure(size = (760, 320))
ax1 = Axis(fig[1, 1]; xlabel = "Fraction of the skin that is wet", ylabel = "Body temperature (°C)", xscale = log10)
ax2 = Axis(fig[1, 2]; xlabel = "Fraction of the skin that is wet", ylabel = "Water loss (g/h)", xscale = log10)
lines!(ax1, wetness, [ustrip(u"°C", r.core_temperature) for r in results]; linewidth = 2)
lines!(ax2, wetness, [ustrip(u"g/hr", r.mass_balance.cutaneous_mass) for r in results]; linewidth = 2)
fig
```

A wet-skinned animal of this size in this sun and dry air is many degrees cooler than a dry one, and pays for it
with a water loss of a large part of its mass each hour. The heat budget and the water budget cannot be
separated (Tracy 1976).

## Beyond the heat budget

The NicheMapR ectotherm model does much more than this: it runs through a year of microclimates, chooses where
the animal is and what it is doing each hour, and can grow the animal with a Dynamic Energy Budget model. Here
those are separate. The microclimate is from
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl), see
[Environments and the ecosystem](../manual/ecosystem.md), and the behaviour from BiophysicalBehaviour.jl. This
package is the calculation at the centre of each hour.
