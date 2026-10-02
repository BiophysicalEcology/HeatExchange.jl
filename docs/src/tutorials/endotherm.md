# An endotherm: metabolic rate

An endotherm holds its core temperature, and the question for its heat budget is what that costs: the metabolic
rate, and the water, needed in a given environment. This tutorial follows the example of the NicheMapR endotherm
model (Kearney et al. 2021), with its default animal, and compares the results with `endoR_devel`.

```@setup endotherm
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## The animal

The default settings of `endoR_devel` are for a 65 kg ellipsoid with a pelt 2 mm deep and no layer of fat,
holding a core temperature of 37 °C. Its long axis is 1.1 times its short axes, so it is nearly a sphere, a
curled-up posture that conserves heat. It roughly approximates a typical-sized, lightly clothed, resting human
(Kearney et al. 2021). The same defaults are in the `example_` functions:

```@example endotherm
using HeatExchange, BiophysicalGeometry, Unitful

function animal(; mass = 65.0u"kg", axis_ratio = 1.1, fur_depth = 2.0u"mm")
    shape = Ellipsoid(mass, 1000.0u"kg/m^3", axis_ratio, axis_ratio)
    insulation_pars = example_insulation_pars(; insulation_depth_dorsal = fur_depth, insulation_depth_ventral = fur_depth,
                                                insulation_depth_compressed = fur_depth)
    fibres = insulation_pars.dorsal
    body = Body(shape, CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density), FatLayer(0.0, 901.0u"kg/m^3")))
    traits = example_heat_exchange_traits(; shape_pars = shape, insulation_pars,
        metabolism_pars = example_metabolism_pars(; metabolic_heat_flow = metabolic_rate(Kleiber(), mass)))
    return Organism(body, traits)
end

mammal = animal()
basal = metabolism_pars(mammal).metabolic_heat_flow
```

The minimum metabolic rate is the basal rate of a mammal of this mass, ``3.39 \, M^{0.75}`` (Kleiber 1947).

The environment is the default as well: the ground and sky at air temperature, a wind of 0.1 m/s, a relative
humidity of 5 % and no sun. It is an indoor environment, or a metabolic chamber:

```@example endotherm
chamber(air_temperature) = (; environment_pars = example_environment_pars(),
                              environment_vars = example_environment_vars(; air_temperature = u"K"(air_temperature)))
nothing # hide
```

## The simplest operation

At an air temperature of 0 °C, with first guesses for the skin and fur temperatures:

```@example endotherm
solve(air_temperature) = solve_metabolic_rate(mammal, chamber(air_temperature), u"K"(34.0u"°C"), u"K"(air_temperature))
out = solve(0.0u"°C")
out.energy_flows.metabolic_heat_flow, out.energy_flows.metabolic_heat_flow / basal
```

The animal must produce about 1.4 times its basal rate to hold 37 °C. The output has four groups, see
[Temperature or metabolic rate](../manual/solvers.md).

`energy_flows` has the components of the energy balance, in W, the `enbal` table of NicheMapR:

```@example endotherm
e = out.energy_flows
(; e.solar_flow, e.longwave_flow_in, e.longwave_flow_out, e.convection_heat_flow, e.conduction_flow, e.evaporation_heat_flow,
   e.metabolic_heat_flow)
```

`thermoregulation` has the temperatures, the `treg` table. The lung temperature is the mean of the core and skin
temperatures:

```@example endotherm
t = out.thermoregulation
u"°C"(t.core_temperature), u"°C"(t.lung_temperature), u"°C"(t.skin_temperature), u"°C"(t.insulation_temperature)
```

```@example endotherm
t.insulation_conductivity_effective, t.dorsal.insulation_conductivity
```

The skin is at about 20 °C and the outer surface of the fur at 13 °C. The second conductivity includes the
radiation within the fur, see [Insulation](../manual/insulation.md).

`mass_flows` has the mass balance, the `masbal` table: the air and oxygen passing through the lungs, and the
water lost in the breath and from the skin:

```@example endotherm
m = out.mass_flows
uconvert(u"l/hr", m.air_flow), uconvert(u"l/hr", m.oxygen_flow_standard), m.respiration_mass_flow, m.m_sweat
```

`morphology` has the lengths, areas and volumes, the `morph` table:

```@example endotherm
g = out.morphology
g.area_skin, g.total_area, g.area_evaporation, g.characteristic_dimension
```

The skin area is 0.783 m². After taking away the area covered by the bases of the hairs, the area that can
evaporate water is 0.766 m², and with a skin wetness of 0.5 % only 0.0038 m² of it acts as a free water surface.

## Across air temperatures

```@example endotherm
air_temperatures = 0.0:1.0:40.0   # °C
sweep = [solve(T * u"°C") for T in air_temperatures]
rates = [ustrip(u"W", r.energy_flows.metabolic_heat_flow) for r in sweep]

fig = Figure(size = (760, 620))
ax1 = Axis(fig[1, 1]; ylabel = "Metabolic rate required (W)")
lines!(ax1, air_temperatures, rates; linewidth = 2)
hlines!(ax1, [ustrip(u"W", basal)]; color = :black, linestyle = :dash)
ax2 = Axis(fig[2, 1]; xlabel = "Air temperature (°C)", ylabel = "Temperature (°C)")
lines!(ax2, air_temperatures, [ustrip(u"°C", r.thermoregulation.skin_temperature) for r in sweep]; linewidth = 2, label = "skin")
lines!(ax2, air_temperatures, [ustrip(u"°C", r.thermoregulation.insulation_temperature) for r in sweep]; linewidth = 2, label = "fur surface")
lines!(ax2, air_temperatures, collect(air_temperatures); color = :black, linestyle = :dot, label = "air")
axislegend(ax2; position = :lt)
linkxaxes!(ax1, ax2)
hidexdecorations!(ax1; grid = false)
fig
```

The metabolic rate required falls in a nearly straight line as the air warms, as Scholander's classic picture of
an endotherm has it, and crosses the basal rate, the dashed line, at the *lower critical temperature*:

```@example endotherm
lower_critical = air_temperatures[findfirst(<(ustrip(u"W", basal)), rates)]
```

Below it the animal is spending energy to keep warm. Above it the result is below the basal rate, which no animal
can achieve: it makes at least its basal heat, and cannot lose it all in this state. The line above the lower
critical temperature is therefore not a prediction of what the animal does. It is a statement that the animal must
change.

## What the animal does next

In `endoR_devel` a sequence of responses follows, each tried until it is used up (Kearney et al. 2021): the
animal uncurls, which raises its surface area; it sends more blood to the skin, which raises the conductivity of
the flesh; it lets its core temperature rise; it pants; and it sweats. Each is a change to a trait followed by
another solve of the heat budget, and they are in
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), by rules as in NicheMapR
or by optimisation. The effect of each on the heat budget can be seen here by making the change by hand.

**Posture.** The default animal is curled into a ball. Stretched out, with its long axis 4 times its short axes,
it has more surface, and must produce more heat in the cold. Its lower critical temperature is higher:

```@example endotherm
shape_gallery(("axis ratio $ratio" => body(animal(; axis_ratio = ratio, fur_depth = 20.0u"mm")) for ratio in (1.1, 2.0, 4.0))...) # hide
```

The animals are drawn with 20 mm of fur so that it can be seen.

```@example endotherm
required(organism, T) = ustrip(u"W", solve_metabolic_rate(organism, chamber(T * u"°C"), u"K"(34.0u"°C"), u"K"(T * u"°C")).energy_flows.metabolic_heat_flow)

fig, ax = figure_axis("Air temperature (°C)", "Metabolic rate required (W)")
for ratio in (1.1, 2.0, 4.0)
    lines!(ax, air_temperatures, [required(animal(; axis_ratio = ratio), T) for T in air_temperatures]; linewidth = 2,
           label = "axis ratio $ratio")
end
hlines!(ax, [ustrip(u"W", basal)]; color = :black, linestyle = :dash)
axislegend(ax; position = :rt)
fig
```

**Fur.** A deeper coat lowers the cost of the cold and the lower critical temperature together:

```@example endotherm
fig, ax = figure_axis("Air temperature (°C)", "Metabolic rate required (W)")
for depth in (2.0, 10.0, 30.0)
    lines!(ax, air_temperatures, [required(animal(; fur_depth = depth * u"mm"), T) for T in air_temperatures]; linewidth = 2,
           label = "$depth mm of fur")
end
hlines!(ax, [ustrip(u"W", basal)]; color = :black, linestyle = :dash)
axislegend(ax; position = :rt)
fig
```

**Size.** The same in animals from 100 g to 1000 kg, as a multiple of each animal's basal rate. A small animal
has more surface for its mass and less room for fur, and its lower critical temperature is high:

```@example endotherm
fig, ax = figure_axis("Air temperature (°C)", "Metabolic rate required / basal rate")
for mass in (0.1, 1.0, 10.0, 100.0, 1000.0)
    organism = animal(; mass = mass * u"kg")
    lines!(ax, air_temperatures, [required(organism, T) for T in air_temperatures] ./ ustrip(u"W", metabolic_rate(Kleiber(), mass * u"kg"));
           linewidth = 2, label = "$mass kg")
end
hlines!(ax, [1.0]; color = :black, linestyle = :dash)
axislegend(ax; position = :rt)
fig
```

## Compared with NicheMapR

The reference values are from `endoR_devel` of NicheMapR 3.3.3 with thermoregulation turned off
(`THERMOREG = 0`), written by the script `docs/src/data/nichemapr_reference.R`:

```@example endotherm
using DelimitedFiles

table, header = readdlm(joinpath(pkgdir(HeatExchange), "docs", "src", "data", "endoR_reference.csv"), ','; header = true)
column(name, rows) = Float64.(table[rows, findfirst(==(name), vec(header))])
exhaled_at_lung = table[:, 2] .== 100     # DELTAR = 100
exhaled_at_air = table[:, 2] .== 0        # DELTAR = 0, the default of NicheMapR
reference_temperatures = column("air_temperature_C", exhaled_at_lung)
here = [solve(T * u"°C") for T in reference_temperatures]

fig, ax = figure_axis("Air temperature (°C)", "Metabolic rate required (W)")
lines!(ax, reference_temperatures, column("metabolic_W", exhaled_at_lung); linewidth = 2, color = :grey50,
       label = "NicheMapR, air exhaled at lung temperature")
lines!(ax, reference_temperatures, column("metabolic_W", exhaled_at_air); linewidth = 2, color = :grey50, linestyle = :dash,
       label = "NicheMapR default, air exhaled at air temperature")
scatter!(ax, reference_temperatures, [ustrip(u"W", r.energy_flows.metabolic_heat_flow) for r in here]; markersize = 9,
         label = "HeatExchange.jl")
axislegend(ax; position = :rt)
fig
```

```@example endotherm
row = findfirst(==(0.0), reference_temperatures)
markdown_table(["At 0 °C", "HeatExchange.jl", "NicheMapR, exhaled at lung temperature", "NicheMapR default"], [ # hide
    ("Metabolic rate", e.metabolic_heat_flow, column("metabolic_W", exhaled_at_lung)[row] * u"W", column("metabolic_W", exhaled_at_air)[row] * u"W"), # hide
    ("Skin temperature", celsius(t.skin_temperature), column("skin_C", exhaled_at_lung)[row] * u"°C", column("skin_C", exhaled_at_air)[row] * u"°C"), # hide
    ("Fur surface temperature", celsius(t.insulation_temperature), column("fur_C", exhaled_at_lung)[row] * u"°C", column("fur_C", exhaled_at_air)[row] * u"°C"), # hide
    ("Fur conductivity", t.dorsal.insulation_conductivity, column("fur_conductivity", exhaled_at_lung)[row] * u"W/m/K", column("fur_conductivity", exhaled_at_air)[row] * u"W/m/K"), # hide
    ("Convection", e.convection_heat_flow, column("convection_W", exhaled_at_lung)[row] * u"W", column("convection_W", exhaled_at_air)[row] * u"W"), # hide
    ("Evaporation, with respiration", e.evaporation_heat_flow, column("evaporation_W", exhaled_at_lung)[row] * u"W", column("evaporation_W", exhaled_at_air)[row] * u"W"), # hide
    ("Respiratory water loss", m.respiration_mass_flow, column("respiratory_water_g_h", exhaled_at_lung)[row] * u"g/hr", column("respiratory_water_g_h", exhaled_at_air)[row] * u"g/hr"), # hide
    ("Skin area", g.area_skin, column("skin_area_m2", exhaled_at_lung)[row] * u"m^2", column("skin_area_m2", exhaled_at_air)[row] * u"m^2"), # hide
]) # hide
```

Two things are to be seen.

With air exhaled at lung temperature in both, the metabolic rate agrees to within about 1 %, and the surface
temperatures to within a quarter of a degree. The remaining difference is in the solution of the surface: this
package balances the surface energy budget to convergence, where `SIMULSOL` of NicheMapR stops when successive
guesses of the temperatures agree within a tolerance, see [For NicheMapR users](../manual/nichemapr.md).

The default of NicheMapR gives a lower metabolic rate in the cold, the 101.9 W at 0 °C quoted by Kearney et al.
(2021). In NicheMapR air leaves the nose at the air temperature plus an offset, `DELTAR`, zero by default, and so
gives back most of its heat and water on the way out. This package exhales air at lung temperature in this
version, and the animal at 0 °C then loses about eight times as much water, and 10 W more heat, in its breath, see
[Evaporation and respiration](../manual/evaporation_respiration.md). Longwave radiation absorbed and emitted are
also reported differently by the two, though their difference, the net loss, is the same.

## In a natural environment

Outdoors the sky, the ground and the air are at different temperatures, the wind is stronger and there is sun. All
of these are fields of [`EnvironmentalVars`](@ref), and the animal can have a different coat on its back and its
belly. Here is the default animal on a clear, cold night and on a sunny day:

```@example endotherm
night = example_environment_vars(; air_temperature = u"K"(0.0u"°C"), wind_speed = 2.0u"m/s", relative_humidity = 0.6)
night = EnvironmentalVars(; (; (name => getfield(night, name) for name in fieldnames(typeof(night)))...,
                               sky_temperature = u"K"(-25.0u"°C"), ground_temperature = u"K"(-3.0u"°C"))...)
day = example_environment_vars(; air_temperature = u"K"(15.0u"°C"), wind_speed = 2.0u"m/s", relative_humidity = 0.3,
                                 global_radiation = 800.0u"W/m^2", zenith_angle = 30.0u"°")
day = EnvironmentalVars(; (; (name => getfield(day, name) for name in fieldnames(typeof(day)))...,
                             sky_temperature = u"K"(-5.0u"°C"), ground_temperature = u"K"(30.0u"°C"))...)

outdoors = map((night, day)) do environment_vars
    solve_metabolic_rate(mammal, (; environment_pars = example_environment_pars(), environment_vars), u"K"(34.0u"°C"), environment_vars.air_temperature)
end
markdown_table(["", "Cold night", "Sunny day"], [ # hide
    ("Metabolic rate required", outdoors[1].energy_flows.metabolic_heat_flow, outdoors[2].energy_flows.metabolic_heat_flow), # hide
    ("Solar radiation absorbed", outdoors[1].energy_flows.solar_flow, outdoors[2].energy_flows.solar_flow), # hide
    ("Dorsal fur surface", celsius(outdoors[1].thermoregulation.dorsal.insulation_temperature), celsius(outdoors[2].thermoregulation.dorsal.insulation_temperature)), # hide
    ("Ventral fur surface", celsius(outdoors[1].thermoregulation.ventral.insulation_temperature), celsius(outdoors[2].thermoregulation.ventral.insulation_temperature)), # hide
]) # hide
```

On the cold night the back, facing the sky, is colder than the belly. In the sun it is much the warmer, and the
animal is given more heat than it can use. Hourly conditions like these for a whole year come from
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl), see
[Environments and the ecosystem](../manual/ecosystem.md), and Kearney et al. (2021) show the result for an animal
in central Australia.
