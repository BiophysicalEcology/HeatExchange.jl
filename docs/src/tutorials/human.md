# A human of many parts

HomoTherm (Kearney et al. 2026) is the human model of NicheMapR: an ellipsoidal head and a cylindrical trunk,
arms and legs, with the endotherm model solved for each part and the parts added together. This tutorial builds
the same person with the functions of [Bodies of many parts](../manual/multipart.md) and compares the result
with HomoTherm, without thermoregulation: the person here does not vasodilate, sweat or let the core
temperature rise. Those responses are added in
[A human that thermoregulates](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/tutorials/human)
in the documentation of BiophysicalBehaviour.jl.

The geometry of this person is built and checked against NicheMapR in the tutorial
[A human: comparison with NicheMapR](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl) of
BiophysicalGeometry.jl. This is its counterpart for heat.

```@setup human
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## The person

The proportions are those of HomoTherm, held in
[BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl): the fraction of the body
mass in each part, and the ratio of length to width of each. The person is 70 kg, with a density of 1050 kg/m³:

```@example human
using HeatExchange, BiophysicalGeometry, Unitful
import BiologicalScaling
import HeatExchange: EnvironmentTemperatures, ViewFactors, AtmosphericConditions, GasFractions

proportions = BiologicalScaling.body_part_proportions(BiologicalScaling.Human())
part_names = (:head, :trunk, :arm, :leg)
mass = NamedTuple{part_names}(Tuple(proportions.mass_fraction .* 70.0u"kg"))
ratio = NamedTuple{part_names}(Tuple(proportions.aspect_ratio))
density = 1050.0u"kg/m^3"

shape(part) = part == :head ? Ellipsoid(mass[part], density, ratio[part], ratio[part]) :
                              Cylinder(mass[part], density, ratio[part])
nothing # hide
```

Each part has its own physiology, the defaults of HomoTherm for a resting person. The limbs have a lower flesh
conductivity than the head and trunk, standing for less blood flow, and a slightly lower core temperature:

```@example human
core_temperature = map(T -> u"K"(T * u"°C"), (head = 36.8, trunk = 36.8, arm = 36.5, leg = 36.7))
flesh_conductivity = map(k -> k * u"W/m/K", (head = 1.1, trunk = 0.9, arm = 0.5, leg = 0.5))
fat_fraction = (head = 0.035, trunk = 0.252, arm = 0.07, leg = 0.161)
sky_view = (head = 0.50, trunk = 0.42, arm = 0.35, leg = 0.35)
ground_view = (head = 0.38, trunk = 0.42, arm = 0.35, leg = 0.35)
bare_skin = (head = 0.6, trunk = 0.0, arm = 0.0, leg = 0.0)      # the face
joined = (head = 0.02666667, trunk = 0.08088128, arm = 0.02, leg = 0.03333333)   # fraction of area in joins
nothing # hide
```

The trunk, arms and legs are clothed in 6 mm of fabric, treated as a layer of very fine fibres. The head has
10 mm of hair on top and none on the face, and is given the mean of 5 mm:

```@example human
clothing = FibreProperties(; diameter = 1.0u"μm", length = 50.0u"mm", density = 3.0e8u"m^-2", depth = 6.0u"mm",
                             reflectance = 0.3, conductivity = 0.209u"W/m/K")
hair = FibreProperties(; diameter = 75.0u"μm", length = 50.0u"mm", density = 3.0e8u"m^-2", depth = 5.0u"mm",
                         reflectance = 0.3, conductivity = 0.209u"W/m/K")
coat = (head = hair, trunk = clothing, arm = clothing, leg = clothing)

body(part) = Body(shape(part), CompositeInsulation(FibrousLayer(coat[part].depth, coat[part].diameter, coat[part].density),
                                                   FatLayer(fat_fraction[part], density)))
nothing # hide
```

## One description for each part

The surface of each part is solved from a description of that part in its surroundings, see
[Bodies of many parts](../manual/multipart.md#Describing-a-part). The room has still air at 50 % relative
humidity, with its walls at air temperature:

```@example human
atmosphere = AtmosphericConditions(0.5, 0.1u"m/s", 101325.0u"Pa")   # relative humidity, wind speed, pressure

function part_setup(part, air_temperature; view_factors = ViewFactors(sky_view[part], ground_view[part], 0.0, 0.0),
                    covered_area = joined[part] * total_area(body(part)))
    T = air_temperature
    return (;
        body = body(part),
        insulation_pars = InsulationParameters(; dorsal = coat[part], ventral = coat[part],
                                                 depth_compressed = coat[part].depth, longwave_depth_fraction = 1.0),
        traits = (; core_temperature = core_temperature[part], flesh_conductivity = flesh_conductivity[part],
                    fat_conductivity = 0.23u"W/m/K", ϵ_body = 0.98, skin_wetness = 0.01, insulation_wetness = 0.0,
                    bare_skin_fraction = bare_skin[part], eye_fraction = 0.0),
        environment_vars = (; temperature = EnvironmentTemperatures(T, T, T, T, T, T), view_factors,
                              atmos = atmosphere, fluid = Air(), solar_flow = 0.0u"W", gas_fractions = GasFractions(),
                              convection_enhancement = 1.0),
        conduction_fraction = 0.0, conductance_coefficient = 0.0u"W/K", ventral_fraction = 0.5,
        longwave_depth_fraction = 1.0,
        covered_area,
        characteristic_dim = characteristic_dimension(VolumeCubeRoot(), body(part)),
    )
end
nothing # hide
```

The area of each part hidden in its joins exchanges no heat with the room. For now it is the fixed fraction
that HomoTherm uses.

## Solving the person

[`solve_coupled_metabolic_rate`](@ref) solves the surface of each part at that part's core temperature, adds
up the heat each must be supplied with, and finds the metabolic rate that supplies it after the heat lost in
breathing. This is what HomoTherm does: it calls `endoR` for each part with thermoregulation and respiration
off, sums the heat generation, and closes the respiration balance once with `ZBRENT_ENDO`. The arms and legs
are entered twice. The oxygen extraction efficiency and the humidity of exhaled air are those of HomoTherm,
and the resting metabolic rate is 105 W:

```@example human
import ModelParameters: stripparams

breathing = stripparams(RespirationParameters(; oxygen_extraction_efficiency = 0.25, exhaled_relative_humidity = 0.92))
resting = 105.0u"W"
order = (:head, :trunk, :arm, :arm, :leg, :leg)

function person(air_temperature; setups = map(part -> part_setup(part, air_temperature), order), kw...)
    solve_coupled_metabolic_rate(;
        part_surface_setups = setups,
        core_temperature = core_temperature.trunk,
        skin_temperature = u"K"(30.0u"°C"), insulation_temperature = air_temperature + 10.0u"K",
        temperature_tolerance = 1e-3u"K",
        respire = true, respiration_pars = breathing, lung_mass = 70.0u"kg",
        air_temperature, atmos = atmosphere, gas_fractions = GasFractions(),
        metabolic_heat_flow_setpoint = resting, resp_tolerance = 1e-5, kw...)
end

cold = person(u"K"(10.0u"°C"))
cold.metabolic_heat_flow
```

At 10 °C in light clothing this person must produce about one and a half times the resting rate. Each part
contributes its own share, at its own skin and clothing temperature:

```@example human
part_of(result) = (head = result.parts[1], trunk = result.parts[2], arm = result.parts[3], leg = result.parts[5])
parts = part_of(cold)
markdown_table(["Part", "Core", "Skin", "Clothing or hair surface", "Heat conducted to the skin"], # hide
    [(string(name), celsius(core_temperature[name]), celsius(p.skin_temperature), celsius(p.insulation_temperature), # hide
      p.net_metabolic) for (name, p) in pairs(parts)]) # hide
```

The arms are thin, and lose the most heat for their mass. Their skin is nonetheless the warmest, because their
core is close beneath it. The trunk has the coolest skin, under its thick layer of fat.

## Through the layers of the trunk

The temperature at each boundary within a part follows from the heat flowing through its layers, see
[Layers as a radial graph](../manual/radial_layers.md#The-temperature-through-the-layers). For the trunk, from
the core through flesh and fat to the skin, and through the clothing to its surface (compare Fig. 3 of Kearney
et al. 2026):

```@example human
import HeatExchange: ThermalConductivities

trunk = body(:trunk)
flesh, fat = core_to_skin_stack(trunk, ThermalConductivities(flesh_conductivity.trunk, 0.23u"W/m/K", nothing))
fat_boundary = core_temperature.trunk - parts.trunk.net_metabolic * stack_resistance((flesh,), trunk)

radial_stack_diagram(["flesh", "fat", "clothing"], # hide
    [flesh_radius(trunk), skin_radius(trunk), insulation_radius(trunk)], # hide
    [core_temperature.trunk, fat_boundary, parts.trunk.skin_temperature, parts.trunk.insulation_temperature]) # hide
```

Fat and clothing each hold a large part of the difference between the core and the air.

## Compared with HomoTherm

The reference values were written by `docs/src/data/nichemapr_reference.R`, which runs `HomoTherm` of NicheMapR
3.3.3 with its defaults at a series of air temperatures, with a wind speed of 0.1 m/s and 50 % relative
humidity:

```@example human
using DelimitedFiles

data(file) = readdlm(joinpath(pkgdir(HeatExchange), "docs", "src", "data", file), ','; header = true)
whole, whole_header = data("homotherm_whole.csv")
by_part, part_header = data("homotherm_parts.csv")
column(table, header, name) = table[:, findfirst(==(name), vec(header))]

reference_temperature = Float64.(column(whole, whole_header, "air_temperature_C"))
reference_rate = Float64.(column(whole, whole_header, "metabolic_W"))
cold_range = reference_temperature .<= 18.0     # above this, HomoTherm begins to thermoregulate

results = [person(u"K"(T * u"°C")) for T in reference_temperature[cold_range]]
rates = [ustrip(u"W", r.metabolic_heat_flow) for r in results]

fig, ax = figure_axis("Air temperature (°C)", "Metabolic rate (W)")
lines!(ax, reference_temperature, reference_rate; linewidth = 2, color = :grey50, label = "HomoTherm")
scatter!(ax, reference_temperature[cold_range], rates; markersize = 9, label = "HeatExchange.jl")
hlines!(ax, [ustrip(u"W", resting)]; color = :black, linestyle = :dash)
axislegend(ax; position = :rt)
fig
```

The dashed line is the resting metabolic rate. Below about 19 °C the person is below the lower critical
temperature, and must raise the metabolic rate to stay warm. Above it, the rate required falls below the
resting rate, and HomoTherm begins to dilate the blood vessels of the skin, and then to sweat. That is where
this tutorial stops.

The two agree to about 1 % for the whole person:

```@example human
maximum(abs.(rates ./ reference_rate[cold_range] .- 1))
```

and part by part, at 10 °C:

```@example human
at_10 = Float64.(column(by_part, part_header, "air_temperature_C")) .== 10.0
reference(part, name) = Float64(only(column(by_part, part_header, name)[at_10 .& (column(by_part, part_header, "part") .== string(part))]))
mean_of(part, a, b) = (reference(part, a) + reference(part, b)) / 2

markdown_table(["Part", "Heat, here", "Heat, HomoTherm", "Skin, here", "Skin, HomoTherm"], # hide
    [(string(name), p.net_metabolic, reference(name, "heat_generated_W") * u"W", celsius(p.skin_temperature), # hide
      mean_of(name, "skin_dorsal_C", "skin_ventral_C") * u"°C") for (name, p) in pairs(parts)]) # hide
```

The differences have known causes:

- **Joins.** HomoTherm removes the joined area of a part by treating it as lying on a substrate at core
  temperature. Here it is taken out of the area that exchanges heat. The trunk, with the most joined area,
  differs the most.
- **Sides.** HomoTherm solves each part as a dorsal and a ventral side and averages them, so the top of the
  head has hair and the face none. Here each part has one surface, and the head a coat of the mean depth.
- **Breathing.** HomoTherm takes the lung temperature from the temperature profile of the trunk, and exhales
  air at a temperature that depends on the air temperature. Here the lung temperature is the mean of the core
  and mean skin temperatures, and air is exhaled at it.

## A person with a place in space

In HomoTherm the parts have no positions. The fraction of each part that is joined, and the fraction of its
view that is sky and ground, are numbers given to the model. With a `CompositeBody` of BiophysicalGeometry.jl
the parts are joined at named places, and those numbers follow from where the parts are:

```@example human
trunk, head, arm, leg = body(:trunk), body(:head), body(:arm), body(:leg)
trunk_length = trunk.geometry.length.length_skin
r_arm, r_leg = skin_radius(arm), skin_radius(leg)
shoulder(side) = Attachment(Lateral(trunk_length - r_arm, side), Disc(r_arm))
hip(side) = Attachment(EndA(1.1r_leg, side), Disc(r_leg))
leg_top = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))

human = CompositeBody(;
    parts = (; trunk, head, arm_left = arm, arm_right = arm, leg_left = leg, leg_right = leg),
    joins = (
        Join(trunk = Attachment(EndB(0.0u"m", 0.0), Disc(r_arm)), head = Attachment(PoleB(), Disc(r_arm))),
        Join(trunk = shoulder(0.0), arm_left = Attachment(Lateral(r_arm, π), Disc(r_arm)); twist = π),
        Join(trunk = shoulder(π), arm_right = Attachment(Lateral(r_arm, 0.0), Disc(r_arm)); twist = π),
        Join(trunk = hip(0.0), leg_left = leg_top),
        Join(trunk = hip(π), leg_right = leg_top),
    ),
)
names = keys(human.parts)
kind = (trunk = :trunk, head = :head, arm_left = :arm, arm_right = :arm, leg_left = :leg, leg_right = :leg)
nothing # hide
```

```@example human
composite_views(human; views = (:oblique, :front, :side), titles = ["", "side", "front"], size = (700, 340)) # hide
```

```@example human
body_graph!(Axis(Figure(size = (520, 300))[1, 1]), human); current_figure() # hide
```

The area of each part that is covered is the sum of the areas of its joins:

```@example human
covered = map(names) do name
    sum(join_area(join, human) for join in human.joins if name in join_partners(join))
end
map(area -> uconvert(u"cm^2", area), NamedTuple{names}(covered))
```

and each part's view of the sky, the ground and each other part comes from `silhouette_factors`:

```@example human
views = silhouette_factors(human, Sky(0.5))
markdown_table(["Part", "Sky", "Ground", "Other parts"], # hide
               [(string(name), v.sky, v.ground, sum(v.neighbours)) for (name, v) in pairs(views)]) # hide
```

A fifth of what an arm sees is the trunk beside it. That part of its view exchanges longwave radiation with
the clothing of the trunk, which is warmer than the walls. The `neighbour_topology` gives, for each part, the
parts it sees and how much of its view each takes, see
[Bodies of many parts](../manual/multipart.md#Parts-that-see-each-other):

```@example human
index = NamedTuple{names}(Tuple(1:length(names)))
topology = map(names) do name
    Tuple((; index = index[other], fraction) for (other, fraction) in pairs(views[name].neighbours))
end

placed(air_temperature) = person(air_temperature;
    setups = map(names, covered) do name, covered_area
        part_setup(kind[name], air_temperature; covered_area,
                   view_factors = ViewFactors(views[name].sky, views[name].ground, 0.0, 0.0))
    end,
    neighbour_topology = topology)

in_place = placed(u"K"(10.0u"°C"))
in_place.metabolic_heat_flow, cold.metabolic_heat_flow
```

The skin and clothing temperatures of the person in place, at 10 °C:

```@example human
skin = NamedTuple{names}(Tuple(ustrip(u"°C", part.skin_temperature) for part in in_place.parts)) # hide
temperature_views(human, skin; views = (:oblique, :front), titles = ["", "skin"]) # hide
```

```@example human
surface = NamedTuple{names}(Tuple(ustrip(u"°C", part.insulation_temperature) for part in in_place.parts)) # hide
temperature_views(human, surface; views = (:oblique, :front), titles = ["", "clothing and hair"], # hide
                  label = "Surface temperature (°C)", colormap = :viridis) # hide
```

```@example human
placed_rates = [ustrip(u"W", placed(u"K"(T * u"°C")).metabolic_heat_flow) for T in reference_temperature[cold_range]]

fig, ax = figure_axis("Air temperature (°C)", "Metabolic rate (W)")
lines!(ax, reference_temperature[cold_range], reference_rate[cold_range]; linewidth = 2, color = :grey50, label = "HomoTherm")
scatter!(ax, reference_temperature[cold_range], rates; markersize = 9, label = "parts as given to HomoTherm")
scatter!(ax, reference_temperature[cold_range], placed_rates; markersize = 9, marker = :utriangle, label = "parts in place")
hlines!(ax, [ustrip(u"W", resting)]; color = :black, linestyle = :dash)
axislegend(ax; position = :rt)
fig
```

The joins made here hide less area than the fractions of HomoTherm, which were chosen so that the share of the
whole area in each part matches that measured on people. And the parts now face each other in place of the
walls over part of their view. The two changes work in opposite directions, and in this room they nearly
cancel.

The value of placing the parts is not in this comparison, in a room with walls at one temperature. It is that
the same description can be turned. With the arms held out, the person curled up, or the sun to one side, the
joined areas, the views and the shade of one part on another all change, and HomoTherm has no way to know it.

## What is left out

This person holds every core temperature fixed and has one state. A real person in the cold constricts the
blood vessels of the limbs and lets them cool, and in the heat dilates them, sweats, and lets the core
temperature rise, and HomoTherm does all of these (Kearney et al. 2026). With the functions of this package:

- cooler limbs are a change of `core_temperature` and `flesh_conductivity` for those parts, or a limb that is
  its own compartment joined to the trunk, see
  [Bodies of many parts](../manual/multipart.md#A-torso-and-four-legs);
- sweating is a change of `skin_wetness`.

Choosing those changes is the work of BiophysicalBehaviour.jl, see
[A human that thermoregulates](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/tutorials/human).
