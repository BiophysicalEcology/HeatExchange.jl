# Back and belly: two halves

The back of an animal faces the sun and sky, and its belly the ground, often with thinner fur. There are two
ways to put this into a heat budget. The solvers of this package, like NicheMapR, solve one body twice, as if
it were all back and then all belly, and average. Or the body can be two real halves, each with its own coat
and surroundings, sharing one core. This tutorial does both and compares them. It is the companion of the
tutorial of the same name in
[BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl), which builds the
halves.

```@setup two
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## One body, two sides

The animal is a furred cylinder of 1 kg. With [`solve_metabolic_rate`](@ref) it is one body, and the dorsal and
ventral coats are the two [`FibreProperties`](@ref) of its [`InsulationParameters`](@ref), see
[Insulation](../manual/insulation.md#Back-and-belly):

```@example two
using HeatExchange, BiophysicalGeometry, Unitful
import HeatExchange: EnvironmentTemperatures, ViewFactors, AtmosphericConditions, Absorptivities, SolarConditions

mass, density, ratio = 1.0u"kg", 1000.0u"kg/m^3", 3.0
coat(depth) = FibreProperties(; diameter = 30.0u"μm", length = 23.9u"mm", density = 3000.0u"cm^-2", depth,
                                reflectance = 0.2, conductivity = 0.209u"W/m/K")
fat = FatLayer(0.0, 901.0u"kg/m^3")
whole_shape = Cylinder(mass, density, ratio)
metabolism = example_metabolism_pars(; metabolic_heat_flow = metabolic_rate(Kleiber(), mass))
core_temperature = metabolism.core_temperature

function one_body(back, belly, environment)
    mean_depth = (back.depth + belly.depth) / 2
    body = Body(whole_shape, CompositeInsulation(FibrousLayer(mean_depth, back.diameter, back.density), fat))
    traits = example_heat_exchange_traits(; shape_pars = whole_shape, metabolism_pars = metabolism,
        insulation_pars = InsulationParameters(; dorsal = back, ventral = belly, depth_compressed = belly.depth,
                                                 longwave_depth_fraction = 1.0))
    return solve_metabolic_rate(Organism(body, traits), environment, core_temperature - 5.0u"K",
                                environment.environment_vars.air_temperature + 2.0u"K")
end
nothing # hide
```

## Two halves, one core

As two parts, the animal is two `HalfCylinder`s of half the mass. Each has the radius and length of the whole
cylinder and half its volume and curved surface. The flat face of each is against the other, and is its
`covered_area`. The halves have the outline of the whole cylinder, so the length that governs convection is
that of the whole:

```@example two
function half(fibres, view_factors, solar_flow, environment)
    vars, pars = environment.environment_vars, environment.environment_pars
    body = Body(HalfCylinder(mass / 2, density, ratio), CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density), fat))
    whole = Body(whole_shape, CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density), fat))
    return (;
        body,
        insulation_pars = InsulationParameters(; dorsal = fibres, ventral = fibres, depth_compressed = fibres.depth,
                                                 longwave_depth_fraction = 1.0),
        traits = (; core_temperature, flesh_conductivity = 0.9u"W/m/K", fat_conductivity = 0.23u"W/m/K", ϵ_body = 0.99,
                    skin_wetness = 0.005, insulation_wetness = 0.0, bare_skin_fraction = 0.0, eye_fraction = 0.0),
        environment_vars = (; temperature = EnvironmentTemperatures(vars), view_factors, atmos = AtmosphericConditions(vars),
                              fluid = pars.fluid, solar_flow, gas_fractions = pars.gas_fractions,
                              convection_enhancement = pars.convection_enhancement),
        conduction_fraction = 0.0, conductance_coefficient = 0.0u"W/K", ventral_fraction = 0.5, longwave_depth_fraction = 1.0,
        covered_area = 2 * body.geometry.length.radius_skin * body.geometry.length.length_skin,
        characteristic_dim = characteristic_dimension(VolumeCubeRoot(), whole),
    )
end
nothing # hide
```

The dorsal half sees only sky and the ventral half only ground. Sunlight is divided the same way: the direct
beam and the light scattered from the sky fall on the back, and the light reflected from the ground on the
belly:

```@example two
function two_halves(back, belly, environment)
    vars, pars = environment.environment_vars, environment.environment_pars
    body = Body(whole_shape, CompositeInsulation(FibrousLayer((back.depth + belly.depth) / 2, back.diameter, back.density), fat))
    sun = solar(body, Absorptivities(example_radiation_pars(), pars), ViewFactors(0.5, 0.5, 0.0, 0.0), SolarConditions(vars),
                silhouette(body, Intermediate()), 0.0u"m^2")
    dorsal = half(back, ViewFactors(1.0, 0.0, 0.0, 0.0), sun.solar_direct_flow + sun.solar_sky_flow, environment)
    ventral = half(belly, ViewFactors(0.0, 1.0, 0.0, 0.0), sun.solar_substrate_flow, environment)
    return solve_coupled_metabolic_rate(;
        part_surface_setups = (dorsal, ventral), core_temperature,
        skin_temperature = core_temperature - 5.0u"K", insulation_temperature = vars.air_temperature + 2.0u"K",
        temperature_tolerance = 1e-3u"K", respire = true, respiration_pars = example_respiration_pars(), lung_mass = mass,
        air_temperature = vars.air_temperature, atmos = AtmosphericConditions(vars), gas_fractions = pars.gas_fractions,
        metabolic_heat_flow_setpoint = metabolism.metabolic_heat_flow, resp_tolerance = 1e-5)
end
nothing # hide
```

The two halves are one compartment: they are joined by a [`SharedCore`](@ref), and the heat the core must
supply is the sum of what each half conducts to its skin, see
[Bodies of many parts](../manual/multipart.md#Compartments-and-couplings).

## The symmetric case

With the same coat on both halves, sky and ground at air temperature and no sun, nothing distinguishes back
from belly, and the two formulations must agree:

```@example two
indoors = (; environment_pars = example_environment_pars(), environment_vars = example_environment_vars())
even = coat(2.0u"mm")
a = one_body(even, even, indoors)
b = two_halves(even, even, indoors)
a.energy_flows.metabolic_heat_flow, b.metabolic_heat_flow
```

```@example two
b.metabolic_heat_flow / a.energy_flows.metabolic_heat_flow - 1
```

This agreement is a test of the package. It holds because a half cylinder is treated exactly as half of a
whole one: its exposed area leaves out the flat face, its fur covers half of the shell, and its convection is
that of the whole outline.

## Different coats, different surroundings

Now the back has 10 mm of fur and the belly 3 mm, and the animal is outdoors on a clear, cold night, and then
in the sun:

```@example two
function outdoors(; air = 5.0u"°C", sky = -20.0u"°C", ground = 2.0u"°C", sun = 0.0u"W/m^2")
    base = example_environment_vars(; air_temperature = u"K"(air), wind_speed = 1.0u"m/s", global_radiation = sun, zenith_angle = 30.0u"°")
    vars = EnvironmentalVars(; (; (name => getfield(base, name) for name in fieldnames(typeof(base)))...,
                                  sky_temperature = u"K"(sky), ground_temperature = u"K"(ground))...)
    return (; environment_pars = example_environment_pars(), environment_vars = vars)
end

back, belly = coat(10.0u"mm"), coat(3.0u"mm")
cases = ("Indoors, even coat" => (even, even, indoors), "Indoors" => (back, belly, indoors),
         "Clear cold night" => (back, belly, outdoors()),
         "Sunny day" => (back, belly, outdoors(; air = 15.0u"°C", sky = -5.0u"°C", ground = 30.0u"°C", sun = 600.0u"W/m^2")))

compared = map(cases) do (label, (dorsal_coat, ventral_coat, environment))
    sides = one_body(dorsal_coat, ventral_coat, environment)
    halves = two_halves(dorsal_coat, ventral_coat, environment)
    (label, sides.energy_flows.metabolic_heat_flow, halves.metabolic_heat_flow,
     celsius(sides.thermoregulation.dorsal.insulation_temperature), celsius(halves.parts[1].insulation_temperature),
     celsius(sides.thermoregulation.ventral.insulation_temperature), celsius(halves.parts[2].insulation_temperature))
end
markdown_table(["", "Metabolic rate, two sides", "two halves", "Back surface, two sides", "two halves", # hide
                "Belly surface, two sides", "two halves"], collect(compared)) # hide
```

The two agree in every case, to a tenth of a percent, with different coats, a cold sky and sun. Two sides, each
solved as a whole animal and averaged, are the same model as two halves that share a core, as long as the back
sees only sky and the belly only ground. The halves are a generalisation of the two sides, and the method of
NicheMapR is recovered exactly.

## What each half really sees

The assumption that the back sees only sky is the weak point of both. The flanks of the dorsal half face
outwards, and see the ground as well. With the halves joined into a `CompositeBody`, BiophysicalGeometry.jl
computes what each sees:

```@example two
dorsal_body = Body(HalfCylinder(mass / 2, density, ratio), CompositeInsulation(FibrousLayer(back.depth, back.diameter, back.density), fat))
ventral_body = Body(HalfCylinder(mass / 2, density, ratio), CompositeInsulation(FibrousLayer(belly.depth, belly.diameter, belly.density), fat))
lying = Pose((0.0u"m", 0.0u"m", 0.0u"m"), [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0])   # the long axis horizontal
composite = CompositeBody(;
    parts = (; dorsal = dorsal_body, ventral = ventral_body),
    joins = (Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover()); twist = -π / 2),),
    root_pose = lying,
)
views = silhouette_factors(composite, Sky(0.5))
(; dorsal = (; views.dorsal.sky, views.dorsal.ground), ventral = (; views.ventral.sky, views.ventral.ground))
```

```@example two
composite_views(composite; views = (:oblique, :side, :front)) # hide
```

The back has the deeper fur, so its half is the wider. Seen by the sun, overhead and low to the side:

```@example two
fig = Figure(size = (520, 230)) # hide
for (i, (label, direction)) in enumerate(("Sun overhead" => (0.0, 0.0, 1.0), "Sun low, to the side" => (0.0, 1.0, 0.3))) # hide
    ax = Axis(fig[1, i]; title = label, titlesize = 12) # hide
    silhouette_panel!(ax, composite, direction) # hide
    hidedecorations!(ax); hidespines!(ax) # hide
end # hide
fig # hide
```

These are fractions of the whole surface of each half, and the rest is its flat face, against the other half.
As fractions of the exposed surface:

```@example two
exposed(v) = ViewFactors(v.sky / (v.sky + v.ground), v.ground / (v.sky + v.ground), 0.0, 0.0)
exposed(views.dorsal), exposed(views.ventral)
```

Nearly a third of what the back sees is ground, and a fifth of what the belly sees is sky. On the cold night,
with these views in place of all sky and all ground:

```@example two
function two_halves_in_place(back, belly, environment)
    vars, pars = environment.environment_vars, environment.environment_pars
    dorsal = half(back, exposed(views.dorsal), 0.0u"W", environment)
    ventral = half(belly, exposed(views.ventral), 0.0u"W", environment)
    return solve_coupled_metabolic_rate(;
        part_surface_setups = (dorsal, ventral), core_temperature,
        skin_temperature = core_temperature - 5.0u"K", insulation_temperature = vars.air_temperature + 2.0u"K",
        temperature_tolerance = 1e-3u"K", respire = true, respiration_pars = example_respiration_pars(), lung_mass = mass,
        air_temperature = vars.air_temperature, atmos = AtmosphericConditions(vars), gas_fractions = pars.gas_fractions,
        metabolic_heat_flow_setpoint = metabolism.metabolic_heat_flow, resp_tolerance = 1e-5)
end

night = outdoors()
assumed = two_halves(back, belly, night)
in_place = two_halves_in_place(back, belly, night)
markdown_table(["Views", "Metabolic rate", "Back surface", "Belly surface"], [ # hide
    ("back sees sky, belly sees ground", assumed.metabolic_heat_flow, celsius(assumed.parts[1].insulation_temperature), celsius(assumed.parts[2].insulation_temperature)), # hide
    ("computed from the geometry", in_place.metabolic_heat_flow, celsius(in_place.parts[1].insulation_temperature), celsius(in_place.parts[2].insulation_temperature)), # hide
]) # hide
```

```@example two
temperature_views(composite, (dorsal = ustrip(u"°C", in_place.parts[1].insulation_temperature), # hide
                              ventral = ustrip(u"°C", in_place.parts[2].insulation_temperature)); # hide
                  views = (:oblique, :front), label = "Fur surface temperature (°C)") # hide
```

The back is warmer and the belly colder than the simple division gives, because each sees some of what the
other does. The difference grows with the contrast between sky and ground.

## What two halves allow

The comparison above is on ground both formulations can stand on. The halves go further:

- Views, and shade from other parts, computed from the positions of the parts of a `CompositeBody`, with a
  head, limbs and tail joined to either half.
- One half on the ground while the other is in the wind, with `conduction_fraction` set for that half alone.
- Its own skin wetness and flesh conductivity for each half: an animal that sweats from its back, or sends
  blood to a bare belly.
- A half as a compartment of its own, joined by a [`ConductiveCoupling`](@ref) in place of a shared core: a
  shell or a hump.

These are assembled from a `CompositeBody` by
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), see
[Bodies of many parts](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/multipart) in
its documentation. The next tutorial, [A human of many parts](human.md), builds a body of six parts by hand.
