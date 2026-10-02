# Bodies of many parts

A real animal is not one shape. Its limbs, ears and tail are thin and lose heat fast, its back and belly have
different coats and face different surroundings, and it can be warm at the core and cold at the feet. Kearney et
al. (2021) showed how the endotherm model of NicheMapR could be arranged for several body parts, and HomoTherm
(Kearney et al. 2026) did so for a human. This page describes the functions here that solve a heat budget for a
body of several parts.

The parts themselves, their joins, their areas and what each one sees, come from a `CompositeBody` of
[BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl). The functions on this page
take one description for each part, and are the layer beneath
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), which builds those
descriptions from a `CompositeBody` and solves the whole animal. They can also be called directly, as here.

To build a body of several parts, see the documentation of BiophysicalGeometry.jl: its manual pages on bodies as
graphs and on joins, its tutorials, which build a dog, a cow and a human, and its interactive page
[Build an animal](https://biophysicalecology.github.io/BiophysicalGeometry.jl/dev/builder), where an animal is
assembled with sliders and the code that builds it can be copied. The tutorial
[A human of many parts](../tutorials/human.md) here takes such a body through to its heat budget.

```@setup multipart
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## The scheme

The heat budget of a body of several parts is solved in three steps, the first two of which are the two sides of
[Solving a heat balance](heat_balance.md) for one part:

1. **Each part's surface.** For a given core temperature of a part, find its skin and fur surface temperatures
   and the heat ``Q_{gen,net}`` that must be conducted from its core to hold them. This is
   [`solve_part_surface`](@ref).
2. **The cores.** Parts exchange heat with each other inside the body. Find the core temperature of each part that
   is not regulated, from the heat it makes, the heat it loses through its own surface, and the heat it is given
   by its neighbours. This is [`solve_regulated_core_temperatures`](@ref).
3. **The whole animal.** Add up the heat that the regulated core must supply, and find the metabolic rate that
   supplies it after the heat lost in breathing. This is [`solve_coupled_metabolic_rate`](@ref).

## Describing a part

A part is described by a NamedTuple with what its surface solve needs. This function builds one, for a furred
shape in a still, shaded environment:

```@example multipart
using HeatExchange, BiophysicalGeometry, Unitful
import HeatExchange: EnvironmentTemperatures, ViewFactors, AtmosphericConditions

vars = example_environment_vars(; air_temperature = u"K"(5.0u"°C"), wind_speed = 1.0u"m/s")
pars = example_environment_pars()

function part_setup(shape, fibres, core_temperature; covered_area = 0.0u"m^2",
                    view_factors = ViewFactors(0.5, 0.5, 0.0, 0.0))   # sky, ground, bush, vegetation
    fur = FibrousLayer(fibres.depth, fibres.diameter, fibres.density)
    body = Body(shape, CompositeInsulation(fur, FatLayer(0.0, 901.0u"kg/m^3")))
    return (;
        body,
        insulation_pars = InsulationParameters(; dorsal = fibres, ventral = fibres, depth_compressed = fibres.depth,
                                                 longwave_depth_fraction = 1.0),
        traits = (; core_temperature, flesh_conductivity = 0.9u"W/m/K", fat_conductivity = 0.23u"W/m/K",
                    ϵ_body = 0.99, skin_wetness = 0.005, insulation_wetness = 0.0, bare_skin_fraction = 0.0,
                    eye_fraction = 0.0),
        environment_vars = (;
            temperature = EnvironmentTemperatures(vars), view_factors, atmos = AtmosphericConditions(vars),
            fluid = pars.fluid, solar_flow = 0.0u"W", gas_fractions = pars.gas_fractions,
            convection_enhancement = pars.convection_enhancement,
        ),
        conduction_fraction = 0.0, conductance_coefficient = 0.0u"W/K", ventral_fraction = 0.5,
        longwave_depth_fraction = 1.0,
        covered_area,
        characteristic_dim = characteristic_dimension(VolumeCubeRoot(), body),
    )
end
nothing # hide
```

| Field | Content |
|:--|:--|
| `body` | the part, with its own layers |
| `insulation_pars` | its coat, see [Insulation](insulation.md). A part has one coat, given as both `dorsal` and `ventral` |
| `traits` | its core temperature, conductivities, emissivity and wetness |
| `environment_vars` | its surroundings: temperatures, **its own view factors**, the air, and **the sunlight it absorbs**, in W |
| `conduction_fraction`, `conductance_coefficient` | its contact with the ground |
| `covered_area` | the area hidden by its joins to other parts, which exchanges no heat with the environment |
| `characteristic_dim` | the length for convection, see [Convection and conduction](convection_conduction.md) |

Each part has its own coat and its own exposure, and the two are independent: a part can have thick fur and face
the ground, or thin fur and face the sky. The view factors and the absorbed sunlight are where the geometry of
the whole body enters. For a `CompositeBody` they come from `silhouette_factors` and `silhouette` of
BiophysicalGeometry.jl, which account for parts shading each other.

## One part's surface

```@example multipart
core_temperature = u"K"(37.0u"°C")
thick = FibreProperties(; diameter = 30.0u"μm", length = 25.0u"mm", density = 3000.0u"cm^-2", depth = 15.0u"mm",
                          reflectance = 0.2, conductivity = 0.209u"W/m/K")
thin = FibreProperties(; diameter = 30.0u"μm", length = 25.0u"mm", density = 3000.0u"cm^-2", depth = 3.0u"mm",
                         reflectance = 0.2, conductivity = 0.209u"W/m/K")
density = 1000.0u"kg/m^3"

torso_shape = Cylinder(8.0u"kg", density, 2.5)
leg_shape = Cylinder(0.5u"kg", density, 6.0)
join_area = π * skin_radius(Body(leg_shape, Naked()))^2      # the top of a leg, against the torso

guess = (; skin_temperature = u"K"(33.0u"°C"), insulation_temperature = u"K"(10.0u"°C"), temperature_tolerance = 1e-3u"K")
torso = part_setup(torso_shape, thick, core_temperature; covered_area = 4 * join_area)
torso_surface = solve_part_surface(; torso..., guess...)
u"°C"(torso_surface.skin_temperature), u"°C"(torso_surface.insulation_temperature), torso_surface.net_metabolic
```

The torso and its four legs as a `CompositeBody`, with its graph of parts and joins. The joins are discs of the
radius of a leg, which is the `join_area` used above:

```@example multipart
torso_body, leg_body = torso.body, part_setup(leg_shape, thin, core_temperature).body # hide
torso_length = torso_body.geometry.length.length_skin # hide
hip(z, φ) = Attachment(Lateral(z * torso_length, φ), Disc(skin_radius(leg_body))) # hide
leg_top = Attachment(EndA(0.0u"m", 0.0), Disc(skin_radius(leg_body))) # hide
lying = Pose((0.0u"m", 0.0u"m", 0.0u"m"), [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0]) # hide
animal = CompositeBody(; # hide
    parts = (; torso = torso_body, leg_fl = leg_body, leg_fr = leg_body, leg_bl = leg_body, leg_br = leg_body), # hide
    joins = (Join(torso = hip(0.8, -π / 2 + 0.35), leg_fl = leg_top), Join(torso = hip(0.8, -π / 2 - 0.35), leg_fr = leg_top), # hide
             Join(torso = hip(0.2, -π / 2 + 0.35), leg_bl = leg_top), Join(torso = hip(0.2, -π / 2 - 0.35), leg_br = leg_top)), # hide
    root_pose = lying) # hide
body_graph(animal) # hide
```

[`solve_part_surface`](@ref) also returns `flesh_conductance`, the heat conducted per degree between the core and
the skin of the part, which the next step uses.

## Compartments and couplings

How heat moves between two joined parts is a [`HeatCoupling`](@ref):

| Coupling | Meaning |
|:--|:--|
| [`SharedCore`](@ref)`()` | the two parts have one core temperature. Blood mixes freely between them |
| [`ConductiveCoupling`](@ref)`()` | heat is conducted through the flesh of each part to the join, over the area of the join |
| [`ConductiveCoupling`](@ref)`(h)` | heat crosses the join with a given conductance per area, `h` in W m⁻² K⁻¹ |
| none | the join is insulated |

Parts joined by [`SharedCore`](@ref) form one *compartment*, with one core temperature.
[`compartment_graph`](@ref) finds the compartments from the names of the parts and the pairs that share a core:

```@example multipart
graph = compartment_graph((:dorsal, :ventral, :head, :leg), ((:dorsal, :ventral),))
num_compartments(graph), parts_in_compartment(graph, 1), compartment_of(graph, :leg)
```

```@example multipart
compartment_diagram(graph; edges = ((:dorsal, :ventral, :shared), (:dorsal, :head, :conductive), (:ventral, :leg, :conductive)), # hide
    positions = Dict(:dorsal => (0.0, 1.0), :ventral => (0.0, 0.0), :head => (2.0, 1.0), :leg => (2.0, 0.0))) # hide
```

The thick line is a shared core and the dashed lines are conductive couplings. The number in each part is its
compartment. The back and the belly of an animal are the usual case of a shared core, and the tutorial
[Back and belly: two halves](../tutorials/two_parts.md) solves it.

The conductance of a conductive join follows from its area, the distance from the centre of each part to the
join, and the conductivity of the flesh of each part, with [`contribution_to_conductance`](@ref). The area and the
distances are `join_area` and `internal_distance` of BiophysicalGeometry.jl:

```@example multipart
leg_length = Body(leg_shape, Naked()).geometry.length.length_skin
by_conduction = contribution_to_conductance(ConductiveCoupling(), join_area, skin_radius(torso.body), leg_length / 2,
                                            0.9u"W/m/K", 0.9u"W/m/K")
with_blood_flow = contribution_to_conductance(ConductiveCoupling(300.0u"W/m^2/K"), join_area, skin_radius(torso.body),
                                              leg_length / 2, 0.9u"W/m/K", 0.9u"W/m/K")
by_conduction, with_blood_flow
```

Conduction through flesh alone carries very little heat along a limb. In a living animal most of the heat that
reaches a limb is carried by blood, and a conductance per area is a simple way to stand for that. A coupling for
blood flow in its own right is planned, and is added by defining a new subtype of [`HeatCoupling`](@ref) with its
contribution to the conductance and to the heat load of each compartment, [`contribution_to_heat_load`](@ref).

## The core temperatures

At the core of each compartment, the heat made there, less respiration, plus heat from any coupling, equals the
heat conducted to the skins of its own parts and to the cores of its neighbours:

```math
Q_{gen,c} - Q_{resp,c} = \sum_{p \in c} G_p \, (T_c - T_{s,p}) + \sum_j G_{cj} \, (T_c - T_j)
```

where ``G_p`` is the `flesh_conductance` of part ``p`` and ``G_{cj}`` the conductance of the join to compartment
``j``. With the skin temperatures held at their current values this is a set of linear equations in the core
temperatures, one for each compartment. [`build_conductance_matrix`](@ref) assembles the matrix of the
conductances between compartments, and [`solve_core_temperatures`](@ref) solves the system:

```@example multipart
two = compartment_graph((:a, :b), ())
solve_core_temperatures(two, ((1, 2, 0.5u"W/K"),),     # the join between compartments 1 and 2
    (5.0u"W", 1.0u"W"),                               # heat made in each, less respiration
    (1.0u"W/K", 1.0u"W/K"),                           # flesh conductance of each
    (1.0u"W/K" * 300.0u"K", 1.0u"W/K" * 300.0u"K"))   # flesh conductance times skin temperature
```

Both skins are at 300 K. The first compartment makes 5 W, loses 4 W through its own skin and passes 1 W to the
second, which makes 1 W and loses 2 W.

In an endotherm one compartment is regulated: its core is held at the set temperature, and the heat it must make
is what is being solved for. [`solve_regulated_core_temperatures`](@ref) holds that core fixed and solves for the
others, which float to whatever balances their own heat budgets.

## A torso and four legs

Here the torso is regulated at 37 °C and each leg is a compartment of its own, joined to the torso with the
conductance for blood flow above. The legs have thin fur and make a little heat of their own. Their surface
temperatures depend on their core temperature and their core temperature on their surface, so the two steps are
repeated until the core temperature of a leg stops changing:

```@example multipart
function solve_leg(conductance; leg_generation = 0.3u"W")
    legs = compartment_graph((:torso, :leg), ())
    leg_core = core_temperature
    leg_surface = nothing
    for _ in 1:500
        leg = part_setup(leg_shape, thin, leg_core; covered_area = join_area)
        leg_surface = solve_part_surface(; leg..., guess...)
        cores = solve_regulated_core_temperatures(legs, 1, ((1, 2, conductance),),
            (0.0u"W", leg_generation),
            (torso_surface.flesh_conductance, leg_surface.flesh_conductance),
            (torso_surface.flesh_conductance * torso_surface.skin_temperature,
             leg_surface.flesh_conductance * leg_surface.skin_temperature),
            core_temperature)
        abs(cores[2] - leg_core) < 1e-4u"K" && break
        leg_core = cores[2]
    end
    return (; core_temperature = leg_core, surface = leg_surface, heat_from_torso = conductance * (core_temperature - leg_core))
end

leg = solve_leg(with_blood_flow)
u"°C"(leg.core_temperature), u"°C"(leg.surface.skin_temperature), leg.heat_from_torso
```

```@example multipart
leg_skin = ustrip(u"°C", leg.surface.skin_temperature) # hide
temperature_views(animal, (torso = ustrip(u"°C", torso_surface.skin_temperature), leg_fl = leg_skin, leg_fr = leg_skin, # hide
                           leg_bl = leg_skin, leg_br = leg_skin); views = (:oblique, :side)) # hide
```

The legs are much cooler than the torso. The heat that the torso passes to them is added to what it must supply to
its own surface, as `extra_net_metabolic`, and the metabolic rate of the animal follows:

```@example multipart
whole_animal(setups; extra = 0.0u"W") = solve_coupled_metabolic_rate(;
    part_surface_setups = setups, core_temperature, guess...,
    respire = true, respiration_pars = example_respiration_pars(), lung_mass = 10.0u"kg",
    air_temperature = vars.air_temperature, atmos = AtmosphericConditions(vars), gas_fractions = pars.gas_fractions,
    metabolic_heat_flow_setpoint = metabolic_rate(Kleiber(), 10.0u"kg"), resp_tolerance = 1e-5,
    extra_net_metabolic = extra)

cool_legs = whole_animal((torso,); extra = 4 * leg.heat_from_torso)
cool_legs.metabolic_heat_flow
```

Compare this with the two limits. If the legs were held at the core temperature of the torso, as parts of one
compartment, each would lose heat as fast as its thin fur allows. If they were joined to the torso by conduction
through flesh alone, they would be given almost nothing and would fall nearly to air temperature:

```@example multipart
warm_leg = part_setup(leg_shape, thin, core_temperature; covered_area = join_area)
warm_legs = whole_animal((torso, warm_leg, warm_leg, warm_leg, warm_leg))
cold_leg = solve_leg(by_conduction)
cold_legs = whole_animal((torso,); extra = 4 * cold_leg.heat_from_torso)

markdown_table(["Legs", "Leg core temperature", "Heat to each leg", "Metabolic rate"], [ # hide
    ("one compartment with the torso", celsius(core_temperature), warm_legs.parts[2].net_metabolic, warm_legs.metabolic_heat_flow), # hide
    ("joined with blood flow", celsius(leg.core_temperature), leg.heat_from_torso, cool_legs.metabolic_heat_flow), # hide
    ("joined by conduction only", celsius(cold_leg.core_temperature), cold_leg.heat_from_torso, cold_legs.metabolic_heat_flow), # hide
]) # hide
```

The basal metabolic rate of an animal of this mass is

```@example multipart
metabolic_rate(Kleiber(), 10.0u"kg")
```

so letting the legs cool is the difference between a large cost of keeping warm and a small one. Animals of cold
climates do exactly this.

## Parts that see each other

A part whose view is partly blocked by another part exchanges longwave radiation with that part's surface over
the blocked fraction, in place of the sky or the ground behind it. The surfaces of the parts then depend on each
other, and [`solve_coupled_metabolic_rate`](@ref) solves them together when given a `neighbour_topology`: for each
part, the parts it sees and the fraction of its view that each takes up. The fractions are those returned by
`silhouette_factors` of BiophysicalGeometry.jl, in which the views of the sky, the ground and the other parts sum
to one.

```@example multipart
sees = (((; index = 2, fraction = 0.1),),    # the torso sees the leg over a tenth of its view
        ((; index = 1, fraction = 0.3),))    # the leg sees the torso over three tenths of its view
torso_seen = part_setup(torso_shape, thick, core_temperature; covered_area = join_area, view_factors = ViewFactors(0.45, 0.45, 0.0, 0.0))
leg_seen = part_setup(leg_shape, thin, core_temperature; covered_area = join_area, view_factors = ViewFactors(0.35, 0.35, 0.0, 0.0))
together = solve_coupled_metabolic_rate(;
    part_surface_setups = (torso_seen, leg_seen), core_temperature, guess..., respire = false,
    respiration_pars = example_respiration_pars(), lung_mass = 10.0u"kg", air_temperature = vars.air_temperature,
    atmos = AtmosphericConditions(vars), gas_fractions = pars.gas_fractions,
    metabolic_heat_flow_setpoint = metabolic_rate(Kleiber(), 10.0u"kg"), resp_tolerance = 1e-5,
    neighbour_topology = sees)
map(part -> u"°C"(part.insulation_temperature), together.parts)
```

A warm neighbour is a warmer thing to face than a cold sky, so parts that are close together lose less heat than
the same parts apart. This is why limbs held against the body, and animals huddled together, save energy.

## Half shapes

The back and the belly of one shape are two parts, each a half: `HalfCylinder`, `HalfEllipsoid` and `HalfSphere`
of BiophysicalGeometry.jl. Every heat-transfer method for a shape family works for its halves, and the fur of a
half covers half of the shell. See the tutorial [Back and belly: two halves](../tutorials/two_parts.md).

## What is and is not here

| | State |
|:--|:--|
| Surface solve of each part, with its own coat, exposure and core temperature | [`solve_part_surface`](@ref) |
| Parts sharing a core, summed to a metabolic rate | [`solve_coupled_metabolic_rate`](@ref) |
| Compartments with their own core temperatures, coupled by conduction | [`solve_regulated_core_temperatures`](@ref), with the iteration shown above |
| Longwave exchange between parts | `neighbour_topology` |
| The residual form of each part for an optimiser | [`part_surface_residuals`](@ref), see [Differentiability and the NLP interface](autodiff.md) |
| Building the parts from a `CompositeBody`, the iteration over compartments, and thermoregulation | BiophysicalBehaviour.jl |
| Heat carried between parts by blood, with counter-current exchange | planned, as a [`HeatCoupling`](@ref) |
| Solving for the temperature of a body of several parts | planned. [`solve_coupled_metabolic_rate`](@ref) solves for metabolic rate |
