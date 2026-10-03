# A leaf

A leaf has a heat budget like any other organism. It absorbs sunlight and longwave radiation, loses heat by
longwave radiation and convection, and cools itself by evaporating water. What differs is how the water leaves:
through stomata that open and close, on a surface otherwise nearly sealed. NicheMapR treats a leaf as an option
of its ectotherm model (`leaf = 1`). Here a leaf is an [`Organism`](@ref) whose evaporation parameters are
[`LeafEvaporationParameters`](@ref), and everything else is shared with animals.

```@setup leaf
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## The leaf

The body is a thin square plate, 0.2 g with the density of leaf tissue. The second ratio makes its sides 100
times its thickness:

```@example leaf
using HeatExchange, BiophysicalGeometry, Unitful

shape = Plate(0.2u"g", 700.0u"kg/m^3", 1.0, 100.0)
body = Body(shape, Naked())
map(x -> uconvert(u"mm", x), body.geometry.length), uconvert(u"cm^2", total_area(body))
```

```@example leaf
fig = Figure(size = (420, 300)) # hide
draw_parts!(body_axis(fig[1, 1]), body; colors = [RGBf(0.35, 0.62, 0.30)]) # hide
fig # hide
```

Its traits are built by position, with the leaf parameters in place of those of an animal:

```@example leaf
leaf_traits(; abaxial = 0.3u"mol/m^2/s", adaxial = 0.0u"mol/m^2/s", absorptivity = 0.5) = HeatExchangeTraits(
    shape,
    InsulationParameters(),
    ExternalConductionParameters(; conduction_fraction = 0.0),
    InternalConductionParameters(; flesh_conductivity = 100.0u"W/m/K"),
    RadiationParameters(; body_absorptivity_dorsal = absorptivity, body_absorptivity_ventral = absorptivity,
                          body_emissivity_dorsal = 0.97, body_emissivity_ventral = 0.97,
                          sky_view_factor = 0.5, ground_view_factor = 0.5, solar_orientation = Intermediate()),
    ConvectionParameters(; characteristic_dimension_formula = ScaledDimension(0.7, :width_skin)),
    LeafEvaporationParameters(; abaxial_vapour_conductance = abaxial, adaxial_vapour_conductance = adaxial,
                                cuticular_conductance = 0.01u"mol/m^2/s"),
    HydraulicParameters(; water_potential = 0.0u"J/kg"),
    RespirationParameters(),
    MetabolismParameters(; metabolic_heat_flow = 0.0u"W", model = PlantDarkRespiration()),
    SolveMetabolicRateOptions(),
)
leaf = Organism(body, leaf_traits())
nothing # hide
```

Four choices make this a leaf:

| Choice | Reason |
|:--|:--|
| [`LeafEvaporationParameters`](@ref) | water leaves through stomata, with a vapour conductance for the lower (abaxial) and upper (adaxial) surface, and a small cuticular conductance that remains when they are closed, see [Evaporation and respiration](../manual/evaporation_respiration.md#From-a-leaf) |
| [`ScaledDimension`](@ref)`(0.7, :width_skin)` | the boundary layer of a flat leaf is set by 0.7 of its width (Campbell and Norman 1998), not by the cube root of its volume, see [Convection and conduction](../manual/convection_conduction.md#Size) |
| a very high `flesh_conductivity` | a leaf is too thin to have a temperature difference between its inside and its surface |
| [`PlantDarkRespiration`](@ref) | the metabolic heat of a leaf is its dark respiration, negligible in the heat budget |

## Leaf temperature

On a sunny day with a light wind:

```@example leaf
conditions(; air_temperature = 25.0u"°C", wind_speed = 1.0u"m/s", relative_humidity = 0.4, global_radiation = 800.0u"W/m^2") =
    (; environment_pars = example_environment_pars(),
       environment_vars = EnvironmentalVars(;
           air_temperature = u"K"(air_temperature), sky_temperature = u"K"(10.0u"°C"),
           ground_temperature = u"K"(35.0u"°C"), substrate_temperature = u"K"(air_temperature),
           relative_humidity, wind_speed, atmospheric_pressure = 101325.0u"Pa", zenith_angle = 30.0u"°",
           substrate_conductivity = 0.5u"W/m/K", global_radiation, diffuse_fraction = 0.15, shade = 0.0))

out = solve_temperature(leaf, conditions())
u"°C"(out.core_temperature)
```

```@example leaf
flow_table(out.energy_balance) # hide
```

Evaporation, here transpiration, is a large term for a leaf, as it never is for a lizard. The water transpired
is in the mass balance:

```@example leaf
uconvert(u"mmol/m^2/s", out.mass_balance.transpiration_mass / 18.0u"g/mol" / total_area(body))
```

## Stomata

When the stomata close, the leaf loses its cooling and warms:

```@example leaf
conductances = 0.0:0.02:0.6   # mol m⁻² s⁻¹
results = [solve_temperature(Organism(body, leaf_traits(; abaxial = g * u"mol/m^2/s")), conditions()) for g in conductances]

fig = Figure(size = (760, 320))
ax1 = Axis(fig[1, 1]; xlabel = "Stomatal conductance (mol m⁻² s⁻¹)", ylabel = "Leaf temperature (°C)")
ax2 = Axis(fig[1, 2]; xlabel = "Stomatal conductance (mol m⁻² s⁻¹)", ylabel = "Transpiration (mmol m⁻² s⁻¹)")
lines!(ax1, conductances, [ustrip(u"°C", r.core_temperature) for r in results]; linewidth = 2)
hlines!(ax1, [25.0]; color = :black, linestyle = :dash)
lines!(ax2, conductances, [ustrip(u"mmol/m^2/s", r.mass_balance.transpiration_mass / 18.0u"g/mol" / total_area(body)) for r in results];
       linewidth = 2)
fig
```

The dashed line is the air temperature. Transpiration does not rise in proportion to stomatal conductance, for
two reasons that the heat budget supplies:

- the stomata are in series with the boundary layer, which limits the flow when they are wide open;
- a leaf that transpires more is cooler, which lowers the vapour density at its surface.

This feedback between leaf temperature and transpiration is why the two must be solved together. In the terms
of [Gradients, resistances and flows](../manual/gradients.md), the stomata are a resistor that the plant
controls.

A plant under water stress closes its stomata, and the effect on leaf temperature is one of the ways drought
harms it. Stomatal conductance is an input here. A model of how stomata respond to light, humidity and soil
water supplies it, and the leaf temperature found here feeds back to that model.

## Wind and leaf size

The boundary layer is thinner on a small leaf and in a strong wind, and a leaf with a thin boundary layer is
held close to air temperature:

```@example leaf
wind_speeds = 10 .^ range(-1, 1; length = 30)

fig, ax = figure_axis("Wind speed (m/s)", "Leaf temperature − air temperature (°C)"; xscale = log10)
for mass in (0.02u"g", 0.2u"g", 2.0u"g")
    sized = Plate(mass, 700.0u"kg/m^3", 1.0, 100.0)
    width = uconvert(u"cm", Body(sized, Naked()).geometry.length.width_skin)
    traits = HeatExchangeTraits(sized, InsulationParameters(), conduction_pars_external(leaf), conduction_pars_internal(leaf),
        radiation_pars(leaf), convection_pars(leaf), evaporation_pars(leaf), hydraulic_pars(leaf), respiration_pars(leaf),
        metabolism_pars(leaf), SolveMetabolicRateOptions())
    organism = Organism(Body(sized, Naked()), traits)
    excess = [ustrip(u"K", solve_temperature(organism, conditions(; wind_speed = v * u"m/s")).core_temperature - u"K"(25.0u"°C"))
              for v in wind_speeds]
    lines!(ax, wind_speeds, excess; linewidth = 2, label = "$(round(u"cm", width; digits = 1)) wide")
end
hlines!(ax, [0.0]; color = :black, linestyle = :dash)
axislegend(ax; position = :rt)
fig
```

Large leaves in still air run well above air temperature in the sun, one reason that the leaves of plants of
hot, dry places are small.

## Thick leaves and stems

The same organism with a realistic tissue conductivity and a shape with some depth is a succulent leaf or a
cactus stem, in which the surface in the sun is hotter than the inside. Nothing else needs to change.

## Other leaf models

The leaf temperature vignette of NicheMapR compares this calculation with other leaf energy balance models, and
couples it to a model of photosynthesis and stomatal conductance. Those are outside this package.
