# Radiation

An organism absorbs shortwave radiation from the sun, directly and after it has been scattered by the sky and
reflected by the ground. It absorbs longwave radiation emitted by the sky, the ground and the vegetation around it,
and emits longwave radiation from its own surface. Outdoors by day these are usually the largest terms of the heat
budget.

```@setup radiation
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## Solar radiation

[`solar`](@ref) computes ``Q_{sol}`` as the sum of three parts:

```math
\begin{aligned}
Q_{sol,direct} &= \alpha_d \, A_{sil} \, \frac{Q_{direct}}{\cos z} \, (1 - s) \\
Q_{sol,sky} &= \alpha_d \, F_{sky} \, A \, Q_{diffuse} \, (1 - s) \\
Q_{sol,ground} &= \alpha_v \, F_{ground} \, (A - A_{cond}) \, \rho_g \, Q_{global} \, (1 - s)
\end{aligned}
```

| Symbol | Meaning | Where it is set |
|:--|:--|:--|
| ``Q_{global}`` | solar radiation on a horizontal surface | `global_radiation` in [`EnvironmentalVars`](@ref) |
| ``Q_{diffuse}``, ``Q_{direct}`` | its scattered and direct parts | `diffuse_fraction` |
| ``z`` | zenith angle of the sun | `zenith_angle` |
| ``s`` | shade | `shade` |
| ``\rho_g`` | fraction of sunlight reflected by the ground | `ground_albedo` in [`EnvironmentalPars`](@ref) |
| ``\alpha_d``, ``\alpha_v`` | solar absorptivity of the dorsal and ventral surface | `body_absorptivity_dorsal`, `body_absorptivity_ventral` in [`RadiationParameters`](@ref) |
| ``F_{sky}``, ``F_{ground}`` | fraction of the surface that faces the sky and the ground | `sky_view_factor`, `ground_view_factor` |
| ``A``, ``A_{cond}`` | total area, and area in contact with the ground | the body, and `conduction_fraction` |
| ``A_{sil}`` | silhouette area, the area of the shadow cast on a plane normal to the beam | the body, and `solar_orientation` |

The direct beam is measured on a horizontal surface, and dividing by ``\cos z`` gives its strength on a surface
that faces the sun. The silhouette area is from BiophysicalGeometry.jl, and depends on how the body is turned:

```@example radiation
using HeatExchange, BiophysicalGeometry, Unitful
import HeatExchange: Absorptivities, DorsalVentral, ViewFactors, SolarConditions

shape = Cylinder(0.04u"kg", 1000.0u"kg/m^3", 6.0)
body = Body(shape, Naked())
map(orientation -> uconvert(u"cm^2", silhouette(body, orientation)), (NormalToSun(), Intermediate(), ParallelToSun()))
```

The silhouette of this body as the sun sees it, lying on the ground, with the sun overhead, low to its side and
low along its length. Each panel is drawn to its own scale:

```@example radiation
lying = Pose((0.0u"m", 0.0u"m", 0.0u"m"), [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0]) # hide
on_ground = CompositeBody(; parts = (; body), joins = (), root_pose = lying) # hide
fig = Figure(size = (640, 230)) # hide
for (i, (label, direction)) in enumerate(("Sun overhead" => (0.0, 0.0, 1.0), "Sun low, to the side" => (0.0, 1.0, 0.3), # hide
                                          "Sun low, along the body" => (1.0, 0.0, 0.3))) # hide
    ax = Axis(fig[1, i]; title = label, titlesize = 12) # hide
    silhouette_panel!(ax, on_ground, direction) # hide
    hidedecorations!(ax); hidespines!(ax) # hide
end # hide
fig # hide
```

A lizard that turns its side to the sun on a cold morning, and its head to the sun at midday, moves between the
first and the last of these. With the structs that hold the inputs:

```@example radiation
absorptivities = Absorptivities(; body = DorsalVentral(0.85, 0.85), ground = 0.8)
view_factors = ViewFactors(0.4, 0.4, 0.0, 0.0)   # sky, ground, bush, vegetation
sun = SolarConditions(; zenith_angle = 30.0u"°", global_radiation = 800.0u"W/m^2", diffuse_fraction = 0.15, shade = 0.0)

absorbed = map((NormalToSun(), Intermediate(), ParallelToSun())) do orientation
    solar(body, absorptivities, view_factors, sun, silhouette(body, orientation), 0.0u"m^2")
end
markdown_table(["Orientation", "Direct", "From the sky", "From the ground", "Total"], # hide
    [(o, a.solar_direct_flow, a.solar_sky_flow, a.solar_substrate_flow, a.solar_flow) # hide
     for (o, a) in zip(("Normal to the sun", "Intermediate", "Parallel to the sun"), absorbed)]) # hide
```

Here `ground` in `Absorptivities` is the solar absorptivity of the ground, one less its albedo.

Shade reduces all three parts in proportion. In the insulated solvers it also turns the shaded fraction of the
view of the sky into a view of vegetation, which changes the longwave radiation received, see below.

## Longwave radiation

Every surface emits radiation in proportion to the fourth power of its absolute temperature. The organism receives
it from the sky and the ground, [`radiation_in`](@ref):

```math
Q_{IR,in} = \epsilon_d \, F_{sky} \, A \, \epsilon_{sky} \, \sigma T_{sky}^4
          + \epsilon_v \, F_{ground} \, (A - A_{cond}) \, \epsilon_{ground} \, \sigma T_{ground}^4
```

and emits it from its own surface, [`radiation_out`](@ref):

```math
Q_{IR,out} = \epsilon_d \, F_{sky} \, A \, \sigma T_{s}^4 + \epsilon_v \, F_{ground} \, (A - A_{cond}) \, \sigma T_{s}^4
```

where ``\epsilon`` are emissivities, ``\sigma`` is the Stefan–Boltzmann constant and ``T_s`` is the temperature
of the surface. The sky temperature is the temperature of a black body that would emit what the sky does. Under a
clear sky it is well below air temperature, and this is why animals and leaves in the open at night are colder
than the air.

```@example radiation
import HeatExchange: Emissivities, EnvironmentTemperatures

emissivities = Emissivities(; body = DorsalVentral(0.95, 0.95), ground = 1.0, sky = 1.0)
#                                             air, sky, ground, vegetation, bush, substrate
clear_night = EnvironmentTemperatures(u"K"(10.0u"°C"), u"K"(-15.0u"°C"), u"K"(8.0u"°C"), u"K"(10.0u"°C"), u"K"(10.0u"°C"), u"K"(8.0u"°C"))
surface_temperature = u"K"(10.0u"°C")
gained = radiation_in(body, view_factors, emissivities, clear_night).longwave_flow_in
lost = radiation_out(body, view_factors, emissivities, 0.0, surface_temperature, surface_temperature).longwave_flow_out
uconvert(u"W", gained), uconvert(u"W", lost)
```

A lizard at air temperature under this sky loses more longwave radiation than it gains, and cools below the air.

## View factors

A view factor is the fraction of the radiation leaving a surface that reaches another. For an animal in the open,
about half of its surface faces the sky and half the ground. For a single body they are given in
[`RadiationParameters`](@ref):

| Field | Faces | At temperature |
|:--|:--|:--|
| `sky_view_factor` | open sky | `sky_temperature` |
| `ground_view_factor` | ground | `ground_temperature` |
| `bush_view_factor` | vegetation beside and below the organism | `bush_temperature` |
| `vegetation_view_factor` | vegetation overhead | `vegetation_temperature` |

For a body of several parts they are computed for each part from the positions of the parts, by
`silhouette_factors` of BiophysicalGeometry.jl, and the part of a view that is taken up by another part exchanges
radiation with that part, see [Bodies of many parts](multipart.md).

## Through fur

With bare skin, radiation is absorbed and emitted at the skin. With fur or feathers, sunlight is absorbed at the
outer surface of the coat, with an absorptivity of one less the `reflectance` of the fibres. Longwave radiation is
exchanged at the outer surface too by default, and is written in a linear form,

```math
Q_{rad} = \sum_i 4 \, \epsilon \, \sigma \, F_i \, A \left(\frac{T_{rad} + T_i}{2}\right)^3 (T_{rad} - T_i)
```

summed over the sky, the ground, bushes and overhead vegetation, so that it can be solved together with
conduction through the coat. This is within about 1 % of the fourth-power form for the temperature differences
that occur (Kearney et al. 2021). Radiation also carries heat within the coat, from fibre to fibre, which is
included as part of the conductivity of the fur, see [Insulation](insulation.md).
