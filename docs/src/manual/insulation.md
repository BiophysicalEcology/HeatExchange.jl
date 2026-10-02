# Insulation

Fur and feathers are porous: fibres with still air between them. Heat moves through the layer by conduction
along the fibres and through the air, and by radiation between the fibres, all at once (Porter et al. 1994,
Kearney et al. 2021). This page describes how a coat is specified and how its conductivity is found. "Fur"
here covers feathers, hair and clothing too.

```@setup insulation
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## Specifying a coat

A coat is described in two places, which must agree.

The **body** has a `FibrousLayer` from BiophysicalGeometry.jl, with a depth, a fibre diameter and a fibre
density. From these come the outer radius and area of the animal, and the area of skin left bare between the
fibres.

The **traits** have an [`InsulationParameters`](@ref), holding a [`FibreProperties`](@ref) for the dorsal and
for the ventral surface:

| Field of [`FibreProperties`](@ref) | Meaning | NicheMapR |
|:--|:--|:--|
| `depth` | depth of the coat, from skin to outer surface | `ZFURD`, `ZFURV` |
| `length` | length of a fibre. Fibres longer than the coat is deep lie at an angle | `LHAIRD`, `LHAIRV` |
| `diameter` | diameter of a fibre | `DHAIRD`, `DHAIRV` |
| `density` | number of fibres per area of skin | `RHOD`, `RHOV` |
| `reflectance` | fraction of sunlight reflected | `REFLD`, `REFLV` |
| `conductivity` | thermal conductivity of the fibre material, 0.209 W m⁻¹ K⁻¹ for keratin | `KHAIR` |

and two values for the whole coat:

| Field of [`InsulationParameters`](@ref) | Meaning | NicheMapR |
|:--|:--|:--|
| `depth_compressed` | depth of the coat where the animal lies on it | `ZFURCOMP` |
| `longwave_depth_fraction` | fraction of the depth at which longwave radiation is exchanged, 1 for the outer surface | `XR` |

```@example insulation
using HeatExchange, BiophysicalGeometry, Unitful

fibres = FibreProperties(; diameter = 30.0u"μm", length = 23.9u"mm", density = 3000.0u"cm^-2", depth = 9.0u"mm",
                           reflectance = 0.3, conductivity = 0.209u"W/m/K")
insulation_pars = InsulationParameters(; dorsal = fibres, ventral = fibres, depth_compressed = 9.0u"mm",
                                         longwave_depth_fraction = 1.0)
nothing # hide
```

The body takes the depth, diameter and density from the same values. Where dorsal and ventral coats differ, the
body is given their mean, weighted by the fraction of the surface that is ventral, as in NicheMapR:

```@example insulation
fur = FibrousLayer(fibres.depth, fibres.diameter, fibres.density)
body = Body(Ellipsoid(10.0u"kg", 1000.0u"kg/m^3", 2.7, 2.7), CompositeInsulation(fur, FatLayer(0.0, 901.0u"kg/m^3")))
nothing # hide
```

```@example insulation
shape_gallery("" => body; size = (420, 300)) # hide
```

The coat drawn to scale, with the fraction of the skin that the bases of its fibres cover:

```@example insulation
plot_insulation_properties(fur)
```

## Conduction through the coat

[`insulation_properties`](@ref) computes what the heat budget needs of the coat, by the theory of Conley and
Porter (1986). It corresponds to `IRPROP` and `GETKFUR` of NicheMapR.

```@example insulation
properties = insulation_properties(insulation_pars, u"K"(20.0u"°C"), 0.3)   # temperature in the coat, ventral fraction
properties.conductivities.dorsal
```

This *effective conductivity* is that of fibres and air together. The fibres are treated as a regular array of
cylinders. Heat is conducted in parallel through fibre and air along the fibres, and in series across them, and
the effective conductivity is the mean of the two, held between that of air and that of the fibre. The air is
taken to be still, and its conductivity rises with temperature, so the result depends on the temperatures in
the coat. The fibres conduct much better than air, so anything that puts more fibre into the layer raises the
conductivity:

```@example insulation
air_conductivity = 0.0257u"W/m/K"   # at 20 °C
vary(; kw...) = ustrip(u"W/m/K", insulation_thermal_conductivity(FibreProperties(; diameter = fibres.diameter,
    length = fibres.length, density = fibres.density, depth = fibres.depth, reflectance = 0.3,
    conductivity = fibres.conductivity, kw...), air_conductivity).effective_conductivity)

fig = Figure(size = (760, 300))
diameters = 2.0:2.0:100.0
densities = 500.0:500.0:35000.0
depths = 1.0:1.0:50.0
for (i, (x, y, label)) in enumerate((
        (diameters, [vary(diameter = d * u"μm") for d in diameters], "Fibre diameter (μm)"),
        (densities, [vary(density = d * u"cm^-2") for d in densities], "Fibre density (cm⁻²)"),
        (depths, [vary(depth = d * u"mm") for d in depths], "Coat depth (mm)")))
    ax = Axis(fig[1, i]; xlabel = label, ylabel = i == 1 ? "Effective conductivity (W m⁻¹ K⁻¹)" : "")
    lines!(ax, x, y; linewidth = 2)
    hlines!(ax, [ustrip(u"W/m/K", air_conductivity)]; color = :grey50, linestyle = :dash)
end
fig
```

The dashed line is still air. The ranges stop short of the point at which the fibres touch. A deeper coat of the
same fibres has a lower conductivity, because fibres of fixed length stand more upright and less of each lies
within a given depth.

[`insulation_properties`](@ref) returns an [`InsulationProperties`](@ref), with values for the dorsal surface,
the ventral surface and their average, each as a [`BodyRegionValues`](@ref):

| Field | Content |
|:--|:--|
| `fibres` | the fibre properties |
| `conductivities` | effective conductivities |
| `absorption_coefficients`, `optical_thickness` | how far longwave radiation penetrates the coat |
| `conductivity_compressed` | effective conductivity of the coat pressed to `depth_compressed` |
| `insulation_test` | zero if there is no coat, which sends the solvers down the bare-skin path |

See [The endotherm, piece by piece](../tutorials/components.md) for these functions called one at a time.

## Radiation within the coat

Each fibre radiates to its neighbours, and where there is a temperature gradient this carries heat down it. The
effect is an added conductivity, which rises with the cube of the temperature and falls as the fibres become
denser (Conley and Porter 1986):

```math
k_{rad} = \frac{16 \, \sigma \, T^3}{3 \, \beta}, \qquad \beta = \frac{0.67}{\pi} \, \rho_{eff} \, d
```

where ``\beta`` is the absorption coefficient of the coat, ``\rho_{eff}`` the density of fibres allowing for
their angle and ``d`` their diameter. The conductivity of the fur is the sum

```math
k_{fur} = k_{eff} + k_{rad}
```

It depends on the skin and fur surface temperatures, which are being solved for, and is recomputed at each
step. The value at the solution is in the output as `insulation_conductivity`, beside the effective
conductivity alone. Measured conductivities of fur include both parts.

## The surface of the coat

Heat conducted through a coat on a cylinder of length ``L``, from the skin at radius ``R_s`` to the outer
surface at ``R_{fa}``, is

```math
Q_{fur} = \frac{2 \pi \, k_{fur} \, L \, (T_s - T_{fa})}{\ln(R_{fa} / R_s)}
```

with corresponding forms for a sphere and an ellipsoid. At the outer surface this heat, with the sunlight
absorbed there, leaves by convection, longwave radiation and the evaporation of any water on the fur. These are
the equations of [`solve_part_heat_balance`](@ref), and the skin and surface temperatures are found together,
see [Solving a heat balance](heat_balance.md#With-insulation). The method is that of Mathewson and Porter
(2013), with the surface balance solved to convergence by a Newton iteration.

[`radiant_temperature`](@ref) gives the temperature at which the coat exchanges longwave radiation. With a
`longwave_depth_fraction` of 1 it is that of the outer surface. A smaller fraction places the exchange within
the coat. That option is inherited from NicheMapR and is not yet reliable here: the closed form is singular as
the depth approaches the surface, see
[Layers as a radial graph](radial_layers.md#Where-the-chain-is-not-a-chain).

## Lying on the ground

Where an animal lies on the ground its coat is pressed flat. The fraction in contact, `conduction_fraction`,
loses heat by a second path: from the skin through the compressed coat, of depth `depth_compressed` and
conductivity `conductivity_compressed`, to the substrate. The two paths are in parallel, and
[`compressed_radiant_temperature`](@ref) and [`insulation_radiant_temperature`](@ref) hold the solutions for
each shape. It is not a small correction: for a 1 kg animal under a cold sky on warm ground, a contact of 30 %
of the surface carried about a quarter of the heat exchanged, see
[Layers as a radial graph](radial_layers.md#Where-the-chain-is-not-a-chain).

```@example insulation
radial_network_diagram() # hide
```

## Back and belly

Dorsal and ventral coats can differ in every property. The solvers treat the animal as two sides, each with its
own coat and surroundings: the dorsal side sees the sky and any vegetation overhead and takes the direct
sunlight, and the ventral side sees the ground and takes the sunlight reflected from it. `ventral_fraction` in
[`RadiationParameters`](@ref) is the fraction of the surface that is ventral. The output has the temperatures,
the conductivity and every heat flow for each side:

```@example insulation
thin = FibreProperties(; diameter = 30.0u"μm", length = 23.9u"mm", density = 3000.0u"cm^-2", depth = 3.0u"mm",
                         reflectance = 0.3, conductivity = 0.209u"W/m/K")
traits = example_heat_exchange_traits(;
    shape_pars = body.shape,
    insulation_pars = InsulationParameters(; dorsal = fibres, ventral = thin, depth_compressed = 3.0u"mm",
                                             longwave_depth_fraction = 1.0),
    metabolism_pars = example_metabolism_pars(; metabolic_heat_flow = metabolic_rate(Kleiber(), 10.0u"kg")),
)
mean_fur = FibrousLayer(6.0u"mm", 30.0u"μm", 3000.0u"cm^-2")
animal = Organism(Body(body.shape, CompositeInsulation(mean_fur, FatLayer(0.0, 901.0u"kg/m^3"))), traits)
environment = (; environment_pars = example_environment_pars(),
                 environment_vars = example_environment_vars(; air_temperature = u"K"(5.0u"°C")))
out = solve_metabolic_rate(animal, environment, u"K"(34.0u"°C"), u"K"(5.0u"°C"))

markdown_table(["", "Dorsal", "Ventral"], [ # hide
    ("Coat depth", out.thermoregulation.dorsal.insulation_depth, out.thermoregulation.ventral.insulation_depth), # hide
    ("Skin temperature", celsius(out.thermoregulation.dorsal.skin_temperature), celsius(out.thermoregulation.ventral.skin_temperature)), # hide
    ("Fur surface temperature", celsius(out.thermoregulation.dorsal.insulation_temperature), celsius(out.thermoregulation.ventral.insulation_temperature)), # hide
    ("Fur conductivity", out.thermoregulation.dorsal.insulation_conductivity, out.thermoregulation.ventral.insulation_conductivity), # hide
    ("Heat conducted to the skin", out.energy_flows.dorsal.net_metabolic, out.energy_flows.ventral.net_metabolic), # hide
]) # hide
```

The heat conducted to the skin is given for each side as if the whole animal were that side, and the two are
averaged by the view of each. For a body in which back and belly are separate parts, see
[Back and belly: two halves](../tutorials/two_parts.md) and [Bodies of many parts](multipart.md).

## Fat

Fat under the skin is a second layer of insulation, inside the skin where fur is outside it. It is given to the
body as a `FatLayer`, a fraction of the body mass at a density, and to the traits as `fat_fraction`,
`fat_density` and `fat_conductivity` in [`InternalConductionParameters`](@ref). See
[Layers as a radial graph](radial_layers.md).

## Changing the coat

An animal raises and flattens its coat. That is a change to `depth` followed by another solve, and is the
piloerection response of BiophysicalBehaviour.jl, see
[Endotherm thermoregulation by rules](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/endotherm_rules).
