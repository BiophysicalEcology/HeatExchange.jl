# Convection and conduction

Heat passes between the surface of an organism and the air or water around it by convection, and between the
surface and the ground it touches by conduction.

```@setup convection
using Main.FigureHelpers
using CairoMakie
```

## Convection

Convective heat loss is proportional to the area exposed and to the difference between surface and fluid
temperature:

```math
Q_{conv} = h_c \, A_{conv} \, (T_s - T_a)
```

Everything else is in the heat transfer coefficient ``h_c``, which depends on the size and shape of the body
and the speed and properties of the fluid. [`convection`](@ref) finds it from standard correlations between
dimensionless numbers (Gates 1980, Bird et al. 2002):

```math
h_c = \frac{Nu \, k}{D}
```

where ``Nu`` is the Nusselt number, ``k`` the thermal conductivity of the fluid and ``D`` the characteristic
dimension of the body.

Two processes move the fluid:

- **Forced convection**: wind. ``Nu`` is a function of the Reynolds number ``Re = \rho v D / \mu``.
- **Free convection**: the buoyancy of fluid warmed or cooled by the surface. ``Nu`` is a function of the
  Grashof number, which rises with the temperature difference and with ``D^3``, and of the Prandtl number.

They are combined as (Bird et al. 2002):

```math
Nu = \left(Nu_{free}^3 + Nu_{forced}^3\right)^{1/3}
```

The correlations are chosen by the family of the shape:

| Shape family | Free convection, [`nusselt_free`](@ref) | Forced convection, [`nusselt_forced`](@ref) |
|:--|:--|:--|
| cylinders | piecewise in the Rayleigh number (McAdams 1954, in Kreith 1965) | piecewise in the Reynolds number (McAdams 1954) |
| spheres and ellipsoids | ``2 + 0.6 \, Gr^{1/4} Pr^{1/3}`` (Bird et al. 2002) | ``0.35 \, Re^{0.6}`` (McAdams 1954) |
| plates | ``0.13 \, (Gr \, Pr)^{1/3}`` (Gates 1980) | ``0.032 \, Re^{0.8}`` |
| `DesertIguana`, `LeopardFrog` | as cylinders | as spheres |

These are fitted relations between dimensionless numbers, see
[Units, dimensions and functional traits](units_traits.md#Dimensionless-numbers).

```@example convection
using HeatExchange, BiophysicalGeometry, Unitful

body = Body(Cylinder(1.0u"kg", 1000.0u"kg/m^3", 3.0), Naked())
out = convection(; body, area = total_area(body), air_temperature = u"K"(20.0u"°C"),
                   surface_temperature = u"K"(30.0u"°C"), wind_speed = 1.0u"m/s",
                   atmospheric_pressure = 101325.0u"Pa", fluid = Air())
out.convection_flow, out.heat_transfer_coefficient
```

The coefficient is returned as a [`TransferCoefficients`](@ref), with the combined, free and forced values. In
still air free convection sets the rate, and with wind it soon ceases to matter:

```@example convection
wind_speeds = 10 .^ range(-2, 1; length = 60) .* u"m/s"
coefficients = map(wind_speeds) do wind_speed
    convection(; body, area = total_area(body), air_temperature = u"K"(20.0u"°C"), surface_temperature = u"K"(30.0u"°C"),
                 wind_speed, atmospheric_pressure = 101325.0u"Pa", fluid = Air()).heat_transfer_coefficient
end
fig, ax = figure_axis("Wind speed (m/s)", "Heat transfer coefficient (W m⁻² K⁻¹)"; xscale = log10)
for (field, label) in ((:combined, "combined"), (:free, "free"), (:forced, "forced"))
    lines!(ax, ustrip.(wind_speeds), [ustrip(u"W/m^2/K", getfield(c, field)) for c in coefficients]; linewidth = 2, label)
end
axislegend(ax; position = :lt)
fig
```

### Size

The characteristic dimension is the one length that stands for the size of the body. Because ``h_c`` falls as
``D`` rises, a small animal is coupled closely to air temperature and a large one much less. The formula is a
[`CharacteristicDimFormula`](@ref), set in [`ConvectionParameters`](@ref):

| Formula | Dimension |
|:--|:--|
| [`VolumeCubeRoot`](@ref), the default | the cube root of the volume (Mitchell 1976), plus the depth of the fur |
| [`ScaledDimension`](@ref)`(factor, name)` | `factor` times a named length of the shape, such as `:width_skin` |

```@example convection
leaf = Body(Plate(0.2u"g", 700.0u"kg/m^3", 1.0, 100.0), Naked())
characteristic_dimension(VolumeCubeRoot(), leaf), characteristic_dimension(ScaledDimension(0.7, :width_skin), leaf)
```

For a thin leaf the cube root of the volume says little about how air flows over it, and 0.7 of the width is
used (Gates 1980, Campbell and Norman 1998), see [A leaf](../tutorials/leaf.md).

### Outdoors, and in water

The correlations are for smooth flow in a wind tunnel. Natural wind is turbulent, which raises convection:
`convection_enhancement` in [`EnvironmentalPars`](@ref) multiplies the forced Nusselt number, 1 by default and
about 1.4 outdoors (Kowalski and Mitchell 1976).

`fluid` in [`EnvironmentalPars`](@ref) is [`Air`](@ref) or [`Water`](@ref), with properties from
[FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl). The composition of the air can
be changed with `gas_fractions`, for a burrow or an atmosphere of the past.

### Mass transfer

Water vapour is carried from a wet surface by the same motion of air that carries heat. [`convection`](@ref)
therefore also returns a mass transfer coefficient, from the Sherwood number, which follows from the Nusselt
number by the Chilton–Colburn analogy:

```math
Sh = Nu \left(\frac{Sc}{Pr}\right)^{1/3}, \qquad h_d = \frac{Sh \, D_v}{D}
```

where ``Sc`` is the Schmidt number and ``D_v`` the diffusivity of water vapour in air. It is passed to
[`evaporation`](@ref), so that evaporation and convection are always computed for the same body and wind, see
[Evaporation and respiration](evaporation_respiration.md) and [Flows of mass](gradients.md#Flows-of-mass).

## Conduction

An animal lying on the ground exchanges heat through the area in contact, [`conduction`](@ref):

```math
Q_{cond} = A_{cond} \, \frac{k_{sub}}{x_{sub}} \, (T_s - T_{sub})
```

| Symbol | Meaning | Where it is set |
|:--|:--|:--|
| ``A_{cond}`` | area in contact | `conduction_fraction` in [`ExternalConductionParameters`](@ref), times the total area |
| ``k_{sub}`` | thermal conductivity of the substrate | `substrate_conductivity` in [`EnvironmentalVars`](@ref) |
| ``x_{sub}`` | depth into the substrate over which the temperature difference is taken | `conduction_depth` in [`EnvironmentalPars`](@ref), 2.5 cm by default |
| ``T_{sub}`` | temperature of the substrate | `substrate_temperature` |

```@example convection
conduction(; conduction_area = 0.1 * total_area(body), L = 2.5u"cm", surface_temperature = u"K"(30.0u"°C"),
             substrate_temperature = u"K"(45.0u"°C"), substrate_conductivity = 0.5u"W/m/K")
```

The value is negative: the animal gains heat from the hot ground. The area in contact is taken out of the area
for convection and for radiation to the ground.

With fur, the coat under an animal is pressed flat. Heat passes from the skin through the compressed fur, of
depth `depth_compressed` in [`InsulationParameters`](@ref), to the substrate, as a second path beside the one
through the uncompressed fur to the air, see [Insulation](insulation.md).

## Conduction inside the body

For conduction from the core to the skin, see [Layers as a radial graph](radial_layers.md). For conduction
between the parts of a body, see [Bodies of many parts](multipart.md).
