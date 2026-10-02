# Metabolism

Metabolic heat is the one term of the heat budget that the organism generates. For an ectotherm it is small and
follows from body temperature. For an endotherm it is large, and is the quantity that the organism adjusts to
hold its temperature.

```@setup metabolism
using Main.FigureHelpers
using CairoMakie
```

!!! note "Moving to BiologicalScaling.jl"
    The metabolic rate equations on this page are defined in HeatExchange.jl in this version. They are allometric
    equations, and their home is [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl),
    which already has a wider set (for example `standard_metabolic_rate(Squamate(), mass, temperature)`). A coming
    version of this package will take its metabolic rates from there, and the types below will be replaced by
    those of BiologicalScaling.jl. The role of metabolism in the heat budget, described here, will not change.
    These equations are descriptive and not physical, and `allometric` is there to mark that, see
    [Units, dimensions and functional traits](units_traits.md).

## The parameters

[`MetabolismParameters`](@ref) holds:

| Field | Meaning | Used when |
|:--|:--|:--|
| `model` | an equation for metabolic rate from mass and temperature | solving for temperature |
| `core_temperature` | the core temperature to be held | solving for metabolic rate |
| `metabolic_heat_flow` | the minimum metabolic rate, basal or resting | solving for metabolic rate |
| `q10` | the factor by which the metabolic rate changes for 10 °C of core temperature | by thermoregulation, when the core temperature is allowed to rise |

When the temperature is the unknown, [`solve_temperature`](@ref) calls
[`metabolic_rate`](@ref)`(model, mass, core_temperature)` at each trial temperature. When the metabolic rate is the
unknown, no equation is used for it. `metabolic_heat_flow` is then the floor below which a solution means that the
animal cannot lose its heat, and the rate at which it still breathes, see
[Temperature or metabolic rate](solvers.md).

## The equations

Each is a subtype of [`MetabolicRateEquation`](@ref), called with [`metabolic_rate`](@ref):

| Equation | For | Source |
|:--|:--|:--|
| [`AndrewsPough2`](@ref) | squamate reptiles, from mass and body temperature, standard or resting | Andrews and Pough (1985), their equation 2 |
| [`Kleiber`](@ref) | basal rate of mammals, ``3.39 \, M^{0.75}`` W with ``M`` in kg | Kleiber (1947) |
| [`McKechnieWolf`](@ref) | basal rate of birds | McKechnie and Wolf (2004) |
| [`PlantDarkRespiration`](@ref) | dark respiration of leaves, from mass and temperature | Reich et al. (2006) |
| `nothing` | no metabolic heat | |

```@example metabolism
using HeatExchange, Unitful

metabolic_rate(Kleiber(), 65.0u"kg"), metabolic_rate(McKechnieWolf(), 30.0u"g"),
metabolic_rate(AndrewsPough2(), 40.0u"g", u"K"(30.0u"°C"))
```

The rate of an ectotherm is far below that of an endotherm of the same mass, and it rises steeply with
temperature:

```@example metabolism
temperatures = 5.0:1.0:45.0   # °C
fig, ax = figure_axis("Body temperature (°C)", "Metabolic rate (mW)")
for (state, label) in ((0.0, "standard"), (1.0, "resting"))
    model = AndrewsPough2(; metabolic_state = state)
    lines!(ax, temperatures, [ustrip(u"mW", metabolic_rate(model, 40.0u"g", u"K"(T * u"°C"))) for T in temperatures];
           linewidth = 2, label)
end
axislegend(ax; position = :lt)
fig
```

The equation of Andrews and Pough is

```math
\dot{V}_{O_2} = M_1 \, m^{M_2} \, 10^{M_3 T_b} \, 10^{M_4}
```

in ml of oxygen per hour, with mass ``m`` in g and body temperature ``T_b`` in °C, held between 1 and 50 °C.
The four constants are the fields `mass_normalisation`, `mass_exponent`, `thermal_sensitivity` and
`metabolic_state` of [`AndrewsPough2`](@ref), with a `metabolic_state` of 0 for the standard rate and 1 for the
resting rate. Another equation is added by defining a subtype of [`MetabolicRateEquation`](@ref) and a method of
[`metabolic_rate`](@ref) for it.

## Oxygen and heat

Metabolic rate is measured as oxygen consumed and is needed here as heat produced. [`O2_to_Joules`](@ref) and
[`Joules_to_O2`](@ref) convert between them with an [`OxygenJoulesConversion`](@ref):

| Conversion | Energy per volume of oxygen |
|:--|:--|
| [`Typical`](@ref) | 20.1 J/ml |
| [`Kleiber1961`](@ref) | 5.0, 4.5 or 4.7 kcal/l for a respiratory quotient of 1 or more, between 0.7 and 1, or 0.7 or less (Kleiber 1961) |

```@example metabolism
oxygen = Joules_to_O2(Typical(), 100.0u"W", 0.8)
uconvert(u"ml/hr", oxygen), uconvert(u"ml/hr", Joules_to_O2(Kleiber1961(), 100.0u"W", 0.8))
```

The oxygen consumed sets the air breathed, and so the heat and water lost in respiration, see
[Evaporation and respiration](evaporation_respiration.md).

## Metabolism and temperature

The relation of metabolic rate to temperature can be given in more detail than a single equation allows, with
thermal performance curves that fall away at high temperature, from
[ThermalPhysiology.jl](https://github.com/BiophysicalEcology/ThermalPhysiology.jl). And the metabolic rate of an
animal over its life, as it grows and reproduces, is the subject of Dynamic Energy Budget theory, which the
ectotherm model of NicheMapR includes (Kearney and Porter 2020) and which is outside this package.
