# Solving a heat balance

The heat budget that the package solves, what it means to solve it, and how the solution is found: first for
bare skin, then for fur or feathers.

```@setup balance
using Main.FigureHelpers
using CairoMakie
```

## The budget

At steady state the heat an organism gains equals the heat it loses (Porter and Gates 1969):

```math
Q_{sol} + Q_{IR,in} + Q_{met} = Q_{IR,out} + Q_{conv} + Q_{cond} + Q_{evap} + Q_{resp}
```

```@example balance
heat_budget_diagram() # hide
```

| Term | Heat flow | Function | Depends on |
|:--|:--|:--|:--|
| ``Q_{sol}`` | solar radiation absorbed | [`solar`](@ref) | sunlight, silhouette and total area, absorptivity |
| ``Q_{IR,in}`` | longwave radiation absorbed | [`radiation_in`](@ref) | sky and ground temperature, view factors, emissivity |
| ``Q_{met}`` | metabolic heat produced | [`metabolic_rate`](@ref) | mass, core temperature |
| ``Q_{IR,out}`` | longwave radiation emitted | [`radiation_out`](@ref) | **surface temperature**, area, emissivity |
| ``Q_{conv}`` | convection to the air | [`convection`](@ref) | **surface temperature**, air temperature, wind, size and shape |
| ``Q_{cond}`` | conduction to the substrate | [`conduction`](@ref) | **surface temperature**, substrate temperature, area of contact |
| ``Q_{evap}`` | evaporation from the skin and eyes | [`evaporation`](@ref) | **surface temperature**, humidity, wet area |
| ``Q_{resp}`` | heat lost in breathing | [`respiration`](@ref) | metabolic rate, **lung temperature**, humidity |

Terms on the right are losses, and a negative value is a gain: convection from warmer air, conduction from a
warmer rock, condensation. See [Radiation](radiation.md),
[Convection and conduction](convection_conduction.md),
[Evaporation and respiration](evaporation_respiration.md) and [Metabolism](metabolism.md).

The first two terms are set by the environment. Every other term depends on the temperature of the organism,
which is why the budget must be solved and not just added up.

Many of the terms are flows down a difference in potential against a resistance, and the rest are sources, see
[Gradients, resistances and flows](gradients.md).

## From the core to the surface

Heat is exchanged at the surface, but the temperature that matters is that of the core. Metabolic heat is made
throughout the flesh and conducted outwards, so the surface is cooler. For a cylinder of radius ``R`` producing
``q'''`` per volume, with flesh of conductivity ``k`` (Bird et al. 2002):

```math
T_s = T_c - \frac{q''' R^2}{4 k}
```

[`surface_and_lung_temperature`](@ref) holds this for each family of shape, and returns a lung temperature,
midway through the flesh, at which air is exhaled. The difference between core and surface is a fraction of a
degree in a small ectotherm, and several degrees in a large animal with a high metabolic rate.

## The residual

Given a core temperature, [`heat_balance`](@ref) computes the metabolic rate and respiration, then the surface
temperature, then every flow at the surface, and returns the difference between the two sides of the budget:

```@example balance
using HeatExchange, BiophysicalGeometry, Unitful

shape = DesertIguana(40.0u"g", 1000.0u"kg/m^3")
lizard = Organism(Body(shape, Naked()), example_ectotherm_heat_exchange_traits(; shape_pars = shape))
environment = (;
    environment_pars = example_environment_pars(),
    environment_vars = example_environment_vars(; air_temperature = u"K"(20.0u"°C"), wind_speed = 1.0u"m/s",
                                                  global_radiation = 800.0u"W/m^2", zenith_angle = 30.0u"°"),
)

u"W"(heat_balance(u"K"(20.0u"°C"), lizard, environment).heat_balance)
```

This is the *residual*. A lizard at 20 °C here gains that much more heat than it loses, so it would warm: the
residual is the rate at which it would store heat. At a higher temperature it loses more, and the residual
falls:

```@example balance
temperatures = 10.0:1.0:60.0   # °C
residuals = [ustrip(u"W", heat_balance(u"K"(T * u"°C"), lizard, environment).heat_balance) for T in temperatures]
solution = solve_temperature(lizard, environment)

fig, ax = figure_axis("Core temperature (°C)", "Residual of the heat budget (W)")
lines!(ax, temperatures, residuals; linewidth = 2)
hlines!(ax, [0.0]; color = :black, linewidth = 1)
scatter!(ax, [ustrip(u"°C", solution.core_temperature)], [0.0]; color = :black, markersize = 12)
fig
```

The steady-state body temperature is where the curve crosses zero. To *solve* the heat budget is to find that
crossing.

## Finding the root

The residual falls steadily as temperature rises, so there is one crossing. [`solve_temperature`](@ref) finds it
with Brent's method (Brent 2002), in [`zbrent`](@ref), between 270 K and 370 K by default, to 0.001 K:

```@example balance
solution = solve_temperature(lizard, environment)
u"°C"(solution.core_temperature), u"W"(solution.heat_balance)
```

The solution is the output of [`heat_balance`](@ref) at the root, with every term of the budget:

```@example balance
flow_table(solution.energy_balance) # hide
```

If no root lies in the bracket, the heat budget at air temperature is returned. Give another bracket with the
keyword `temperature_bracket`.

## With insulation

Fur and feathers add unknowns. The heat the body generates, less what leaves in the breath, is conducted to the
skin. Less what evaporates there, it passes through the fur. It then leaves the outer surface of the fur
(Kearney et al. 2021):

```math
Q_{gen,net} - Q_{evap} = Q_{fur} = Q_{env}
```

where

```math
Q_{gen,net} = Q_{gen} - Q_{resp}
\qquad \text{and} \qquad
Q_{env} = Q_{rad} + Q_{conv} + Q_{cond} + Q_{evap,fur} - Q_{sol}
```

```@example balance
layer_budget_diagram() # hide
```

There are now three temperatures:

| Symbol | In the code | Acted on by |
|:--|:--|:--|
| ``T_c``, core | `core_temperature` | conduction through flesh and fat, with ``T_s`` |
| ``T_s``, skin under the fur | `skin_temperature` | evaporation from the skin |
| ``T_{fa}``, fur–air interface | `insulation_temperature` | convection and radiation |

Heat moves through fur by conduction along the fibres and the air between them, and by radiation from fibre to
fibre. Both depend on the temperatures in the fur, so the conductivity of the fur is part of the solution
(Conley and Porter 1986), see [Insulation](insulation.md).

### Three residuals

[`solve_part_heat_balance`](@ref) is the insulated counterpart of [`heat_balance`](@ref). It takes ``T_c``,
``T_s``, ``T_{fa}`` and ``Q_{gen}`` and returns three residuals:

| Residual | Is zero when |
|:--|:--|
| `residual_energy_balance` | ``Q_{gen} + Q_{sol} = Q_{resp} + Q_{evap} + Q_{rad} + Q_{conv} + Q_{cond} + Q_{evap,fur}``: the whole budget closes |
| `residual_internal_conduction` | ``Q_{gen} - Q_{resp}`` equals the heat conducted from core to skin |
| `residual_skin_temperature` | ``T_s`` equals the skin temperature implied by the heat passing through the fur |

The third is continuity of heat flow at the skin: the skin temperatures implied from the core side and from the
fur side are averaged, and the residual is the distance from that mean, see
[Layers as a radial graph](radial_layers.md#Where-the-chain-is-not-a-chain).

Subtracting the second from the first removes ``Q_{gen}`` and ``Q_{resp}`` and leaves a balance of the surface
alone. This `surface_balance` and `residual_skin_temperature` are two equations in ``T_s`` and ``T_{fa}``, for
a given ``T_c``.

Here they are for a furred cylinder of 1 kg with its core at 37 °C in air at 20 °C. Each residual is zero along
a line, and the solution is where the lines cross:

```@example balance
using HeatExchange: EnvironmentTemperatures, ViewFactors, AtmosphericConditions

fibres = example_insulation_pars().dorsal
body = Body(Cylinder(1.0u"kg", 1000.0u"kg/m^3", 3.0),
            CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density), FatLayer(0.0, 901.0u"kg/m^3")))
vars = example_environment_vars(; air_temperature = u"K"(20.0u"°C"))
pars = example_environment_pars()
core_temperature = u"K"(37.0u"°C")

part = (;
    body,
    insulation_pars = example_insulation_pars(),
    traits = (; core_temperature, flesh_conductivity = 0.9u"W/m/K", fat_conductivity = 0.23u"W/m/K", ϵ_body = 0.99,
                skin_wetness = 0.005, insulation_wetness = 0.0, bare_skin_fraction = 0.0, eye_fraction = 0.0),
    environment_vars = (;
        temperature = EnvironmentTemperatures(vars),
        view_factors = ViewFactors(0.5, 0.5, 0.0, 0.0),
        atmos = AtmosphericConditions(vars),
        fluid = pars.fluid, solar_flow = 0.0u"W", gas_fractions = pars.gas_fractions,
        convection_enhancement = pars.convection_enhancement,
    ),
    conduction_fraction = 0.0, conductance_coefficient = 0.0u"W/K", ventral_fraction = 0.5,
    longwave_depth_fraction = 1.0, covered_area = 0.0u"m^2",
    characteristic_dim = characteristic_dimension(VolumeCubeRoot(), body),
)

part_residuals(skin, surface) = part_surface_residuals(part, core_temperature, u"K"(skin * u"°C"), u"K"(surface * u"°C"), 5.0u"W";
    k_flesh = 0.9u"W/m/K", pant = 1.0, skin_wetness = 0.005, resp_pars = example_respiration_pars())

skin = 30.0:0.1:37.0      # °C
surface = 20.0:0.1:34.0   # °C
surface_balance = [ustrip(u"W", part_residuals(s, f).surface_balance) for s in skin, f in surface]
skin_residual = [ustrip(u"K", part_residuals(s, f).residual_skin_temperature) for s in skin, f in surface]
solved = solve_part_surface(; part..., skin_temperature = u"K"(34.0u"°C"), insulation_temperature = u"K"(25.0u"°C"),
                              temperature_tolerance = 1e-3u"K")

fig, ax = figure_axis("Skin temperature (°C)", "Fur surface temperature (°C)")
contour!(ax, skin, surface, surface_balance; levels = [0.0], linewidth = 2, color = :firebrick)
contour!(ax, skin, surface, skin_residual; levels = [0.0], linewidth = 2, color = :steelblue)
scatter!(ax, [ustrip(u"°C", solved.skin_temperature)], [ustrip(u"°C", solved.insulation_temperature)];
         color = :black, markersize = 12)
axislegend(ax, [LineElement(color = :firebrick, linewidth = 2), LineElement(color = :steelblue, linewidth = 2)],
           ["surface balance = 0", "skin temperature residual = 0"]; position = :lt)
fig
```

[`solve_part_surface`](@ref) finds the point by a Newton iteration on the two residuals, and with it the heat
that must be conducted from the core to hold that state:

```@example balance
u"°C"(solved.skin_temperature), u"°C"(solved.insulation_temperature), u"W"(solved.net_metabolic)
```

### Closing the budget

What remains is the first equality: the heat generated, less respiration, must equal that ``Q_{gen,net}``. The
two solvers differ in which quantity is adjusted to make it so.

**Solving for metabolic rate.** ``T_c`` is given. The surface is solved once. Respiratory heat loss depends on
the air breathed, which depends on the metabolic rate, so ``Q_{gen}`` is the root of

```math
Q_{gen} - Q_{resp}(Q_{gen}) - Q_{gen,net} = 0
```

found with [`zbrent`](@ref). This residual is the `balance` returned by [`respiration`](@ref).

**Solving for temperature.** ``Q_{gen}`` is given by the metabolic rate equation at each temperature. The same
residual is then a function of ``T_c``, with the surface solved anew at each trial, and ``T_c`` is its root.

So an insulated heat budget is one search inside another. For bare skin the inner search disappears, as ``T_s``
follows directly from ``T_c``. See [Temperature or metabolic rate](solvers.md).

## The back and the belly

The fur and surroundings of the back differ from those of the belly. The solvers solve the surface twice, for
dorsal fur facing the sky and ventral fur facing the ground, and combine the two values of ``Q_{gen,net}`` by
the view of each, as NicheMapR does. A body can instead have real parts, each with its own surface solve, see
[Bodies of many parts](multipart.md) and [Back and belly: two halves](../tutorials/two_parts.md).

## Steady state and storage

All of the above is steady state. Away from it, the residual of [`heat_balance`](@ref) is the rate of heat
storage, and divided by the heat capacity of the body it is the rate of change of body temperature.

Transient heat budgets belong to this package and are to be added in that way: the residual tracked through
time, in an environment that may vary, and turned into body temperature by the heat capacity. That is the
equivalent of `onelump_var` of NicheMapR, and is not yet in the version documented here. What an organism does
while its temperature changes is for
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl).

## Residuals for other solvers

A root finder is one way to drive a residual to zero. Another is to give the residuals to an optimiser as
constraints, with the temperatures and metabolic rate among its variables.
[`solve_part_heat_balance`](@ref) and [`part_surface_residuals`](@ref) contain no search of their own for that
reason, see [Differentiability and the NLP interface](autodiff.md).
