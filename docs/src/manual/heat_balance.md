# Solving a heat balance

This page sets out the heat budget that the package solves, what it means to solve it, and how the solution is
found, first for an organism with bare skin and then for one with fur or feathers.

```@setup balance
using Main.FigureHelpers
using CairoMakie
```

## The budget

At steady state the heat that an organism gains equals the heat that it loses (Porter and Gates 1969):

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

By convention the terms on the right are losses, and a negative value is a gain: convection from warmer air,
conduction from a warmer rock, or condensation. Each is described on its own page of this manual.

The first two terms on the left are set by the environment. Every other term depends on the temperature of the
organism, and that is what makes the budget something to be solved and not simply added up.

Almost every term has one form, a flow down a gradient of a potential against a resistance, see
[Gradients, resistances and flows](gradients.md).

## From the core to the surface

Heat is exchanged with the environment at the surface, but the temperature that matters to the organism is that
of its core. Metabolic heat is produced throughout the flesh and conducted outwards, so the surface is cooler than
the core. For a cylinder of radius ``R`` producing heat at ``q'''`` per volume, with flesh of conductivity ``k``
(Bird et al. 2002):

```math
T_s = T_c - \frac{q''' R^2}{4 k}
```

[`surface_and_lung_temperature`](@ref) holds this for each family of shape. It also returns a lung temperature,
midway through the flesh, at which air is exhaled. For a small ectotherm the difference between core and surface
is a small fraction of a degree. For a large animal with a high metabolic rate it is several degrees.

## The residual

[`heat_balance`](@ref) puts these together. Given a core temperature, it computes the metabolic rate and the
respiration at that temperature, then the surface temperature, then every flow at the surface, and returns the
difference between the two sides of the budget:

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

This number is the *residual* of the heat budget. A lizard at 20 °C in this environment gains that much more heat
than it loses, so it is not at steady state: it would warm up. The residual is the rate at which it would store
heat. At a higher temperature it loses more, and the residual falls:

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

The residual falls steadily as the temperature rises, so there is one crossing, and it can be found by a
bracketing root finder. [`solve_temperature`](@ref) uses Brent's method (Brent 2002), in [`zbrent`](@ref), between
270 K and 370 K by default, to a tolerance of 0.001 K:

```@example balance
solution = solve_temperature(lizard, environment)
u"°C"(solution.core_temperature), u"W"(solution.heat_balance)
```

The solution is the output of [`heat_balance`](@ref) at the root, with every term of the budget:

```@example balance
flow_table(solution.energy_balance) # hide
```

If no root lies in the bracket, the heat budget at air temperature is returned. A different bracket is given with
the keyword `temperature_bracket`.

## With insulation

Fur and feathers put a layer between the skin and the environment, and that adds unknowns. The heat that the body
generates, less what leaves in the breath, must be conducted to the skin. Less what evaporates from the skin, it
must pass through the fur. And it must then leave the outer surface of the fur (Kearney et al. 2021):

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

There are now three temperatures: the core temperature ``T_c``, the skin temperature ``T_s`` under the fur, and
the temperature of the outer surface of the fur, the fur–air interface, ``T_{fa}``. In the code these are
`core_temperature`, `skin_temperature` and `insulation_temperature`. Convection and radiation act on
``T_{fa}``, evaporation of sweat on ``T_s``, and conduction through the flesh and fat on ``T_c - T_s``.

Heat moves through fur by conduction along the fibres and through the air between them, and by radiation from
fibre to fibre. Both depend on the temperatures in the fur, so the conductivity of the fur is itself part of the
solution (Conley and Porter 1986), see [Insulation](insulation.md).

### Three residuals

[`solve_part_heat_balance`](@ref) is the insulated counterpart of [`heat_balance`](@ref). It takes ``T_c``,
``T_s``, ``T_{fa}`` and ``Q_{gen}`` and returns three residuals, one for each equality that must hold:

| Residual | Is zero when |
|:--|:--|
| `residual_energy_balance` | ``Q_{gen} + Q_{sol} = Q_{resp} + Q_{evap} + Q_{rad} + Q_{conv} + Q_{cond} + Q_{evap,fur}``: the whole budget closes |
| `residual_internal_conduction` | ``Q_{gen} - Q_{resp}`` equals the heat conducted from core to skin through flesh and fat |
| `residual_skin_temperature` | ``T_s`` equals the skin temperature implied by the heat passing through the fur |

The third is continuity of heat flow at the skin. The skin temperature implied from the core side and that
implied from the fur side are averaged, and the residual is the distance of the skin temperature from that
mean. At a solution the two estimates agree to rounding error, see
[Layers as a radial graph](radial_layers.md#Where-the-chain-is-not-a-chain).

Subtracting the second from the first removes ``Q_{gen}`` and ``Q_{resp}``, and leaves a balance of the surface
alone: the heat conducted to the skin, plus sunlight, equals the losses from the skin and the fur. This
`surface_balance` and `residual_skin_temperature` are two equations in the two unknowns ``T_s`` and ``T_{fa}``,
for a given ``T_c``.

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

The point is found by a Newton iteration on the two residuals, in [`solve_part_surface`](@ref), and with it the
heat that must be conducted from the core to hold that state:

```@example balance
u"°C"(solved.skin_temperature), u"°C"(solved.insulation_temperature), u"W"(solved.net_metabolic)
```

### Closing the budget

What remains is the first equality: the heat generated, less respiration, must equal that ``Q_{gen,net}``. Which
quantity is adjusted to make it so is the difference between the two solvers.

**Solving for metabolic rate.** ``T_c`` is given. The surface is solved once for ``T_s``, ``T_{fa}`` and
``Q_{gen,net}``. Respiratory heat loss depends on the air breathed, which depends on the metabolic rate, so
``Q_{gen}`` is the root of

```math
Q_{gen} - Q_{resp}(Q_{gen}) - Q_{gen,net} = 0
```

found with [`zbrent`](@ref). This residual is the `balance` returned by [`respiration`](@ref).

**Solving for temperature.** ``Q_{gen}`` is given, by the metabolic rate equation at each temperature. The same
residual is then a function of ``T_c``, with the surface solved anew at each trial value, and ``T_c`` is its
root.

So an insulated heat budget is one search inside another: a pair of surface temperatures for each trial of the
outer unknown. For a body with bare skin the inner search disappears, as ``T_s`` follows directly from ``T_c``.

## The back and the belly

The fur and the surroundings of the back of an animal differ from those of its belly. The solvers handle this
by solving the surface twice, once for the dorsal fur facing the sky and once for the ventral fur facing the
ground, and combining the two values of ``Q_{gen,net}`` in proportion to the view of each, as NicheMapR does.
A body can instead be given real parts, each with its own surface solve, see
[Bodies of many parts](multipart.md).

## Steady state and storage

All of the above is for steady state. The residual of [`heat_balance`](@ref) at a temperature that is not the
solution is the rate of heat storage at that temperature, so the same function gives the rate of change of body
temperature when divided by the heat capacity of the body. Transient heat budgets built in this way are in
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl).

## Residuals for other solvers

A root finder is one way to drive a residual to zero. Another is to hand the residuals to an optimiser as
constraints, with the temperatures and the metabolic rate among its variables. The functions that return
residuals, [`solve_part_heat_balance`](@ref) and [`part_surface_residuals`](@ref), contain no search of their own
for that reason, see [Differentiability and the NLP interface](autodiff.md).
