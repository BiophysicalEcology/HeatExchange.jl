# Temperature or metabolic rate

There are two solvers, for the two quantities that can be unknown, and each works with bare skin or insulation:

| | Bare skin | Insulated |
|:--|:--|:--|
| [`solve_temperature`](@ref): metabolic rate known, core temperature unknown | a lizard, a frog, an insect, a leaf | a torpid mammal, a nestling, a carcass, a bumblebee |
| [`solve_metabolic_rate`](@ref): core temperature known, metabolic rate unknown | a naked mammal, a thermogenic flower | a mammal, a bird |

[Solving a heat balance](heat_balance.md) explains what each does. This page is about calling them and reading
their output.

```@setup solvers
using Main.FigureHelpers
using Main.GeometryFigures
using CairoMakie
```

## The organisms

Two bodies of the same shape and mass, one bare and one furred:

```@example solvers
using HeatExchange, BiophysicalGeometry, Unitful

shape = Ellipsoid(2.0u"kg", 1000.0u"kg/m^3", 2.0, 2.0)
fibres = FibreProperties(; diameter = 30.0u"μm", length = 25.0u"mm", density = 3000.0u"cm^-2", depth = 15.0u"mm",
                           reflectance = 0.2, conductivity = 0.209u"W/m/K")
fat = FatLayer(0.0, 901.0u"kg/m^3")

bare_body = Body(shape, Naked())
furred_body = Body(shape, CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density), fat))

metabolism_pars = MetabolismParameters(; core_temperature = u"K"(38.0u"°C"),
                                         metabolic_heat_flow = metabolic_rate(Kleiber(), 2.0u"kg"), model = Kleiber())
bare_fibres = FibreProperties(; diameter = 30.0u"μm", length = 0.0u"mm", density = 0.0u"cm^-2", depth = 0.0u"mm",
                                reflectance = 0.2, conductivity = 0.209u"W/m/K")
bare_traits = example_heat_exchange_traits(; shape_pars = shape, metabolism_pars,
    insulation_pars = InsulationParameters(; dorsal = bare_fibres, ventral = bare_fibres,
                                             depth_compressed = 0.0u"mm", longwave_depth_fraction = 1.0),
    evaporation_pars = example_evaporation_pars(; bare_skin_fraction = 1.0))
furred_traits = example_heat_exchange_traits(; shape_pars = shape, metabolism_pars,
    insulation_pars = InsulationParameters(; dorsal = fibres, ventral = fibres, depth_compressed = fibres.depth,
                                             longwave_depth_fraction = 1.0))

bare = Organism(bare_body, bare_traits)
furred = Organism(furred_body, furred_traits)
environment = (; environment_pars = example_environment_pars(),
                 environment_vars = example_environment_vars(; air_temperature = u"K"(10.0u"°C"), wind_speed = 1.0u"m/s"))
nothing # hide
```

```@example solvers
shape_gallery("bare" => bare_body, "furred" => furred_body) # hide
```

The metabolic rate in [`MetabolismParameters`](@ref) has two roles. `model` is the equation used when
temperature is solved for. `metabolic_heat_flow` is the minimum when metabolic rate is solved for, here the
basal rate.

## Solving for temperature

```@example solvers
out = solve_temperature(bare, environment)
u"°C"(out.core_temperature)
```

For bare skin the output is that of [`heat_balance`](@ref) at the solution, a NamedTuple:

| Field | Content |
|:--|:--|
| `core_temperature`, `surface_temperature`, `lung_temperature` | temperatures |
| `heat_balance` | the residual, near zero |
| `energy_balance` | each term of the heat budget |
| `mass_balance` | oxygen consumed, and water evaporated from the lungs, skin and eyes |
| `solar_out`, `longwave_gain_out`, `longwave_loss_out`, `convection_out`, `evaporation_out`, `respiration_out` | the full output of each heat-flow function |

For an insulated body the same call returns a [`ThermoregulationOutput`](@ref), described below:

```@example solvers
out = solve_temperature(furred, environment)
u"°C"(out.thermoregulation.core_temperature)
```

With the same metabolic rate, the furred animal settles well above the bare one. The keyword
`temperature_bracket` sets the range searched, 270 K to 370 K by default.

## Solving for metabolic rate

[`solve_metabolic_rate`](@ref) holds the core at `metabolism_pars.core_temperature`. It needs first guesses of
the skin temperature and of the outer surface temperature, which for a bare body is also the skin:

```@example solvers
core_temperature = metabolism_pars.core_temperature
cold_furred = solve_metabolic_rate(furred, environment, core_temperature - 3u"K", environment.environment_vars.air_temperature)
cold_bare = solve_metabolic_rate(bare, environment, core_temperature - 3u"K", core_temperature - 3u"K")
cold_furred.energy_flows.metabolic_heat_flow, cold_bare.energy_flows.metabolic_heat_flow
```

Against a basal rate of

```@example solvers
metabolism_pars.metabolic_heat_flow
```

the furred animal needs about a third more than basal at 10 °C, and the bare one more than five times as much.

### A result below the minimum

The metabolic rate returned is what the heat budget requires, not what an animal can do. In a warm environment
it can be less than basal, and in a hot one negative:

```@example solvers
warm = (; environment_pars = example_environment_pars(),
          environment_vars = example_environment_vars(; air_temperature = u"K"(35.0u"°C"), wind_speed = 1.0u"m/s"))
solve_metabolic_rate(furred, warm, core_temperature - 3u"K", u"K"(35.0u"°C")).energy_flows.metabolic_heat_flow
```

A value below the minimum means the animal makes more heat than it can lose in that state. It must change
something: flatten its fur, send blood to its skin, pant, sweat, or let its core warm (Kearney et al. 2021).
Each is a change to the organism followed by another solve. The sequence is the work of
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), see
[Endotherm thermoregulation by rules](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/endotherm_rules).
This package reports the requirement.

## The output

A [`ThermoregulationOutput`](@ref) has four groups, the tables `treg`, `morph`, `enbal` and `masbal` of the
NicheMapR endotherm model:

::: tabs

== thermoregulation

Temperatures and the state of the traits that thermoregulation changes, with each side in `dorsal` and
`ventral`.

```@example solvers
flow_table(getfield(cold_furred.thermoregulation, :data); header = ["Field", "Value"]) # hide
```

== morphology

Areas, volumes and lengths of the body.

```@example solvers
flow_table(getfield(cold_furred.morphology, :data); header = ["Field", "Value"]) # hide
```

== energy_flows

The terms of the heat budget, in W. `evaporation_heat_flow` includes respiration. `heat_balance` is the residual
of the whole budget and `balance` that of the respiration balance. `success` and `ntry` report the surface
solve. Each side is in `dorsal` and `ventral`.

```@example solvers
flow_table(getfield(cold_furred.energy_flows, :data); header = ["Field", "Value"]) # hide
```

== mass_flows

Air breathed, oxygen consumed and water evaporated, with the molar flows of each gas in `molar_fluxes_in` and
`molar_fluxes_out`, see [Flows of mass](gradients.md#Flows-of-mass).

```@example solvers
flow_table(getfield(cold_furred.mass_flows, :data); header = ["Field", "Value"]) # hide
```

:::

Printing a result at the REPL shows all four groups in full.

## Which path is taken

Whether a body is treated as bare or insulated follows from its layers, through [`evaluation_strategy`](@ref):

```@example solvers
evaluation_strategy(bare), evaluation_strategy(furred)
```

[`SingleBody`](@ref) solves one surface. [`MultiSided`](@ref) solves a dorsal and a ventral side and combines
them. A fibrous layer of zero depth is solved as bare skin: the test is in [`insulation_properties`](@ref), on
the fibres in the traits.

Two limits of the current version:

- [`heat_balance`](@ref) for a whole organism is defined for bare bodies only. For an insulated body the
  residuals come from [`solve_part_heat_balance`](@ref).
- The layers of a body must be `Naked()`, a `FibrousLayer`, or a `CompositeInsulation` of fibres and fat. Fat
  with no fur is a `CompositeInsulation` whose fibrous layer has zero depth.

## Many solves in a row

A simulation solves the same animal for many hours or states in turn, and the skin and fur temperatures of one
solve are good first guesses for the next. The [CommonSolve.jl](https://github.com/SciML/CommonSolve.jl)
interface keeps them:

```@example solvers
problem = HeatBalanceProblem(furred, environment)
solver = init(problem)
first_solve = solve!(solver)
solver.state
```

```@example solvers
air_temperatures = -20.0:5.0:30.0   # °C
metabolic_rates = map(air_temperatures) do air_temperature
    vars = example_environment_vars(; air_temperature = u"K"(air_temperature * u"°C"), wind_speed = 1.0u"m/s")
    reinit!(solver, HeatBalanceProblem(furred, (; environment_pars = example_environment_pars(), environment_vars = vars)))
    solve!(solver).energy_flows.metabolic_heat_flow
end

fig, ax = figure_axis("Air temperature (°C)", "Metabolic rate required (W)")
scatterlines!(ax, air_temperatures, ustrip.(u"W", metabolic_rates); linewidth = 2)
hlines!(ax, [ustrip(u"W", metabolism_pars.metabolic_heat_flow)]; color = :black, linestyle = :dash)
fig
```

The dashed line is the basal rate. `solve(problem)` does one solve without keeping a solver, and
`reinit!(solver, problem; warm_start = false)` discards the stored temperatures.
