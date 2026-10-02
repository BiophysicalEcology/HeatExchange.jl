# Differentiability and the NLP interface

[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl) can find the
thermoregulatory response of an animal by optimisation. Vasodilation, panting, sweating and the temperatures of
the body are the variables of a nonlinear program, the heat budget is its constraints, and
[IPOPT](https://github.com/jump-dev/Ipopt.jl) solves it. IPOPT needs first and second derivatives of the
constraints, computed by automatic differentiation of this package with
[Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl).

This page describes the features of HeatExchange.jl that make that possible. None needs to be used, or
understood, to solve a heat budget with [`solve_temperature`](@ref) or [`solve_metabolic_rate`](@ref). For the
other side of the boundary, see
[How the optimisation is built](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/nlp)
in the documentation of BiophysicalBehaviour.jl.

```@setup autodiff
using Main.FigureHelpers
using CairoMakie
using HeatExchange
```

## Residuals in place of solutions

[`solve_metabolic_rate`](@ref) contains two searches, one inside the other: for the skin and fur temperatures
of each side, and for the metabolic rate that closes the respiration balance, see
[Solving a heat balance](heat_balance.md#Closing-the-budget). An optimiser cannot use it as a constraint. If
the optimiser proposes a skin temperature, a function that searches for its own ignores the proposal. And a
derivative taken through an iteration is costly and fragile.

So the physics is also written with no search in it. [`solve_part_heat_balance`](@ref) takes the core, skin and
fur surface temperatures and the metabolic rate as arguments, and returns how far each balance is from
satisfied:

```julia
balance = solve_part_heat_balance(core_temperature, skin_temperature, insulation_temperature, metabolic_heat_flow;
    body, geometry, insulation_pars, insulation, geometry_vars, environment_vars, traits, resp_pars)

balance.residual_energy_balance        # W, zero when heat gained equals heat lost
balance.residual_internal_conduction   # W, zero when heat generated equals heat conducted to the skin
balance.residual_skin_temperature      # K, zero when the skin temperature is consistent with both
```

The search is done by the caller. The rule-based solver here drives the residuals to zero with a Newton
iteration. The optimiser of BiophysicalBehaviour.jl holds them at zero as equality constraints. Both call the
same function, so the two cannot drift apart:

```@example autodiff
solver_paths_diagram() # hide
```

[`part_surface_residuals`](@ref) is the form for a body of several parts. It takes the same fixed description
of a part as the iterative [`solve_part_surface`](@ref), and both build their geometry and insulation with one
shared function, so they agree residual for residual at a solution.

Variables that an optimiser changes are positional arguments, and everything fixed during a solve is a keyword.
Derivatives are taken with respect to what comes before the semicolon.

## Smoothing

A heat budget has corners. Free convection depends on the *magnitude* of the temperature difference. Respiratory
heat loss uses the *larger* of the metabolic rate and its minimum. The conductivity of fur is *clamped* between
that of air and of keratin. Wet fur evaporates *only if* there is fur. In code these are `abs`, `max`, `clamp`
and `if`.

For a forward calculation a corner is exact and harmless. For a derivative it is not:

- At the corner of `abs(x)` the slope jumps from −1 to 1, and an optimiser that steps across it gets a gradient
  that points the wrong way.
- A number raised to a fractional power has an infinite slope at zero, and infinity times a zero elsewhere in
  the chain rule is `NaN`, which spreads to every derivative. IPOPT starts with the fur surface at air
  temperature, which is exactly such a point.

Each corner is therefore written with a function that takes a [`SmoothingStrategy`](@ref):

| Function | [`HardBound`](@ref) | [`SmoothBound`](@ref) |
|:--|:--|:--|
| [`safe_abs`](@ref) | `abs(x)` | ``\sqrt{x^2 + (\varepsilon s)^2}`` |
| [`safe_relu`](@ref) | `max(x, 0)` | ``(x + \mathrm{safe\_abs}(x)) / 2`` |
| [`safe_step`](@ref) | `x > 0 ? 1 : 0` | ``(1 + x / \mathrm{safe\_abs}(x)) / 2`` |
| [`safe_max`](@ref), [`safe_min`](@ref) | `max(a, b)`, `min(a, b)` | ``(a + b \pm \mathrm{safe\_abs}(a - b)) / 2`` |
| [`safe_clamp`](@ref) | `clamp(x, lo, hi)` | `safe_max(lo, safe_min(hi, x))` |

[`HardBound`](@ref) is the default, and gives the exact corner. Every result in the rest of this documentation,
and every comparison with NicheMapR, uses it. `SmoothBound(ε)` rounds each corner over a width set by `ε`:

```@example autodiff
smoothing_figure() # hide
```

The width of the rounding is `ε` times a scale, `s` above, given where the function is called, because only
there is it known what a small value is. A skin wetness lies between 0 and 0.05, a temperature difference
between 0 and 50 K, a fur depth is a few millimetres:

```julia
temperature_difference = safe_abs(smoothing, ΔT; scale = 1.0u"K")
mask = safe_step(smoothing, insulation_wetness; scale = 0.01)
```

So `ε` is one dimensionless number for the whole model, 10⁻³ by default, and each corner is rounded in
proportion to its own quantity. Away from a corner the two strategies agree to many digits. The strategy is
passed down through every function as the keyword `smoothing`:

```julia
solve_metabolic_rate(organism, environment, skin_temperature, insulation_temperature; smoothing = SmoothBound(1e-3))
```

Two corners are not smoothed. A division whose numerator and denominator both vanish is kept off the path by an
ordinary `if`, since blending the branches would put the division back. And the temperature difference in free
convection has a floor of 10⁻⁶ K under either strategy, which removes the infinite slope without changing any
result.

## Types that do not change

Automatic differentiation, and fast code in general, needs every function to return the same type whichever
branch is taken.

- **Branches return the same type.** A function returning `0.0u"W"` in one branch and a computed power in the
  other returns one of two types. Such results are multiplied by a mask from [`safe_step`](@ref) instead.
- **Choices fixed for a solve are types.** The side of the body ([`Dorsal`](@ref), [`Ventral`](@ref)), the
  fluid ([`Air`](@ref), [`Water`](@ref)), the family of the shape and the kind of [`HeatCoupling`](@ref) are
  types, so the compiler picks the method once and the branch not taken is not compiled.
- **Units are made canonical where they meet.** Unitful.jl treats `mm/m` and a plain number as different types.
  Quantities that could arrive in different units are converted to one before they are combined.
- **No part costs nothing.** A part with no neighbours gets an empty tuple of neighbours, and the method for an
  empty tuple returns its input unchanged. A body of one part compiles to exactly the code it had before bodies
  of many parts existed.

## Geometry outside the loop

The areas of a body do not change while its temperatures are solved. They are computed once, by
BiophysicalGeometry.jl, and passed to [`solve_part_heat_balance`](@ref) as the `geometry` keyword. The function
that is differentiated then contains no geometry, which keeps the derivative code small.

Changes of posture and fur depth do change the geometry. They are made between solves, by rebuilding the body.

## A fixed number of compartments

For a body of several parts, the core temperatures of the compartments solve a small linear system, see
[Bodies of many parts](multipart.md#The-core-temperatures). The number of compartments is a type parameter of
[`CompartmentGraph`](@ref), found once when the body is described. The matrix is then a
[StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) `SMatrix` of known size, solved without heap
memory.

## Units

Every quantity carries units, and a derivative of a temperature with respect to a conductivity is not a
temperature. The two are kept apart at one boundary. The optimiser works on a plain vector of `Float64`. The
calling package attaches units on the way in and removes them from the residuals on the way out, and between
those steps the physics runs with units as it does anywhere else. The one place this package removes units
itself is the linear solve for compartment core temperatures, in W/K and W. See
[Units, dimensions and functional traits](units_traits.md#Where-units-are-stripped-in-this-package).

## The interface

This package owns the physics. BiophysicalBehaviour.jl owns the optimisation: the bounds on each variable, the
objective, the order of the variables and the solver. The interface between them is small:

```julia
abstract type NLPStrategy end
function nlp_pack end
```

BiophysicalBehaviour.jl defines a strategy, `MultipartNLP <: NLPStrategy`, and a method of [`nlp_pack`](@ref)
for it. The method does once, before a solve, everything that does not depend on the variables: it builds the
fixed description of each part. Each evaluation of the constraints then calls
[`part_surface_residuals`](@ref) for each part with the current variables, and adds one residual for the
respiration balance of the whole animal. No heat exchange calculation is repeated in the other package.

## Starting from the last solution

Thermoregulation by rules solves the heat budget many times in a row with small changes to the animal between
solves. Each solve starts best from the answer to the one before. The
[CommonSolve.jl](https://github.com/SciML/CommonSolve.jl) interface holds that state:

```julia
solver = init(HeatBalanceProblem(organism, environment))
output = solve!(solver)                     # the skin and fur temperatures found are kept in `solver`
reinit!(solver, HeatBalanceProblem(new_organism, environment))
output = solve!(solver)                     # starts from them
```

See [Temperature or metabolic rate](solvers.md#Many-solves-in-a-row).
