# CommonSolve.jl interface for the endotherm heat-balance solve.
#
# `HeatBalanceProblem` is a self-contained heat-balance problem; `HeatBalanceSolver`
# is a reusable workspace carrying the current temperature state across successive
# solves. This exposes the standard SciML `init` / `solve!` / `reinit!` / `solve`
# surface so downstream code composes with the rest of the ecosystem and so
# `BiophysicalBehaviour`'s control loop can warm-start across effector iterations.
#
# The physics behind `solve!` is currently the validated single-body /
# dorsal-ventral `solve_metabolic_rate`. As the multi-part per-part primitive
# (`solve_part_heat_balance`) and compartment graph land, they slot in behind
# this same interface without changing the surface. The problem stores the
# organism + environment for now; the fuller per-part
# geometry/physiology/couplings/compartment_graph decomposition (plan §3.6) is
# added as multi-part physics arrives.

using CommonSolve: CommonSolve, init, solve!, solve

# `reinit!` is not part of CommonSolve (it is a SciMLBase generic); define our own.
"""
    reinit!(solver, new_problem; warm_start = true)

Reset the solver for a new problem, optionally reusing the current temperature
state (warm start).
"""
function reinit! end

"""
    HeatBalanceProblem(organism, environment)

A self-contained multi-part heat-balance problem: everything required to define
the physics (topology, physiology, environment). No solver state — initial
state flows in through `init` / `reinit!` (SciML convention).
"""
struct HeatBalanceProblem{O,E}
    organism::O
    environment::E
end

"""
    HeatBalanceSolver

Mutable workspace + current temperature state, reused across successive solves.
Per-solve output metadata lives on the returned solution, not here.
"""
mutable struct HeatBalanceSolver{P<:HeatBalanceProblem,S,SM}
    problem::P
    state::S            # current best estimate: (; skin_temperature, insulation_temperature)
    smoothing::SM
    solution::Any       # last ThermoregulationOutput, or nothing before the first solve!
end

"""
    default_initial_state(problem) -> (; skin_temperature, insulation_temperature)

Standard initial guesses: skin 3 K below the setpoint core temperature,
insulation surface at air temperature.
"""
function default_initial_state(problem::HeatBalanceProblem)
    o = problem.organism
    core_temperature = metabolism_pars(o).core_temperature
    return (;
        skin_temperature       = core_temperature - 3u"K",
        insulation_temperature = problem.environment.environment_vars.air_temperature,
    )
end

function CommonSolve.init(problem::HeatBalanceProblem;
                          initial_state = default_initial_state(problem),
                          smoothing::SmoothingStrategy = HardBound(), kwargs...)
    return HeatBalanceSolver(problem, initial_state, smoothing, nothing)
end

function CommonSolve.solve!(solver::HeatBalanceSolver)
    problem = solver.problem
    state = solver.state
    output = solve_metabolic_rate(
        problem.organism, problem.environment,
        state.skin_temperature, state.insulation_temperature;
        smoothing = solver.smoothing,
    )
    # Thread the converged surface temperatures forward for the next solve.
    solver.state = (;
        skin_temperature       = output.thermoregulation.skin_temperature,
        insulation_temperature = output.thermoregulation.insulation_temperature,
    )
    solver.solution = output
    return output
end

function reinit!(solver::HeatBalanceSolver, new_problem::HeatBalanceProblem;
                warm_start::Bool = true)
    solver.problem = new_problem
    warm_start || (solver.state = default_initial_state(new_problem))
    return solver
end
