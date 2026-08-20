using HeatExchange
using BiophysicalGeometry
using Unitful
using Test

# Self-consistency checks for the unified onelump/twolump transient physics: no closed
# form survives the unification (see src/transient.jl), so equilibria are checked against
# solve_temperature's independently-computed root instead — the same pattern
# test/onelump.jl already uses for the insulated branch against solve_metabolic_rate.

function rk4_step(f, u, t, dt)
    k1 = f(u, t)
    k2 = f(u + dt / 2 * k1, t + dt / 2)
    k3 = f(u + dt / 2 * k2, t + dt / 2)
    k4 = f(u + dt * k3, t + dt)
    return u + dt / 6 * (k1 + 2k2 + 2k3 + k4)
end

environment_pars = example_environment_pars()
environment_vars = example_environment_vars(;
    air_temperature=u"K"(20.0u"°C"), global_radiation=500.0u"W/m^2",
    zenith_angle=20.0u"°", wind_speed=1.0u"m/s",
)
e = (; environment_pars, environment_vars)

shapes = (
    Ellipsoid(0.5u"kg", 1000.0u"kg/m^3", 1.1, 1.1),
    Cylinder(0.5u"kg", 1000.0u"kg/m^3", 1.5),
)

core_temperature_init = u"K"(20.0u"°C")
dt = 5.0u"s"
nsteps = 7200
shell_thickness = 1.0e-3u"m"

@testset "onelump (Naked): RK4 converges to solve_temperature's equilibrium" for shape_pars in shapes
    body = Body(shape_pars, Naked())
    traits = example_heat_exchange_traits(; shape_pars)
    organism = Organism(body, traits)

    f(u, t) = onelump(u, t, organism, e).core_temperature_rate
    core_temperature = core_temperature_init
    for i in 1:nsteps
        core_temperature = rk4_step(f, core_temperature, (i - 1) * dt, dt)
    end

    equilibrium = solve_temperature(organism, e).core_temperature
    @test ustrip(u"K", core_temperature) ≈ ustrip(u"K", equilibrium) rtol = 1e-3
end

@testset "twolump: RK4-converged trajectory is physically consistent" for shape_pars in shapes
    body = Body(shape_pars, Naked())
    traits = example_heat_exchange_traits(; shape_pars)
    organism = Organism(body, traits)

    # (b) twolump (LinearizedSurface) vs. onelump, both RK4-converged: a loose
    # self-consistency check, not a tight one — twolump omits evaporation/respiration.
    f_one(u, t) = onelump(u, t, organism, e).core_temperature_rate
    core_one = core_temperature_init
    for i in 1:nsteps
        core_one = rk4_step(f_one, core_one, (i - 1) * dt, dt)
    end

    state = (; core_temperature=core_temperature_init, shell_temperature=core_temperature_init)
    for i in 1:nsteps
        out = twolump(state, (i - 1) * dt, organism, e; shell_thickness)
        state = (;
            core_temperature=state.core_temperature + out.core_temperature_rate * dt,
            shell_temperature=state.shell_temperature + out.shell_temperature_rate * dt,
        )
    end
    @test ustrip(u"K", state.core_temperature) ≈ ustrip(u"K", core_one) atol = 5.0

    # (c) final_core_temperature at the converged state should match it, rates ≈ 0.
    out_converged = twolump(state, nsteps * dt, organism, e; shell_thickness)
    @test ustrip(u"K", out_converged.final_core_temperature) ≈ ustrip(u"K", state.core_temperature) rtol = 1e-2
    @test ustrip(u"K/s", out_converged.core_temperature_rate) ≈ 0.0 atol = 1e-4
    @test ustrip(u"K/s", out_converged.shell_temperature_rate) ≈ 0.0 atol = 1e-4

    # (d) LinearizedSurface vs RootFindSurface should agree near steady state: both give
    # rates ≈ 0 and matching surface_temperature.
    out_rootfind = twolump(state, nsteps * dt, organism, e; shell_thickness, surface_solve=RootFindSurface())
    @test ustrip(u"K/s", out_rootfind.core_temperature_rate) ≈ 0.0 atol = 1e-4
    @test ustrip(u"K/s", out_rootfind.shell_temperature_rate) ≈ 0.0 atol = 1e-4
    @test ustrip(u"K", out_rootfind.surface_temperature) ≈ ustrip(u"K", out_converged.surface_temperature) rtol = 1e-2
end
