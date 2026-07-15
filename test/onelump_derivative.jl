using HeatExchange
using BiophysicalGeometry
using Unitful
using Test

# Self-consistency checks for the transient lumped-capacitance models: no R reference
# data exists for these yet (see test/R/onelump_test.R, a manual step), so these check
# the derivative-form onelump/twolump against the closed-form onelump under a constant
# environment, which is the cross-check the plan specifies in place of an R reference.

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
internal_conduction = example_conduction_pars_internal()
kw = (;
    internal_conduction, posture=Intermediate(),
    body_absorptivity=0.85, emissivity=0.95, sky_view_factor=0.4, ground_view_factor=0.4,
    metabolic_heat_volumetric=0.0u"W/m^3",
)

@testset "onelump: derivative matches closed form" for body in (
    Body(Ellipsoid(0.5u"kg", 1000.0u"kg/m^3", 1.1, 1.1), Naked()),
    Body(Cylinder(0.5u"kg", 1000.0u"kg/m^3", 1.5), Naked()),
)
    core_temperature_init = u"K"(20.0u"°C")
    closed = onelump((1:60:36000)u"s", core_temperature_init, body, environment_pars, environment_vars; kw...)

    initial_rate = onelump(core_temperature_init, 0.0u"s", body, environment_pars, environment_vars; kw...)
    @test ustrip(u"K/s", initial_rate) ≈ ustrip(u"K/s", closed.initial_rate) rtol = 1e-8

    dt = 5.0u"s"
    f(u, t) = onelump(u, t, body, environment_pars, environment_vars; kw...)
    core_temperature = core_temperature_init
    core_temperature_trace = [core_temperature]
    for i in 1:7200
        core_temperature = rk4_step(f, core_temperature, (i - 1) * dt, dt)
        push!(core_temperature_trace, core_temperature)
    end
    @test ustrip(u"K", core_temperature) ≈ ustrip(u"K", closed.final_core_temperature) rtol = 1e-3

    # Uncomment to compare the RK4 trajectory against the closed form visually:
    # using Plots
    # plot((0:7200) * dt, core_temperature_trace; label="RK4 (derivative form)")
    # plot!((1:60:36000)u"s", closed.core_temperature; label="closed form")
end

@testset "twolump: steady state matches onelump" for body in (
    Body(Ellipsoid(0.5u"kg", 1000.0u"kg/m^3", 1.1, 1.1), Naked()),
    Body(Cylinder(0.5u"kg", 1000.0u"kg/m^3", 1.5), Naked()),
)
    core_temperature_init = u"K"(20.0u"°C")
    one = onelump((1:60:36000)u"s", core_temperature_init, body, environment_pars, environment_vars; kw...)

    two_kw = (;
        internal_conduction, shell_thickness=1.0e-3u"m", posture=Intermediate(),
        body_absorptivity=0.85, emissivity=0.95, sky_view_factor=0.4, ground_view_factor=0.4,
        metabolic_heat_volumetric=0.0u"W/m^3",
    )
    two = twolump(
        (; core_temperature=core_temperature_init, shell_temperature=core_temperature_init),
        0.0u"s", body, environment_pars, environment_vars; two_kw...,
    )
    @test ustrip(u"K", two.final_core_temperature) ≈ ustrip(u"K", one.final_core_temperature) rtol = 1e-2
end
