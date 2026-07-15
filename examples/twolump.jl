# Standalone two-lump (core + shell) transient body temperature under a full diurnal
# cycle, mimicking twolump.R's own example (which drives the model through a day of
# real microclimate data). Here the environment is a simple synthetic sinusoid instead
# — the time-varying-forcing/interpolation machinery lives in BiophysicalBehaviour.jl's
# EnvironmentForcing, not here. Demonstrates the two-lump model's core point: the shell
# tracks the environment more closely than the core, so the two diverge over the day.
using HeatExchange
using BiophysicalGeometry
using Unitful

body = Body(Ellipsoid(0.5u"kg", 1000.0u"kg/m^3", 1.1, 1.1), Naked())
environment_pars = example_environment_pars()
internal_conduction = example_conduction_pars_internal()

function environment_at(t)
    hour = ustrip(u"hr", t)
    air_temperature = u"K"((20.0 + 8.0 * sin(2π * (hour - 6) / 24))u"°C")
    zenith_angle = clamp(90.0 - 90.0 * sin(2π * (hour - 6) / 24), 0.0, 90.0)u"°"
    global_radiation = max(0.0, 800.0 * sin(2π * (hour - 6) / 24))u"W/m^2"
    example_environment_vars(; air_temperature, global_radiation, zenith_angle, wind_speed=1.0u"m/s")
end

function diurnal_trajectory(times, core_temperature_init, shell_temperature_init)
    kw = (;
        internal_conduction, shell_thickness=1.0e-3u"m", posture=Intermediate(),
        body_absorptivity=0.85, emissivity=0.95,
        sky_view_factor=0.4, ground_view_factor=0.4,
        metabolic_heat_volumetric=0.0u"W/m^3",
    )
    core_temperature, shell_temperature = core_temperature_init, shell_temperature_init
    core_trace, shell_trace = [core_temperature], [shell_temperature]
    for i in 2:length(times)
        dt = times[i] - times[i - 1]
        environment_vars = environment_at(times[i - 1])
        out = twolump((; core_temperature, shell_temperature), times[i - 1], body, environment_pars, environment_vars; kw...)
        core_temperature += out.core_temperature_rate * dt
        shell_temperature += out.shell_temperature_rate * dt
        push!(core_trace, core_temperature)
        push!(shell_trace, shell_temperature)
    end
    return (; core_trace, shell_trace)
end

times = (0:60:86400)u"s"
result = diurnal_trajectory(times, u"K"(20.0u"°C"), u"K"(20.0u"°C"))
println("core temperature range over the day: ", extrema(u"°C".(result.core_trace)))
println("shell temperature range over the day: ", extrema(u"°C".(result.shell_trace)))
println("max core-shell gap: ", maximum(abs.(result.core_trace .- result.shell_trace)))

# Uncomment to plot the core/shell divergence over the day:
# using Plots
# plot(uconvert.(u"hr", times), uconvert.(u"°C", result.core_trace); label="core")
# plot!(uconvert.(u"hr", times), uconvert.(u"°C", result.shell_trace); label="shell")
