# Standalone one-lump transient body temperature, no behavioral loop. Heating trajectories
# of small vs large ellipsoids show how the thermal time constant scales with body size.
# No closed form (onelump's Naked branch always includes evaporative/respiratory heat
# loss now), so this integrates the derivative form directly (RK4).
using HeatExchange
using BiophysicalGeometry
using Unitful

function rk4_step(f, u, t, dt)
    k1 = f(u, t)
    k2 = f(u + dt / 2 * k1, t + dt / 2)
    k3 = f(u + dt / 2 * k2, t + dt / 2)
    k4 = f(u + dt * k3, t + dt)
    return u + dt / 6 * (k1 + 2k2 + 2k3 + k4)
end

environment_pars = example_environment_pars()
times = (1:60:7200)u"s"

function heating_trajectory(; mass, air_temperature, wind_speed)
    shape_pars = Ellipsoid(mass, 1000.0u"kg/m^3", 1.1, 1.1)
    body = Body(shape_pars, Naked())
    traits = example_heat_exchange_traits(; shape_pars)
    organism = Organism(body, traits)
    environment_vars = example_environment_vars(;
        air_temperature=u"K"(air_temperature), global_radiation=500.0u"W/m^2",
        zenith_angle=20.0u"°", wind_speed,
    )
    e = (; environment_pars, environment_vars)

    f(u, t) = onelump(u, t, organism, e).core_temperature_rate
    core_temperature_init = u"K"(20.0u"°C")
    core_temperature = core_temperature_init
    trace = [core_temperature]
    for i in 2:length(times)
        dt = times[i] - times[i - 1]
        core_temperature = rk4_step(f, core_temperature, times[i - 1], dt)
        push!(trace, core_temperature)
    end
    equilibrium = solve_temperature(organism, e).core_temperature
    return (; trace, equilibrium, core_temperature_init)
end

# small (5 g) vs large (500 g), same wind speed and air temperature
small = heating_trajectory(mass=5.0u"g", air_temperature=20.0u"°C", wind_speed=1.0u"m/s")
large = heating_trajectory(mass=500.0u"g", air_temperature=20.0u"°C", wind_speed=1.0u"m/s")

# fraction of the way from initial to equilibrium temperature reached after 30 minutes —
# the derivative-form analogue of onelump.R's closed-form time_constant comparison
fraction_to_equilibrium(traj, t) = (traj.trace[argmin(abs.(times .- t))] - traj.core_temperature_init) /
                                    (traj.equilibrium - traj.core_temperature_init)
println("5 g: ", round(100 * fraction_to_equilibrium(small, 30u"minute"); digits=1), "% of the way to equilibrium after 30 min")
println("500 g: ", round(100 * fraction_to_equilibrium(large, 30u"minute"); digits=1), "% of the way to equilibrium after 30 min")

# same comparison but with calmer, warmer air for the small body (favours it further)
small_calm = heating_trajectory(mass=5.0u"g", air_temperature=25.0u"°C", wind_speed=0.5u"m/s")
println("5 g, calm+warm: ", round(100 * fraction_to_equilibrium(small_calm, 30u"minute"); digits=1), "% of the way to equilibrium after 30 min")

# Uncomment to reproduce onelump.R's size-comparison plot:
# using Plots
# plot(uconvert.(u"hr", times), uconvert.(u"°C", small.trace); label="5 g")
# plot!(uconvert.(u"hr", times), uconvert.(u"°C", large.trace); label="500 g", linestyle=:dash)
