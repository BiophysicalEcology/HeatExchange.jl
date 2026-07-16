# Standalone one-lump transient body temperature, no behavioral loop. Mimics onelump.R's
# own example: heating trajectories of small vs large ellipsoids under otherwise
# identical conditions show how the thermal time constant scales with body size — small
# objects track air temperature far more closely than large ones.
using HeatExchange
using BiophysicalGeometry
using Unitful

environment_pars = example_environment_pars()
internal_conduction = example_conduction_pars_internal()
times = (1:60:7200)u"s"

function heating_trajectory(; mass, air_temperature, wind_speed)
    body = Body(Ellipsoid(mass, 1000.0u"kg/m^3", 1.1, 1.1), Naked())
    environment_vars = example_environment_vars(;
        air_temperature=u"K"(air_temperature), global_radiation=500.0u"W/m^2",
        zenith_angle=20.0u"°", wind_speed,
    )
    kw = (;
        internal_conduction, posture=Intermediate(),
        body_absorptivity=0.85, emissivity=0.95,
        sky_view_factor=0.4, ground_view_factor=0.4,
        metabolic_heat_volumetric=0.0u"W/m^3",
    )
    ectotherm_onelump(times, u"K"(20.0u"°C"), body, environment_pars, environment_vars; kw...)
end

# small (5 g) vs large (500 g), same wind speed and air temperature
small = heating_trajectory(mass=5.0u"g", air_temperature=20.0u"°C", wind_speed=1.0u"m/s")
large = heating_trajectory(mass=500.0u"g", air_temperature=20.0u"°C", wind_speed=1.0u"m/s")
println("5 g time constant:   ", u"minute"(small.time_constant))
println("500 g time constant: ", u"minute"(large.time_constant))

# same comparison but with calmer, warmer air for the small body (favours it further)
small_calm = heating_trajectory(mass=5.0u"g", air_temperature=25.0u"°C", wind_speed=0.5u"m/s")
println("5 g, calm+warm time constant: ", u"minute"(small_calm.time_constant))

# Uncomment to reproduce onelump.R's size-comparison plot:
# using Plots
# plot(uconvert.(u"hr", times), uconvert.(u"°C", small.core_temperature); label="5 g")
# plot!(uconvert.(u"hr", times), uconvert.(u"°C", large.core_temperature); label="500 g", linestyle=:dash)
