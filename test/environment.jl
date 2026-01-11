using HeatExchange
using ModelParameters
using Unitful

env_pars = EnvironmentalPars(
    α_ground = 0.8,
    ϵ_ground = 1.0,
    ϵ_sky = 1.0,
    elevation = 0.0u"m",
    fluid = 0,
    gas = GasFractions(0.2095, 0.0003, 0.79),
)
env_params = Model(env_pars)
env_pars = stripparams(env_params)

env_vars = EnvironmentalVars(;
    T_air = (273.15+20.0)u"K",
    T_sky = (273.15-5.0)u"K",
    T_ground = (273.15+30.0)u"K",
    T_substrate = (273.15+20.0)u"K",
    T_bush = (273.15+20.0)u"K",
    rh = 0.05,
    wind_speed = 1.0u"m/s",
    P_atmos = 101325.0u"Pa",
    zenith_angle = 20.0u"°",
    global_radiation = 1000.0u"W/m^2",
    k_substrate = 0.5u"W/m/K",
    diffuse_fraction = 0.15,
    shade = 0.0,
)
