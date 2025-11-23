using HeatExchange
using ModelParameters
using Unitful

env_pars = EnvironmentalPars(
    α_substrate = 0.8,
    shade = 0.0,
    ϵ_substrate = 1.0,
    ϵ_sky = 1.0,
    elevation = 0.0u"m",
    fluid = 0,
    fN2 = 0.79,
    fO2 = 0.2095,
    fCO2 = 0.0003,
)
env_params = Model(env_pars)
env_pars = stripparams(env_params)

env_vars = EnvironmentalVars(
    T_air = (273.15+20.0)u"K",
    T_sky = (273.15-5.0)u"K",
    T_substrate = (273.15+30.0)u"K",
    rh = 0.05,
    wind_speed = 1.0u"m/s",
    P_atmos = 101325.0u"Pa",
    zenith_angle = 20.0u"°",
    k_substrate = 0.5u"W/m/K",
    Q_sol = 1000.0u"W/m^2",
    Q_dir = 964.177772475912u"W/m^2",
    Q_dif = 100.0u"W/m^2"
)

env_vars.T_substrate
