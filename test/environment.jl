using HeatExchange
using ModelParameters
using Unitful

environment_pars = EnvironmentalPars(;
    ground_albedo=0.8,
    ground_emissivity=1.0,
    sky_emissivity=1.0,
    elevation=0.0u"m",
    fluid=0,
    gas_fractions=GasFractions(0.2095, 0.0003, 0.79),
)
env_params = Model(environment_pars)
environment_pars = stripparams(env_params)

environment_vars = EnvironmentalVars(;
    air_temperature=(273.15+20.0)u"K",
    sky_temperature=(273.15-5.0)u"K",
    ground_temperature=(273.15+30.0)u"K",
    substrate_temperature=(273.15+20.0)u"K",
    bush_temperature=(273.15+20.0)u"K",
    relative_humidity=0.05,
    wind_speed=1.0u"m/s",
    atmospheric_pressure=101325.0u"Pa",
    zenith_angle=20.0u"°",
    global_radiation=1000.0u"W/m^2",
    substrate_conductivity=0.5u"W/m/K",
    diffuse_fraction=0.15,
    shade=0.0,
)
