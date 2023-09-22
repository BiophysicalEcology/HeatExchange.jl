using HeatExchange
using ModelParameters
using Unitful
using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

env_params = Model(EnvironmentalPars())
env_pars = stripparams(env_params)

env_vars = EnvironmentalVars()
env_vars.T_sub
