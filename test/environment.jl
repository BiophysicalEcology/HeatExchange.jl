using Revise
using ModelParameters
using Unitful
using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

include("../src/environment.jl")

env_params = Model(EnvironmentalPars())
env_pars = stripparams(env_params)
env_pars.P_atmos

env_vars = EnvironmentalVars()
env_vars.Tsub