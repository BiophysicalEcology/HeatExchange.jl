using Revise
using ModelParameters
using Unitful
using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

include("../src/organism.jl")

org_params = Model(OrganismalPars())
org_pars = stripparams(org_params)
org_pars.ϵ_org

org_vars = OrganismalVars()
org_vars.T_surf