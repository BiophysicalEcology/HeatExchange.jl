using Revise
using ModelParameters
using Unitful

org_params = Model(OrganismalPars())
org_pars = stripparams(org_params)
org_pars.ϵ_org

org_vars = OrganismalVars()
org_vars.T_surf