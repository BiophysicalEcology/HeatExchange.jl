# NLP strategy interface for IPOPT and other nonlinear solvers.
#
# HeatExchange owns the strategy *abstraction*; the calling package
# (BiophysicalBehaviour.jl) owns the solver policy (bounds, objective, variable
# unpacking) and supplies the concrete strategy. The one concrete strategy,
# `MultipartNLP`, is defined in BiophysicalBehaviour and extends `nlp_pack` here
# to pack the fixed per-part physics once before the solve; the per-iteration
# residuals go straight through `part_surface_residuals` (part_surface.jl), so no
# physics is duplicated across the package boundary.

"""
    NLPStrategy

Abstract supertype for formulations of the heat budget as a nonlinear program. The concrete strategy and its
[`nlp_pack`](@ref) method are defined by the calling package (BiophysicalBehaviour.jl), which owns the bounds,
objective and layout of variables. The physics of each evaluation comes from [`part_surface_residuals`](@ref).
"""
abstract type NLPStrategy end

"""
    nlp_pack(strategy::NLPStrategy, organism, environment,
             initial_skin_temperature, initial_insulation_temperature; smoothing)

Pack the fixed physics parameters for an NLP `strategy` once, before the solve.
Extended by the calling package for its concrete strategy (`MultipartNLP`);
returns a packed problem object the solver callbacks dispatch on.
"""
function nlp_pack end
