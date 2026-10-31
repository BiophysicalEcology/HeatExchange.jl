# NLP strategy interface for IPOPT and other nonlinear solvers. BiophysicalBehaviour
# defines the one strategy, whose residuals are each part's `heat_balance`.

abstract type NLPStrategy end

"""
    nlp_pack(strategy::NLPStrategy, organism, environment,
             initial_skin_temperature, initial_insulation_temperature; smoothing)

Pack the fixed physics parameters for an NLP `strategy` once, before the solve.
Extended by the calling package for its concrete strategy (`MultipartNLP`);
returns a packed problem object the solver callbacks dispatch on.
"""
function nlp_pack end
