abstract type AbstractPhysiologyModel end

abstract type MetabolicRateEquation <: AbstractPhysiologyModel end
abstract type OxygenJoulesConversion <: AbstractPhysiologyModel end

abstract type AbstractPhysiologyParameters end

abstract type AbstractMorphologyParameters end

abstract type AbstractModelParameters end

abstract type AbstractFunctionalTraits end

shape_pars(t::AbstractFunctionalTraits) = stripparams(t.shape_pars)
insulation_pars(t::AbstractFunctionalTraits) = stripparams(t.insulation_pars)
function conduction_pars_external(t::AbstractFunctionalTraits)
    stripparams(t.conduction_pars_external)
end
function conduction_pars_internal(t::AbstractFunctionalTraits)
    stripparams(t.conduction_pars_internal)
end
convection_pars(t::AbstractFunctionalTraits) = stripparams(t.convection_pars)
radiation_pars(t::AbstractFunctionalTraits) = stripparams(t.radiation_pars)
evaporation_pars(t::AbstractFunctionalTraits) = stripparams(t.evaporation_pars)
hydraulic_pars(t::AbstractFunctionalTraits) = stripparams(t.hydraulic_pars)
respiration_pars(t::AbstractFunctionalTraits) = stripparams(t.respiration_pars)
metabolism_pars(t::AbstractFunctionalTraits) = stripparams(t.metabolism_pars)
options(t::AbstractFunctionalTraits) = stripparams(t.options)

# TODO more specific subtypes
struct HeatExchangeTraits{
    SP<:AbstractShape,
    IN<:AbstractMorphologyParameters,
    CE<:AbstractMorphologyParameters,
    CI<:AbstractPhysiologyParameters,
    RA<:AbstractMorphologyParameters,
    CO<:AbstractMorphologyParameters,
    EV<:AbstractMorphologyParameters,
    HD<:AbstractPhysiologyParameters,
    RE<:AbstractPhysiologyParameters,
    ME<:AbstractPhysiologyParameters,
    OP<:AbstractModelParameters,
} <: AbstractFunctionalTraits
    shape_pars::SP
    insulation_pars::IN
    conduction_pars_external::CE
    conduction_pars_internal::CI
    radiation_pars::RA
    convection_pars::CO
    evaporation_pars::EV
    hydraulic_pars::HD
    respiration_pars::RE
    metabolism_pars::ME
    options::OP
end

"""
    AbstractOrganism

Abstract supertype for organisms.
"""
abstract type AbstractOrganism end

# With some generic methods to get the params and body
body(o::AbstractOrganism) = o.body # gets the body from an object of type AbstractOrganism
traits(o::AbstractOrganism) = o.traits
#shape(o::AbstractOrganism) = shape(body(o)) # gets the shape from an object of type AbstractOrganism
#insulation(o::AbstractOrganism) = insulation(body(o)) # gets the insulation from an object of type AbstractOrganism

# Forwarding methods from organism to traits
shape_pars(o::AbstractOrganism) = shape_pars(traits(o))
insulation_pars(o::AbstractOrganism) = insulation_pars(traits(o))
conduction_pars_external(o::AbstractOrganism) = conduction_pars_external(traits(o))
conduction_pars_internal(o::AbstractOrganism) = conduction_pars_internal(traits(o))
convection_pars(o::AbstractOrganism) = convection_pars(traits(o))
radiation_pars(o::AbstractOrganism) = radiation_pars(traits(o))
evaporation_pars(o::AbstractOrganism) = evaporation_pars(traits(o))
hydraulic_pars(o::AbstractOrganism) = hydraulic_pars(traits(o))
respiration_pars(o::AbstractOrganism) = respiration_pars(traits(o))
metabolism_pars(o::AbstractOrganism) = metabolism_pars(traits(o))
options(o::AbstractOrganism) = options(traits(o))

"""
    Organism <: AbstractOrganism

    Organism(body, traits)

A concrete implementation of `AbstractOrganism`, it accepts an
[`AbstractBody`](@ref) and [`AbstractFunctionalTraits`](@ref) object.
"""
struct Organism{B<:AbstractBody,T<:AbstractFunctionalTraits} <: AbstractOrganism
    body::B
    traits::T
end

# TODO use this as a container for outputs? Or remove?
"""
    AbstractOrganismalVars

Abstract supertype for organismal variables.
"""
abstract type AbstractOrganismalVars end

"""
    OrganismalVars <: AbstractOrganismalVars

    - `water_potential` — Body water potential (determines humidity at skin surface
    and liquid water exchange) (J/kg).
Variables for an [`AbstractOrganism`](@ref) model.
"""
Base.@kwdef mutable struct OrganismalVars{TC,TS,TI,TL,P} <: AbstractOrganismalVars
    core_temperature::TC
    skin_temperature::TS
    insulation_temperature::TI = core_temperature
    lung_temperature::TL
    water_potential::P
end
