abstract type AbstractPhysiologyModel end

abstract type MetabolicRateEquation <: AbstractPhysiologyModel end
abstract type OxygenJoulesConversion <: AbstractPhysiologyModel end
        
abstract type AbstractPhysiologyParameters end

abstract type AbstractMorphologyParameters end
        
abstract type AbstractModelParameters end

abstract type AbstractFunctionalTraits end

shapepars(t::AbstractFunctionalTraits) = stripparams(t.shape_pars)
insulationpars(t::AbstractFunctionalTraits) = stripparams(t.insulation_pars)
conductionpars_external(t::AbstractFunctionalTraits) = stripparams(t.conduction_pars_external)
conductionpars_internal(t::AbstractFunctionalTraits) = stripparams(t.conduction_pars_internal)
convectionpars(t::AbstractFunctionalTraits) = stripparams(t.convection_pars)
radiationpars(t::AbstractFunctionalTraits) = stripparams(t.radiation_pars)
evaporationpars(t::AbstractFunctionalTraits) = stripparams(t.evaporation_pars)
hydraulicpars(t::AbstractFunctionalTraits) = stripparams(t.hydraulic_pars)
respirationpars(t::AbstractFunctionalTraits) = stripparams(t.respiration_pars)
metabolismpars(t::AbstractFunctionalTraits) = stripparams(t.metabolism_pars)

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
shapepars(o::AbstractOrganism) = shapepars(traits(o))
insulationpars(o::AbstractOrganism) = insulationpars(traits(o))
conductionpars_external(o::AbstractOrganism) = conductionpars_external(traits(o))
conductionpars_internal(o::AbstractOrganism) = conductionpars_internal(traits(o))
convectionpars(o::AbstractOrganism) = convectionpars(traits(o))
radiationpars(o::AbstractOrganism) = radiationpars(traits(o))
evaporationpars(o::AbstractOrganism) = evaporationpars(traits(o))
hydraulicpars(o::AbstractOrganism) = hydraulicpars(traits(o))
respirationpars(o::AbstractOrganism) = respirationpars(traits(o))
metabolismpars(o::AbstractOrganism) = metabolismpars(traits(o))

"""
    Organism <: AbstractOrganism

    Organism(body, traits)

A concrete implementation of `AbstractOrganism`, it accepts an
[`AbstractBody`](@ref) and [`AbstractFunctionalTraits`](@ref) object.
"""
struct Organism{
        B<:AbstractBody,
        T<:AbstractFunctionalTraits} <: AbstractOrganism
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

    - `ψ_org` — Body water potential (determines humidity at skin surface 
    and liquid water exchange) (J/kg).
Variables for an [`AbstractOrganism`](@ref) model.
"""
Base.@kwdef mutable struct OrganismalVars{TC,TS,TI,TL,P} <: AbstractOrganismalVars
    T_core::TC
    T_skin::TS
    T_insulation::TI = T_core
    T_lung::TL
    ψ_org::P
end
