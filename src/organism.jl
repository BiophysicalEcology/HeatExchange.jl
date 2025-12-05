abstract type AbstractPhysioModel end

abstract type MetabolicRateEquation <: AbstractPhysioModel end
abstract type OxygenJoulesConversion <: AbstractPhysioModel end
        
abstract type AbstractPhysioParameters end

abstract type AbstractMorphoParameters end
        
abstract type AbstractBehavParameters end
        
#abstract type AbstractBehavThresholds end

abstract type AbstractModelParameters end

abstract type AbstractFunctionalTraits end

struct Traits{
        IN<:AbstractPhysioParameters,
        CE<:AbstractMorphoParameters,
        CI<:AbstractMorphoParameters,
        RA<:AbstractMorphoParameters,
        CO<:AbstractMorphoParameters,
        EV<:AbstractMorphoParameters,
        RE<:AbstractPhysioParameters,
        ME<:AbstractPhysioParameters,
 } <: AbstractFunctionalTraits
    insulationpars::IN
    conductionpars_external::CE
    conductionpars_internal::CI
    radiationpars::RA
    convectionpars::CO
    evaporationpars::EV
    respirationpars::RE
    metabolicpars::ME
end

"""
    AbstractOrganism

Abstract supertype for organisms.
"""
abstract type AbstractOrganism end

# With some generic methods to get the params and body
body(o::AbstractOrganism) = o.body # gets the body from an object of type AbstractOrganism
#morphopars(o::AbstractOrganism) = o.morphopars # gets the morphological parameter traits from an object of type AbstractOrganism
radiationpars(o::AbstractOrganism) = o.radiationpars # gets the radiation parameter traits from an object of type AbstractOrganism
evaporationpars(o::AbstractOrganism) = o.evaporationpars # gets the evaporation parameter traits from an object of type AbstractOrganism
physiopars(o::AbstractOrganism) = o.physiopars # gets the physioloigcal parameter traits from an object of type AbstractOrganism
thermoregpars(o::AbstractOrganism) = o.thermoregpars # gets the physioloigcal parameter traits from an object of type AbstractOrganism
thermoregvars(o::AbstractOrganism) = o.thermoregvars # gets the physioloigcal parameter traits from an object of type AbstractOrganism
shape(o::AbstractOrganism) = shape(body(o)) # gets the shape from an object of type AbstractOrganism
insulation(o::AbstractOrganism) = insulation(body(o)) # gets the insulation from an object of type AbstractOrganism

"""
    Organism <: AbstractOrganism

    Organism(body, traits)

A concrete implementation of `AbstractOrganism`, it accepts an
[`AbstractBody`](@ref) and [`AbstractFunctionalTraits`](@ref) object.
"""
struct Organism{
        B<:AbstractBody,
        IN<:AbstractPhysioParameters,
        CE<:AbstractMorphoParameters,
        CI<:AbstractMorphoParameters,
        RA<:AbstractMorphoParameters,
        CO<:AbstractMorphoParameters,
        EV<:AbstractMorphoParameters,
        RE<:AbstractPhysioParameters,
        ME<:AbstractPhysioParameters,
        } <: AbstractOrganism
    body::B

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
