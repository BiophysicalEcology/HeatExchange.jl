abstract type AbstractPhysioModel end
        
abstract type AbstractPhysioParameters end
        
abstract type AbstractPhysioThresholds end

# """
#     AbstractFunctionalTraits 

# An abstract super type for organism functional traits.
# """
# abstract type AbstractFunctionalTraits end

# """
#     FunctionalTraits <: AbstractFunctionalTraits

# A collection of functional traits for an organism.
# """
# Base.@kwdef struct FunctionalTraits{F,K,M,B} <: AbstractFunctionalTraits
#     α_org_dorsal::F = Param(0.85, bounds=(0.2, 1.0))
#     α_org_ventral::F = Param(0.85, bounds=(0.2, 1.0))
#     ϵ_org_dorsal::F = Param(0.95, bounds=(0.1, 1.0))
#     ϵ_org_ventral::F = Param(0.95, bounds=(0.1, 1.0))
#     F_sky::F = Param(0.4, bounds=(0.3, 0.5))
#     F_sub::F = Param(0.4, bounds=(0.3, 0.5))
#     p_eyes::F = Param(0.0, bounds=(0.0, 4e-4))
#     fO2_extract::F = Param(0.20, bounds=(0.10, 0.30))
#     k_body::K = Param(0.5, bounds=(0.412, 2.8), units=u"W/m/K")
#     rq::F = Param(0.8, bounds=(0.7, 0.9))
#     M1::M = Param(0.013, bounds=(0.01, 0.02))
#     M2::M = Param(0.8, bounds=(0.7, 0.9))
#     M3::M = Param(0.038, bounds=(0.02, 0.04))
#     p_wet::F = Param(0.001, bounds=(0.0, 1.0))
#     p_cond::F = Param(0.1, bounds=(0.0, 1.0))
#     pant::B = Param(1.0, bounds=(1.0, 10.0))
# end


"""
    MorphoPars <: AbstractParameterTrait

A collection of morphological parameter functional traits for an organism.
"""
Base.@kwdef struct MorphoPars{F,K} <: AbstractMorphoParameters
    α_org_dorsal::F = Param(0.85, bounds=(0.2, 1.0))
    α_org_ventral::F = Param(0.85, bounds=(0.2, 1.0))
    ϵ_org_dorsal::F = Param(0.95, bounds=(0.1, 1.0))
    ϵ_org_ventral::F = Param(0.95, bounds=(0.1, 1.0))
    F_sky::F = Param(0.4, bounds=(0.3, 0.5))
    F_sub::F = Param(0.4, bounds=(0.3, 0.5))
    k_body::K = Param(0.5, bounds=(0.412, 2.8), units=u"W/m/K")
    p_eyes::F = Param(0.0, bounds=(0.0, 4e-4))
    p_wet::F = Param(0.001, bounds=(0.0, 1.0))
    p_cond::F = Param(0.1, bounds=(0.0, 1.0))
end

"""
    PhysioPars <: AbstractParameterTrait

A collection of physiological parameter functional traits for an organism.
"""
Base.@kwdef struct PhysioPars{F,M,B} <: AbstractPhysioParameters
    fO2_extract::F = Param(0.20, bounds=(0.10, 0.30))
    rq::F = Param(0.8, bounds=(0.7, 0.9))
    M1::M = Param(0.013, bounds=(0.01, 0.02))
    M2::M = Param(0.8, bounds=(0.7, 0.9))
    M3::M = Param(0.038, bounds=(0.02, 0.04))
    pant::B = Param(1.0, bounds=(1.0, 10.0))
end

"""
    AbstractOrganism

Abstract supertype for organisms.
"""
abstract type AbstractOrganism end

# With some generic methods to get the params and body
body(o::AbstractOrganism) = o.body # gets the body from an object of type AbstractOrganism
morphopars(o::AbstractOrganism) = o.morphopars # gets the morphological parameter traits from an object of type AbstractOrganism
physiopars(o::AbstractOrganism) = o.physiopars # gets the physioloigcal parameter traits from an object of type AbstractOrganism
shape(o::AbstractOrganism) = shape(body(o)) # gets the shape from an object of type AbstractOrganism
insulation(o::AbstractOrganism) = insulation(body(o)) # gets the insulation from an object of type AbstractOrganism

"""
    Organism <: AbstractOrganism

    Organism(body, traits)

A concrete implementation of `AbstractOrganism`, it accepts an
[`AbstractBody`](@ref) and [`AbstractFunctionalTraits`](@ref) object.
"""
struct Organism{B<:Body,M<:AbstractMorphoParameters,P<:AbstractPhysioParameters} <: AbstractOrganism
    body::B
    morphopars::M
    physiopars::P
end

"""
    AbstractOrganismalVars

Abstract supertype for organismal variables.
"""
abstract type AbstractOrganismalVars end

"""
    OrganismalVars <: AbstractOrganismalVars

Variables for an [`AbstractOrganism`](@ref) model.
"""
Base.@kwdef mutable struct OrganismalVars{T,P} <: AbstractOrganismalVars
    T_core::T = K(20.0°C)
    T_surf::T = K(20.0°C)
    T_lung::T = K(20.0°C)
    ψ_org::P = -707.0J/kg
end
