abstract type AbstractPhysioModel end
        
abstract type AbstractPhysioParameters end
        
abstract type AbstractPhysioThresholds end

"""
    MorphoPars <: AbstractMorphoParameters

A collection of morphological parameter functional traits for an organism.
"""
Base.@kwdef struct MorphoPars{F} <: AbstractMorphoParameters
    α_body_dorsal::F = Param(0.85, bounds=(0.0, 1.0))
    α_body_ventral::F = Param(0.85, bounds=(0.0, 1.0))
    ϵ_body_dorsal::F = Param(0.95, bounds=(0.0, 1.0))
    ϵ_body_ventral::F = Param(0.95, bounds=(0.0, 1.0))
    F_sky::F = Param(0.5, bounds=(0.0, 1.0))
    F_substrate::F = Param(0.5, bounds=(0.0, 1.0))
    eye_fraction::F = Param(0.0, bounds=(0.0, 1.0))
    skin_wetness::F = Param(0.0, bounds=(0.0, 1.0))
    conduction_fraction::F = Param(0.0, bounds=(0.0, 1.0))
    ventral_fraction::F = Param(0.5, bounds=(0.0, 1.0))
end

"""
    InsulationPars <: AbstractMorphoParameters

A collection of insulation parameter functional traits for an organism.
"""
Base.@kwdef struct InsulationPars{FDD,FDV,FLD,FLV,IDD,IDV,FRD,FRV,IRD,IRV,IDC,FC,LDF} <: AbstractMorphoParameters
    fibre_diameter_dorsal::FDD
    fibre_diameter_ventral::FDV
    fibre_length_dorsal::FLD
    fibre_length_ventral::FLV
    insulation_depth_dorsal::IDD
    insulation_depth_ventral::IDV
    fibre_density_dorsal::FRD 
    fibre_density_ventral::FRV 
    insulation_reflectance_dorsal::IRD
    insulation_reflectance_ventral::IRV
    insulation_depth_compressed::IDC
    fibre_conductivity::FC
    longwave_depth_fraction::LDF
end

"""
    PhysioPars <: AbstractParameterTrait

A collection of physiological parameter functional traits for an organism.
"""
Base.@kwdef struct PhysioPars{F,M,B,K} <: AbstractPhysioParameters
    fO2_extract::F = Param(0.20, bounds=(0.0, 1.0))
    rq::F = Param(0.8, bounds=(0.0, 1.0))
    M1::M = Param(0.013)
    M2::M = Param(0.8)
    M3::M = Param(0.038)
    M4::M = Param(0.0)
    pant::B = Param(1.0)
    k_body::K = Param(0.9, units="W/m/K")
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
    T_core::T
    T_surf::T
    T_lung::T
    ψ_org::P
end
