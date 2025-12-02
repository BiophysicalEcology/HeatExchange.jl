abstract type AbstractMorphoModel end
        
abstract type AbstractMorphoParameters end

abstract type AbstractMorphoThresholds end

abstract type AbstractPhysioModel end
        
abstract type AbstractPhysioParameters end
        
abstract type AbstractPhysioThresholds end

abstract type AbstractBehavModel end 
        
abstract type AbstractBehavParameters end
        
abstract type AbstractBehavThresholds end

abstract type AbstractModelParameters end

"""
    Shape

Abstract supertype for the shape of the organism being modelled.
"""
abstract type Shape <: AbstractMorphoModel end

"""
    Insulation

Abstract supertype for the insulation of the organism being modelled.
"""
abstract type Insulation <: AbstractMorphoParameters end


"""
    Naked <: Insulation

    Naked()

Insulation trait for an organism without fur.
"""
struct Naked <: Insulation end

"""
    CompositeInsulation <: Insulation

    CompositeInsulation(layers)

A composite of insulation layers (e.g., fur, fat) for an organism.
"""
struct CompositeInsulation{T<:Tuple} <: Insulation
    layers::T
end
CompositeInsulation(i::Insulation) = CompositeInsulation((i,))
CompositeInsulation(is::Insulation...) = CompositeInsulation((is...,))
function geometry(shape, ins::CompositeInsulation)
    geometry(shape, ins.layers...)
end

"""
    Fur <: Insulation
    
    Fur(thickness)

Insulation trait for an organism with fur.
"""
struct Fur{T,D,R} <: Insulation
    thickness::T
    fibre_diameter::D
    fibre_density::R
end
"""
    Fat <: Insulation
    
    Fat(fraction, density)

Insulation trait for an organism with fat.
"""
struct Fat{F,D} <: Insulation
    fraction::F
    density::D
end

"""
    Geometry

    Geometry(volume, characteristic_dimension, length, area)

The geometry of an organism.
"""
struct Geometry{V,C,L,A} <: AbstractMorphoParameters
    volume::V
    characteristic_dimension::C
    length::L
    area::A
end

"""
    AbstractBody

Abstract supertype for organism bodies.
"""
abstract type AbstractBody <: AbstractMorphoParameters end

"""
    Body <: AbstractBody

    Body(shape::Shape, insulation::Insulation)
    Body(shape::Shape, insulation::Insulation, geometry::Geometry)

Physical dimensions of a body or body part that may or may note be insulated.
"""
struct Body{S<:Shape,I<:Insulation,G} <: AbstractBody
    shape::S
    insulation::I
    geometry::G
end
function Body(shape::Shape, insulation::Insulation)
    Body(shape, insulation, geometry(shape, insulation))
end

"""
    AbstractOrganism

Abstract supertype for organisms.
"""
abstract type AbstractOrganism end

# With some generic methods to get the params and body
body(o::AbstractOrganism) = o.body # gets the body from an object of type AbstractOrganism
#morphopars(o::AbstractOrganism) = o.morphopars # gets the morphological parameter traits from an object of type AbstractOrganism
integumentpars(o::AbstractOrganism) = o.integumentpars # gets the morphological parameter traits from an object of type AbstractOrganism
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
struct Organism{B<:Body,M<:AbstractMorphoParameters,P<:AbstractPhysioParameters,
    T<:AbstractBehavParameters, V<:AbstractBehavParameters} <: AbstractOrganism
    body::B
    integumentpars::M
    physiopars::P
    thermoregpars::T
    thermoregvars::V
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
