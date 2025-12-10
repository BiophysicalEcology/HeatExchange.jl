abstract type AbstractGeometryModel end

"""
    AbstractShape

Abstract supertype for the shape of the organism being modelled.
"""
abstract type AbstractShape <: AbstractGeometryModel end

"""
    AbstractInsulation

Abstract supertype for the insulation of the organism being modelled.
"""
abstract type AbstractInsulation <: AbstractGeometryModel end

"""
    Naked <: AbstractInsulation

    Naked()

Insulation trait for an organism without fur.
"""
struct Naked <: AbstractInsulation end

"""
    Fur <: AbstractInsulation
    
    Fur(thickness)

Insulation trait for an organism with fur.
"""
struct Fur{T,D,R} <: AbstractInsulation
    thickness::T
    fibre_diameter::D
    fibre_density::R
end

# TODO
# """
#     Feathers <: AbstractInsulation
    
#     Feathers(thickness)

# Insulation trait for an organism with feathers.
# """
# struct Feathers{T,D,R} <: AbstractInsulation
#     thickness::T
#     fibre_diameter::D
#     fibre_density::R
# end

"""
    Fat <: AbstractInsulation
    
    Fat(fraction, density)

Insulation trait for an organism with fat.
"""
struct Fat{F,D} <: AbstractInsulation
    fraction::F
    density::D
end

"""
    CompositeInsulation <: AbstractInsulation

    CompositeInsulation(layers)

A composite of insulation layers (e.g., fur, fat) for an organism.
"""
struct CompositeInsulation{T<:Tuple} <: AbstractInsulation
    layers::T
end
CompositeInsulation(i::AbstractInsulation) = CompositeInsulation((i,))
CompositeInsulation(is::AbstractInsulation...) = CompositeInsulation((is...,))
function geometry(shape, ins::CompositeInsulation)
    geometry(shape, ins.layers...)
end

abstract type AbstractGeometryPars end

"""
    AbstractBody

Abstract supertype for organism bodies.
"""
abstract type AbstractBody <: AbstractGeometryPars end

"""
    Geometry

    Geometry(volume, characteristic_dimension, length, area)

The geometry of an organism.
"""
struct Geometry{V,C,L,A} <: AbstractGeometryPars
    volume::V
    characteristic_dimension::C
    length::L
    area::A
end

"""
    Body <: AbstractBody

    Body(shape::AbstractShape, insulation::AbstractInsulation)
    Body(shape::AbstractShape, insulation::AbstractInsulation, geometry::AbstractGeometryPars)

Physical dimensions of a body or body part that may or may note be insulated.
"""
struct Body{
        S<:AbstractShape,
        I<:AbstractInsulation,
        G<:AbstractGeometryPars} <: AbstractBody
    shape::S
    insulation::I
    geometry::G
end

abstract type SolarOrientation <: AbstractGeometryPars end

struct NormalToSun <: SolarOrientation end
struct ParallelToSun <: SolarOrientation end
struct Intermediate <: SolarOrientation end

# constructors and functions

function Body(shape::AbstractShape, insulation::AbstractInsulation)
    Body(shape, insulation, geometry(shape, insulation))
end

shape(body::AbstractBody) = body.shape
insulation(body::AbstractBody) = body.insulation
geometry(body::AbstractBody) = body.geometry
surface_area(body::AbstractBody) = surface_area(shape(body), body)

# functions to extract appropriate surface areas from different objects

total_area(body::AbstractBody) = total_area(shape(body), insulation(body), body)
skin_area(body::AbstractBody) = skin_area(shape(body), insulation(body), body)
evaporation_area(body::AbstractBody) = evaporation_area(shape(body), insulation(body), body)

# for composite insulation cases (fat and fur/feathers)
outer_insulation(ins::AbstractInsulation) = ins
outer_insulation(ins::CompositeInsulation) = begin
    # find fur layer if present
    fur_layer = findlast(i -> i isa Fur, ins.layers)
    if fur_layer !== nothing
        ins.layers[fur_layer]
    else
        # otherwise the last layer
        ins.layers[end]
    end
end

total_area(shape, ins::CompositeInsulation, body) =
    total_area(shape, outer_insulation(ins), body)
skin_area(shape, ins::CompositeInsulation, body) =
    skin_area(shape, outer_insulation(ins), body)
evaporation_area(shape, ins::CompositeInsulation, body) =
    evaporation_area(shape, outer_insulation(ins), body)

# functions to get the appropriate radii
skin_radius(body::AbstractBody) = skin_radius(shape(body), insulation(body), body)
insulation_radius(body::AbstractBody) = insulation_radius(shape(body), insulation(body), body)
flesh_radius(body::AbstractBody) = flesh_radius(shape(body), insulation(body), body)

# functions to compute silhouette area as a function of zenith angle (passive) or 
# orientation (thermoregulating)

"""
    silhouette_area(body::AbstractBody, θ)

Calculates the silhouette (projected) area of an object given a solar zenith angle, θ.

"""
silhouette_area(body::AbstractBody, θ) = silhouette_area(shape(body), insulation(body), body, θ)

"""
    silhouette_area(body::AbstractBody, orientation)

Calculates the silhouette (projected) area of an object given an orientation towards the sun (normal, parallel or inbetween).

"""
silhouette_area(body::AbstractBody) = silhouette_area(shape(body), insulation(body), body)

# Orientation-specific implementations
silhouette_area(body::AbstractBody, ::NormalToSun) =
    silhouette_area(body).normal

silhouette_area(body::AbstractBody, ::ParallelToSun) =
    silhouette_area(body).parallel

silhouette_area(body::AbstractBody, ::Intermediate) =
    (silhouette_area(body).normal + silhouette_area(body).parallel) * 0.5

    #TODO make this 'insulation_area'
function hair_area(fibre_diameter, fibre_density, skin)
    π * (fibre_diameter / 2) ^ 2 * (fibre_density * skin)
end