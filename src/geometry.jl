abstract type AbstractMorphoModel end
        
abstract type AbstractMorphoParameters end

abstract type AbstractMorphoThresholds end

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
    Fur <: Insulation
    
    Fur(thickness)

Insulation trait for an organism with fur.
"""
struct Fur{T} <: Insulation
    thickness::T
end

"""
    Geometry

    Geometry(volume, characteristic_dimension, lengths, area)

The geometry of an organism.
"""
struct Geometry{V,C,L,A} <: AbstractMorphoParameters
    volume::V
    characteristic_dimension::C
    lengths::L
    area::A
end

"""
    AbstractBody

Abstract supertype for organism bodies.
"""
abstract type AbstractBody <: AbstractMorphoParameters end

shape(body::AbstractBody) = body.shape
insulation(body::AbstractBody) = body.insulation
geometry(body::AbstractBody) = body.geometry
calc_area(body::AbstractBody) = calc_area(shape(body), body)

"""
    calc_silhouette_area(body::AbstractBody, θ)

Calculates the silhouette (projected) area of a cylinder given a solar zenith angle, θ.
Calculates the silhouette (projected) area of a cylinder.
Equation from Fig. 11.6 in Campbell, G. S., & Norman, J. M.
(1998). Environmental Biophysics. Springer.
"""
calc_silhouette_area(body::AbstractBody, θ) = calc_silhouette_area(shape(body), body, θ)

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
    Plate <: Shape

A flat plate-shaped organism shape.
"""
struct Plate{M,D,B,C} <: Shape
    mass::M
    density::D
    b::B
    c::C
end

function geometry(shape::Plate, ::Naked)
    volume = shape.mass / shape.density
    characteristic_dimension = volume^(1 / 3)
    a = (volume / (shape.b * shape.c))^(1 / 3)
    b = shape.b * a
    c = shape.c * a
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    area = calc_area(shape, a, b, c)
    return Geometry(volume, characteristic_dimension, (length1, length2, length3), area)
end
# function geometry(shape::Plate, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, characteristic_dimension, (length1, length2, length3), area, sil_area)
# end

function calc_area(shape::Plate, body)
    a = body.geometry.lengths[1] / 2
    b = body.geometry.lengths[2] / 2
    c = body.geometry.lengths[3] / 2
    calc_area(shape, a, b, c)
end
calc_area(shape::Plate, a, b, c) = a * b * 2 + a * c * 2 + b * c * 2

"""
    Cylinder <: Shape

A cylindrical organism shape.
"""
struct Cylinder{M,D,B} <: Shape
    mass::M
    density::D
    b::B
end

function geometry(shape::Cylinder, ::Naked)
    volume = shape.mass / shape.density
    characteristic_dimension = volume^(1 / 3)
    r = (volume / (shape.b * π * 2))^(1 / 3)
    length1 = shape.b * r * 2
    length2 = 2 * r
    area = calc_area(shape, r, length1)
    return Geometry(volume, characteristic_dimension, (length1, length2), area)
end
# function geometry(shape::Cylinder, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, characteristic_dimension, (length1, length2, length3), area, sil_area)
# end

function calc_area(shape::Cylinder, body::AbstractBody)
    r = body.geometry.lengths[2] / 2
    l = body.geometry.lengths[1]
    calc_area(shape, r, l)
end
calc_area(shape::Cylinder, r, l) = 2 * π * r * l + 2 * π * r^2

function calc_silhouette_area(shape::Cylinder, body::AbstractBody, θ)
    r = body.geometry.lengths[2] / 2
    l = body.geometry.lengths[1]
    return calc_silhouette_area(shape, r, l, θ)
end
calc_silhouette_area(shape::Cylinder, r, l, θ) = 2 * r * l * sin(θ) + π * r^2 * cos(θ)

"""
    Ellipsoid <: Shape

An ellipsoidal organism shape.
"""
struct Ellipsoid{M,D,B,C} <: Shape
    mass::M
    density::D
    b::B
    c::C
end

function geometry(shape::Ellipsoid, ::Naked)
    volume = shape.mass / shape.density
    characteristic_dimension = volume^(1 / 3)
    a = ((3 / 4)* Unitful.ustrip(volume) / (π * shape.b * shape.c)) ^ (1 / 3)
    b = a * shape.b
    c = a * shape.c
    length1 = a * 2m
    length2 = b * 2m
    length3 = c * 2m
    area = calc_area(shape, a, b, c)m^2
    return Geometry(volume, characteristic_dimension, (length1, length2, length3), area)
end
# function geometry(shape::Ellipsoid, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, characteristic_dimension, (length1, length2, length3), area)
# end

function calc_area(shape::Ellipsoid, body::Body)
    a = body.geometry.lengths[1] / 2
    b = body.geometry.lengths[2] / 2
    c = body.geometry.lengths[3] / 2
    return calc_area(shape, a, b, c)
end
function calc_area(shape::Ellipsoid, a, b, c)
    #e = ((a ^ 2 - c ^ 2) ^ 0.5 ) / a # eccentricity
    #2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)
    p =  1.6075
    return(4 * π * (((a ^ p * b ^ p + a ^ p * c ^ p + b ^ p * c ^ p)) / 3) ^ (1 / p))
end

# TODO: should we use this rather than the one below that ignores theta?
# """
#     calc_silhouette_area

# Calculates the silhouette (projected) area of a prolate spheroid.
# """
# function calc_silhouette_area(shape::Ellipsoid, a, b, c, θ)
#     a2 = cos(90) ^ 2 * (cos(θ) ^ 2 / a ^ 2 + sin(θ) ^ 2 / b ^ 2) + sin(90) ^ 2 / c ^ 2
#     twohh = 2 * cos(90) * sin(90) * cos(θ) * (1 / b ^ 2 - 1 / a ^ 2)
#     b2 = sin(θ) ^ 2 / a ^ 2 + cos(θ) ^ 2 / b ^ 2
#     θ2 = 0.5 * atan(twohh, a2 - b2)
#     sps = sin(θ2)
#     cps = cos(θ2)
#     a3 = cps * (a2 * cps + twohh * sps) + b2 * sps * sps
#     b3 = sps * (a2 * sps - twohh * cps) + b2 * cps * cps
#     semax1 = 1 / sqrt(a3)
#     semax2 = 1 / sqrt(b3)
#     π * semax1  * semax2
# end
# function calc_silhouette_area(shape::Ellipsoid, body, θ)
#     a = body.geometry.lengths[1] / 2
#     b = body.geometry.lengths[2] / 2
#     c = body.geometry.lengths[3] / 2
#     calc_silhouette_area(shape, a, b, c, θ)
# end
function calc_silhouette_area(shape::Ellipsoid, body, θ)
    a = body.geometry.lengths[1] / 2
    b = body.geometry.lengths[2] / 2
    c = body.geometry.lengths[3] / 2
    return calc_silhouette_area(shape, a, b, c, θ)
end
calc_silhouette_area(shape::Ellipsoid, a, b, c, θ) = max(π * a * c, π * b * c)

"""
    LeopardFrog <: Shape

A frog-shaped organism. Based on the leopard frog (Tracy 1976 Ecol. Monog.).
"""
struct LeopardFrog{M,D} <: Shape
    mass::M
    density::D
end

function geometry(shape::LeopardFrog, ::Naked)
    volume = shape.mass / shape.density
    characteristic_dimension = volume^(1 / 3)
    area = calc_area(shape, body)
    return Geometry(volume, characteristic_dimension, nothing, area)
end

calc_area(shape::LeopardFrog, body) = calc_area(shape)
function calc_area(shape::LeopardFrog)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    Unitful.uconvert(u"m^2", (12.79 * Unitful.ustrip(mass_g) ^ 0.606)u"cm^2") # eq in Fig. 5
end

calc_silhouette_area(shape::LeopardFrog, body, θ) = calc_silhouette_area(shape::LeopardFrog, θ)
function calc_silhouette_area(shape::LeopardFrog, θ)
    area = calc_area(shape)
    pct = 1.38171e-6 * θ ^ 4 - 1.93335e-4 * θ ^ 3 + 4.75761e-3 * θ ^ 2 - 0.167912 * θ + 45.8228
    return pct * area / 100
end

"""
    DesertIguana <: Shape

A lizard-shaped organism. Based on the desert iguana (Porter and Tracy 1984).
"""
struct DesertIguana{M,D} <: Shape
    mass::M
    density::D
end

function geometry(shape::DesertIguana, ::Naked)
    volume = shape.mass / shape.density
    characteristic_dimension = volume^(1 / 3)
    area = calc_area(shape)
    return Geometry(volume, characteristic_dimension, nothing, area)
end

function calc_area(shape::DesertIguana)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    return Unitful.uconvert(u"m^2", (10.4713 * Unitful.ustrip(mass_g) ^ 0.688)u"cm^2")
end
