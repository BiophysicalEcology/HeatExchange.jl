using Unitful
using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

"""
    Shape

Abstract type for the shape of the object being modelled.
"""
abstract type Shape end

abstract type Insulation end

abstract type AbstractBody end

struct Naked <: Insulation end

struct Fur{T} <: Insulation
    thickness::T
end

struct Geometry{V,C,L,A}
    volume::V
    characteristic_dimension::C
    lengths::L
    area::A
end

"""
    Body <: AbstractBody

Physical dimensions of a body or body part that may or may note be insulated.
"""
struct Body{S<:Shape,I<:Insulation,G} <: AbstractBody
    shape::S
    insulation::I
    geometry::G
end

# Add a Body constructor from shape, insulation, mass and density
function Body(shape::Shape, insulation::Insulation)
    Body(shape, insulation, geometry(shape, insulation))
end

# struct BodyContext{B<:Body,S,X}
#     body::B
#     silhouette:S
#     #other_time_location_dependent_thing::X
# end
# function BodyContext(body::Body, environment)
#     silhouette = calc_silhouette(body, environment)
#     #other_thing = calc_thing(body, environment)
#     BodyContext(body, silhouette)
# end


"""
    Plate <: Shape

A flat plate-shaped organism.
"""
struct Plate{M,D,B,C} <: Shape
    mass::M
    density::D
    b::B
    c::C
end

function geometry(shape::Plate, ::Naked)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    a = (volume / (shape.b * shape.c))^(1 / 3)
    b = shape.b * a
    c = shape.c * a
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    area = calc_area(shape, a, b, c)
    return Geometry(volume, length, (length1, length2, length3), area)
end

# function geometry(shape::Plate, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, length, (length1, length2, length3), area, sil_area)
# end

calc_area(body::AbstractBody) = calc_area(shape(body), body)
calc_area(shape::Plate, a, b, c) = a * b * 2 + a * c * 2 + b * c * 2
function calc_area(shape::Plate, body)
    a = body.geometry.lengths[1] / 2
    b = body.geometry.lengths[2] / 2
    c = body.geometry.lengths[3] / 2
    a * b * 2 + a * c * 2 + b * c * 2
end

"""
    Cylinder <: Shape

A cylindrical organism.
"""
struct Cylinder{M,D,B} <: Shape
    mass::M
    density::D
    b::B
end

function geometry(shape::Cylinder, ::Naked)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    r = (volume / (shape.b * π * 2))^(1 / 3)
    length1 = shape.b * r * 2
    length2 = 2 * r
    area = calc_area(shape, r, length1)
    return Geometry(volume, length, (length1, length2), area)
end

# function geometry(shape::Cylinder, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, length, (length1, length2, length3), area, sil_area)
# end
calc_area(shape::Cylinder, r, l) = 2 * π * r * l + 2 * π * r^2
function calc_area(shape::Cylinder, body)
    r = body.geometry.lengths[2] / 2
    l = body.geometry.lengths[1]
    2 * π * r * l + 2 * π * r^2
end

"""
    sil_area_of_cylinder

Calculates the silhouette (projected) area of a cylinder.
Equation from Fig. 11.6 in Campbell, G. S., & Norman, J. M.
(1998). Environmental Biophysics. Springer.
"""
calc_silhouette_area(body::AbstractBody, θ) = calc_silhouette_area(shape(body), body, θ)
function calc_silhouette_area(shape::Cylinder, body, θ)
    r = body.geometry.lengths[2] / 2
    l = body.geometry.lengths[1]
    2 * r * l * sin(θ) + π * r^2 * cos(θ)
end

"""
    Ellipsoid <: Shape

An ellipsoidal organism.
"""
struct Ellipsoid{M,D,B,C} <: Shape
    mass::M
    density::D
    b::B
    c::C
end

function geometry(shape::Ellipsoid, ::Naked)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    a = ((3 / 4)* Unitful.ustrip(volume) / (π * shape.b * shape.c)) ^ (1 / 3)
    b = a * shape.b
    c = a * shape.c
    length1 = a * 2m
    length2 = b * 2m
    length3 = c * 2m
    area = calc_area(shape, a, b, c)
    return Geometry(volume, length, (length1, length2, length3), area)
end
# function geometry(shape::Ellipsoid, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, length, (length1, length2, length3), area)
# end
function calc_area(shape::Ellipsoid, a, b, c)
    e = ((a ^ 2 - c ^ 2) ^ 0.5 ) / a # eccentricity
    2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)
end
function calc_area(shape::Ellipsoid, body)
    a = body.geometry.lengths[1] / 2
    b = body.geometry.lengths[2] / 2
    c = body.geometry.lengths[3] / 2
    e = ((a ^ 2 - c ^ 2) ^ 0.5 ) / a # eccentricity
    2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)
end

"""
    sil_area_of_ellipsoid

Calculates the silhouette (projected) area of a prolate spheroid.
"""
function calc_silhouette_area(shape::Ellipsoid, a, b, c, θ)
    a2 = cos(90) ^ 2 * (cos(θ) ^ 2 / a ^ 2 + sin(θ) ^ 2 / b ^ 2) + sin(90) ^ 2 / c ^ 2
    twohh = 2 * cos(90) * sin(90) * cos(θ) * (1 / b ^ 2 - 1 / a ^ 2)
    b2 = sin(θ) ^ 2 / a ^ 2 + cos(θ) ^ 2 / b ^ 2
    θ2 = 0.5 * atan(twohh, a2 - b2)
    sps = sin(θ2)
    cps = cos(θ2)
    a3 = cps * (a2 * cps + twohh * sps) + b2 * sps * sps
    b3 = sps * (a2 * sps - twohh * cps) + b2 * cps * cps
    semax1 = 1 / sqrt(a3)
    semax2 = 1 / sqrt(b3)
    π * semax1  * semax2
end
function calc_silhouette_area(shape::Ellipsoid, body, θ)
    a = body.geometry.lengths[1] / 2
    b = body.geometry.lengths[2] / 2
    c = body.geometry.lengths[3] / 2
    calc_silhouette_area(shape, a, b, c, θ)
end

"""
    LeopardFrog <: Shape

An frog-shaped organism. Based on the leopard frog (Tracy 1976 Ecol. Monog.)
"""
struct LeopardFrog{M,D} <: Shape
    mass::M
    density::D
end

function geometry(shape::LeopardFrog, ::Naked)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    area = calc_area(shape)
    return Geometry(volume, length, nothing, area)
end
# function geometry(shape::Ellipsoid, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, length, (length1, length2, length3), area)
# end
calc_area(shape::LeopardFrog) = begin
    mass_g = Unitful.uconvert(u"g", shape.mass)
    Unitful.uconvert(u"m^2", (12.79 * Unitful.ustrip(mass_g) ^ 0.606)u"cm^2") # eq in Fig. 5
end

"""
    sil_area_of_leopardfrog

Calculates the silhouette (projected) area of a leopard frog.
"""
function calc_silhouette_area(shape::LeopardFrog, θ)
    area = calc_area(shape)
    pct = 1.38171e-6 * θ ^ 4 - 1.93335e-4 * θ ^ 3 + 4.75761e-3 * θ ^ 2 - 0.167912 * θ + 45.8228
    sil_area = pct * area / 100
end


"""
    DesertIguana <: Shape

An lizard-shaped organism. Based on the desert iguana (Porter and Tracy 1984)
"""
struct DesertIguana{M,D} <: Shape
    mass::M
    density::D
end

function geometry(shape::DesertIguana, ::Naked)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    area = calc_area(shape)
    return Geometry(volume, length, nothing, area)
end
# function geometry(shape::Ellipsoid, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, length, (length1, length2, length3), area)
# end

calc_area(shape::DesertIguana) = begin
    mass_g = Unitful.uconvert(u"g", shape.mass)
    area = Unitful.uconvert(u"m^2", (10.4713 * Unitful.ustrip(mass_g) ^ 0.688)u"cm^2")
end

