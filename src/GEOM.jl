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
    area = area_of_cylinder(r, length1)
    return Geometry(volume, length, (length1, length2), area)
end

# function geometry(shape::Cylinder, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, length, (length1, length2, length3), area, sil_area)
# end

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
    a = ((3 / 4)* volume / (π * shape.b * shape.c)) ^ (1 / 3)
    b = a * shape.b
    c = a * shape.c
    p = 1.6075
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    area = area_of_ellipsoid(a, b, c)
    return Geometry(volume, length, (length1, length2, length3), area)
end
# function geometry(shape::Ellipsoid, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, length, (length1, length2, length3), area)
# end

function area_of_ellipsoid(a, b, c)
    e = ((a ^ 2 - c ^ 2) ^ 0.5 ) / a # eccentricity
    2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)
end

"""
    sil_area_of_ellipsoid

Calculates the silhouette (projected) area of a prolate spheroid.
"""
function sil_area_of_ellipsoid(a, b, c, θ)
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

function area_of_cylinder(r, l)
    2 * π * r * l + 2 * π * r^2
end

"""
    sil_area_of_cylinder

Calculates the silhouette (projected) area of a cylinder.
Equation from Fig. 11.6 in Campbell, G. S., & Norman, J. M.
(1998). Environmental Biophysics. Springer.
"""
function sil_area_of_cylinder(r, l, θ)
    2 * r * l * sin(θ) + π * r^2 * cos(θ)
end
