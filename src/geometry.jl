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
    length = volume^(1 / 3)
    a = (volume / (shape.b * shape.c))^(1 / 3)
    b = shape.b * a
    c = shape.c * a
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    total = calc_area(shape, a, b, c)
    skin = total
    convection = total
    fat = 0.0u"m"
    return Geometry(volume, length, (; length1, length2, length3, fat), (; total, skin, convection))
end

function geometry(shape::Plate, fur::Fur)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    a = (volume / (shape.b * shape.c))^(1 / 3)
    b = shape.b * a
    c = shape.c * a
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    total = calc_area(shape, a + fur.thickness, b + fur.thickness, c + fur.thickness)
    skin = calc_area(shape, a, b, c)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    fat = 0.0u"m"
    return Geometry(volume, length, (; length1, length2, length3, fat), (; total, skin, convection))
end

function geometry(shape::Plate, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    r_flesh = (flesh_volume / (shape.b * shape.c))^(1 / 3) / 2
    a = (volume / (shape.b * shape.c))^(1 / 3)
    b = shape.b * a
    c = shape.c * a
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    fat = a - r_flesh
    total = calc_area(shape, a, b, c)
    skin = total
    convection = total
    return Geometry(volume, length, (; length1, length2, length3, fat), (; total, skin, convection))
end

function geometry(shape::Plate, fur::Fur, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    r_flesh = (flesh_volume / (shape.b * shape.c))^(1 / 3) / 2
    a = (volume / (shape.b * shape.c))^(1 / 3)
    b = shape.b * a
    c = shape.c * a
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    fat = a - r_flesh
    total = calc_area(shape, a + fur.thickness, b + fur.thickness, c + fur.thickness)
    skin = calc_area(shape, a, b, c)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    return Geometry(volume, length, (; length1, length2, length3, fat), (; total, skin, convection))
end

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
    length = volume^(1 / 3)
    r = (volume / (shape.b * π * 2))^(1 / 3)
    length1 = shape.b * r * 2
    length2 = 2 * r
    total = calc_area(shape, r, length1)
    skin = total
    convection = total
    fat = 0.0u"m"
    return Geometry(volume, length, (; length1, length2, fat), (; total, skin, convection))
end

function geometry(shape::Cylinder, fur::Fur)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    r_skin = (volume / (shape.b * π * 2))^(1 / 3)
    r_fur = r_skin + fur.thickness
    length1 = shape.b * r_skin * 2
    length2 = 2 * r_fur
    total = calc_area(shape, r_fur, length1)
    skin = calc_area(shape, r_skin, length1)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    fat = 0.0u"m"
    return Geometry(volume, length, (; length1, length2, fat), (; total, skin, convection))
end

function geometry(shape::Cylinder, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    length = volume^(1 / 3)
    r_skin = (volume / (shape.b * π * 2))^(1 / 3)
    length1 = shape.b * r_skin * 2
    length2 = 2 * r_skin
    r_flesh = (flesh_volume / (π * length1))^(1 / 2)
    fat = r_skin - r_flesh
    total = calc_area(shape, r_skin, length1)
    skin = total
    convection = total
    return Geometry(volume, length, (; length1, length2, fat), (; total, skin, convection))
end

function geometry(shape::Cylinder, fur::Fur, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    length = volume^(1 / 3)
    r_skin = (volume / (shape.b * π * 2))^(1 / 3)
    r_fur = r_skin + fur.thickness
    length1 = shape.b * r_skin * 2
    length2 = 2 * r_skin
    r_flesh = (flesh_volume / (π * length1))^(1 / 2)
    fat = r_skin - r_flesh
    total = calc_area(shape, r_fur, length1)
    skin = calc_area(shape, r_skin, length1)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair    
    return Geometry(volume, length, (; length1, length2, fat), (; total, skin, convection))
end

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
    Sphere <: Shape

A spherical organism shape.
"""
struct Sphere{M,D} <: Shape
    mass::M
    density::D
end

function geometry(shape::Sphere, ::Naked)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    r = ((3 / 4) * volume / π) ^ (1 / 3)
    length1 = r * 2u"m"
    total = calc_area(shape, r)u"m^2"
    skin = total
    convection = total
    fat = 0.0u"m"
    return Geometry(volume, length, (; length1, fat), (; total, skin, convection))
end

function geometry(shape::Sphere, fur::Fur)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    r_skin = ((3 / 4)* volume / π) ^ (1 / 3)
    r_fur = r_skin + fur.thickness
    length1 = r_skin * 2u"m"
    total = calc_area(shape, r_fur)u"m^2"
    skin = calc_area(shape, r_skin)u"m^2"
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    fat = 0.0u"m"
    return Geometry(volume, length, (; length1, fat), (; total, skin, convection))
end

function geometry(shape::Sphere, fat::Fat)
    volume = shape.mass / shape.density
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    length = volume^(1 / 3)
    r = ((3 / 4) * volume / π) ^ (1 / 3)
    length1 = r * 2u"m"
    r_flesh = ((3 * flesh_volume) / (4 * π)) ^ (1 / 3)
    fat = r - r_flesh
    total = calc_area(shape, r)u"m^2"
    skin = total
    convection = total
    return Geometry(volume, length, (length1, fat), (; total, skin, convection))
end

function geometry(shape::Sphere, fur::Fur, fat::Fat)
    volume = shape.mass / shape.density
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume    
    length = volume^(1 / 3)
    r_skin = ((3 / 4)* volume / π) ^ (1 / 3)
    r_fur = r_skin + fur.thickness
    length1 = r_skin * 2u"m"
    r_flesh = ((3 * flesh_volume) / (4 * π)) ^ (1 / 3)
    fat = r_skin - r_flesh    
    total = calc_area(shape, r_fur)u"m^2"
    skin = calc_area(shape, r_skin)u"m^2"
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    return Geometry(volume, length, (; length1, fat), (; total, skin, convection))
end

function calc_area_hair(fibre_diameter, fibre_density, skin)
    π * (fibre_diameter / 2) ^ 2 * (fibre_density * skin)
end
function calc_area(shape::Sphere, body::Body)
    r = body.geometry.lengths[1] / 2
    return calc_area(shape, r)
end
function calc_area(shape::Sphere, r)
    4 * π * r ^ 2
end

function calc_silhouette_area(shape::Sphere, body)
    r = body.geometry.lengths[1] / 2
    return calc_silhouette_area(shape, r)
end
calc_silhouette_area(shape::Sphere, r) = max(π * r ^ 2)

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
    length = volume^(1 / 3)
    b = ((3 / 4) * volume / (π * shape.b)) ^ (1 / 3)
    c = b
    a = b * shape.b
    e = ((ustrip(u"m", a) ^ 2 - ustrip(u"m", c) ^ 2) ^ (1 / 2)) / ustrip(u"m", a)
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    total = calc_area(shape, ustrip(u"m", a), ustrip(u"m", b), ustrip(u"m", c), e)
    skin = total
    convection = total
    fat = 0.0u"m"
    return Geometry(volume, length, (; length1, length2, length3, fat), (; total, skin, convection))
end

function geometry(shape::Ellipsoid, fur::Fur)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    b = ((3 / 4) * volume / (π * shape.b)) ^ (1 / 3)
    c = b
    a = b * shape.b
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    e = ((ustrip(u"m", a) ^ 2 - ustrip(u"m", c) ^ 2) ^ (1 / 2)) / ustrip(u"m", a)
    a_fur = a + fur.thickness
    b_fur = b + fur.thickness
    c_fur = c + fur.thickness
    e = ((a ^ 2 - c ^ 2) ^ (1 / 2)) / a
    e_fur = ((a_fur ^ 2 - c_fur ^ 2) ^ (1 / 2)) / a_fur
    total = calc_area(shape, ustrip(u"m", a_fur),  ustrip(u"m", b_fur),  ustrip(u"m", c_fur),  e_fur)
    skin = calc_area(shape, ustrip(u"m", a), ustrip(u"m", b), ustrip(u"m", c), e)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair 
    fat = 0.0u"m"
    return Geometry(volume, length, (; length1, length2, length3, fat), (; total, skin, convection))
end

function prolate_fat_layer(
    flesh_volume,
    fat_volume,
    shape_b,
    semi_minor_flesh
)
    # Flesh is approximated as a prolate spheroid:
    # Volume = 4/3 π * A * B * C
    # C = B, A = shape_b * B
    # → B = ((3 * volume) / (4 * shape_b * π))^(1/3)
    # Fat thickness X is root of cubic: A X^3 + B X^2 + C X + D = 0
    A = 1.0
    B = shape_b * semi_minor_flesh + 2 * semi_minor_flesh
    C = 2 * shape_b * semi_minor_flesh^2 + semi_minor_flesh^2
    D = shape_b * semi_minor_flesh^3 - (( (fat_volume + flesh_volume) * 3.0 ) / (4.0 * π))

    # Components of cubic formula

    T1a = (-B)^3 / (27 * A^3)
    T1b = (B * C) / (6 * A^2)
    T1c = D / (2 * A)
    T1  = T1a + T1b - T1c

    T2a = T1^2
    T2b = ( (C / (3*A)) - (B^2) / (9 * A^2) )^3

    # Prevent sqrt of negative number
    T2 = (T2a + T2b >= 0) ? sqrt(T2a + T2b) : 0.0

    T3 = B / (3*A)

    # Cube roots with sign handling

    function signed_cuberoot(x)
        x < 0 ? -((-x)^(1/3)) : x^(1/3)
    end

    root1 = signed_cuberoot(T1 + T2)
    root2 = signed_cuberoot(T1 - T2)

    fat = (root1 + root2 - T3) * u"m"

    # If negative, not enough fat to cover spheroid
    return max(0.0u"m", fat)
end

function geometry(shape::Ellipsoid, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    b_flesh = (((3 / 4) * flesh_volume) / (π * shape.b)) ^ (1 / 3)
    c_flesh = b_flesh # assuming c = b
    a_flesh = shape.b * b_flesh
    fat = prolate_fat_layer(
        ustrip(u"m^3", flesh_volume), 
        ustrip(u"m^3", fat_volume), 
        shape.b, 
        ustrip(u"m", b_flesh))
    length = volume^(1 / 3)
    a = a_flesh + fat
    b = b_flesh + fat
    c = c_flesh + fat
    e = ((a ^ 2 - c ^ 2) ^ (1 / 2)) / a
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    total = calc_area(shape, ustrip(u"m", a), ustrip(u"m", b), ustrip(u"m", c), e)
    skin = total
    convection = total
    return Geometry(volume, length, (; length1, length2, length3, fat), (; total, skin, convection))
end

function geometry(shape::Ellipsoid, fur::Fur, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    b_flesh = ((3.0 * flesh_volume) / (4.0 * shape.b * π))^(1/3)
    c_flesh = b_flesh # assuming c = b
    a_flesh = shape.b * b_flesh
    fat = prolate_fat_layer(
        ustrip(u"m^3", flesh_volume), 
        ustrip(u"m^3", fat_volume), 
        shape.b, 
        ustrip(u"m", b_flesh))
    length = volume^(1 / 3)
    a = a_flesh + fat
    b = b_flesh + fat
    c = c_flesh + fat
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    a_fur = a + fur.thickness
    b_fur = b + fur.thickness
    c_fur = c + fur.thickness
    e = ((a ^ 2 - c ^ 2) ^ (1 / 2)) / a
    e_fur = ((a_fur ^ 2 - c_fur ^ 2) ^ (1 / 2)) / a_fur
    total = calc_area(shape, ustrip(u"m", a_fur),  ustrip(u"m", b_fur),  ustrip(u"m", c_fur),  e_fur)
    skin = calc_area(shape, ustrip(u"m", a), ustrip(u"m", b), ustrip(u"m", c), e)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair 
    return Geometry(volume, length, (; length1, length2, length3, fat), (; total, skin, convection))
end

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

function calc_area(shape::Ellipsoid, a, b, c, e)
    (2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)) * u"m^2"
end

# TODO: should we use this rather than the one below that ignores theta?
"""
    calc_silhouette_area

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
    return calc_silhouette_area(shape, a, b, c, θ)
end
#calc_silhouette_area(shape::Ellipsoid, a, b, c, θ) = max(π * a * c, π * b * c)

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
    length = volume^(1 / 3)
    length1 = (volume / (4 * π)) ^ (1 / 3)
    total = calc_area(shape, body)
    skin = total
    convection = total
    fat = 0.0u"m"
    return Geometry(volume, length, (; length1, fat), (; total, skin, convection))
end

calc_area(shape::LeopardFrog) = calc_area(shape)
function calc_area(shape::LeopardFrog)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    Unitful.uconvert(u"m^2", (12.79 * Unitful.ustrip(mass_g) ^ 0.606)u"cm^2") # eq in Fig. 5
end

calc_silhouette_area(shape::LeopardFrog, θ) = calc_silhouette_area(shape::LeopardFrog, θ)
function calc_silhouette_area(shape::LeopardFrog, θ)
    total = calc_area(shape)
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
    length = volume^(1 / 3)
    length1 = (volume / (4 * π)) ^ (1 / 3)
    total = calc_area(shape).total
    skin = total
    convection = total
    fat = 0.0u"m"
    return Geometry(volume, length, (; length1, fat), (; total, skin, convection))
end

function calc_area(shape::DesertIguana)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    total = Unitful.uconvert(u"m^2", (10.4713 * Unitful.ustrip(mass_g) ^ 0.688)u"cm^2")
    ventral = Unitful.uconvert(u"m^2", (0.425 * Unitful.ustrip(mass_g) ^ 0.85)u"cm^2")
    return (; total, ventral)
end

function calc_silhouette_area(shape::DesertIguana)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    normal = Unitful.uconvert(u"m^2", (3.798 * Unitful.ustrip(mass_g) ^ 0.683)u"cm^2")
    parallel = Unitful.uconvert(u"m^2", (0.694 * Unitful.ustrip(mass_g) ^ 0.743)u"cm^2")
    return (; normal, parallel)
end

"""
    bird_skin_area

Skin area for a bird based on Walsberg & King (1978)
    Journal of Experimental Biology, 76, 185-189.
Note that skin area is greater than plumage area.
"""
function bird_skin_area(shape)
    return Unitful.uconvert(u"m^2", (10.0 * Unitful.ustrip(u"g", shape.mass) ^ 0.667)u"cm^2")
end

"""
    bird_plumage_area

Plumage area for a bird based on Walsberg & King (1978)
    Journal of Experimental Biology, 76, 185-189.
Note that plumage area is less than skin area.
"""
function bird_plumage_area(shape)
    return Unitful.uconvert(u"m^2", (8.11 * Unitful.ustrip(u"g", shape.mass) ^ 0.667)u"cm^2")
end

"""
    mammal_skin_area

Skin area for a bird based on Walsberg & King (1978)
    Journal of Experimental Biology, 76, 185-189.
"""
function mammal_skin_area(shape)
    return Unitful.uconvert(u"m^2", (1110.0 * Unitful.ustrip(u"kg", shape.mass) ^ 0.65)u"cm^2")
end

"""
    mammal_fur_area

Fur area for a mammal based on Stahl (1967)
    Journal of Applied Physiology, 22, 453-460.
"""
function mammal_fur_area(body)
    fur_multiplier = body.geometry.area.total / body.geometry.area.skin
    return Unitful.uconvert(u"m^2", (1110.0 * Unitful.ustrip(u"kg", body.shape.mass) ^ 0.65 * fur_multiplier)u"cm^2")
end
