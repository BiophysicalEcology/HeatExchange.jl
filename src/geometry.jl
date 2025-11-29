shape(body::AbstractBody) = body.shape
insulation(body::AbstractBody) = body.insulation
geometry(body::AbstractBody) = body.geometry
surface_area(body::AbstractBody) = surface_area(shape(body), body)

# functions to extract appropriate surface areas from different objects

get_total_area(body::AbstractBody) = get_total_area(shape(body), insulation(body), body)
get_skin_area(body::AbstractBody) = get_skin_area(shape(body), insulation(body), body)
get_evaporation_area(body::AbstractBody) = get_evaporation_area(shape(body), insulation(body), body)

# for composite insulation cases (fat and fur/feathers)
outer_insulation(ins::Insulation) = ins
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
get_total_area(shape, ins::CompositeInsulation, body) =
    get_total_area(shape, outer_insulation(ins), body)
get_skin_area(shape, ins::CompositeInsulation, body) =
    get_skin_area(shape, outer_insulation(ins), body)
get_evaporation_area(shape, ins::CompositeInsulation, body) =
    get_evaporation_area(shape, outer_insulation(ins), body)

# functions to get the appropriate radii
get_r_skin(body::AbstractBody) = get_r_skin(shape(body), insulation(body), body)
get_r_insulation(body::AbstractBody) = get_r_insulation(shape(body), insulation(body), body)
get_r_flesh(body::AbstractBody) = get_r_flesh(shape(body), insulation(body), body)

"""
    silhouette_area(body::AbstractBody, θ)

Calculates the silhouette (projected) area of a cylinder given a solar zenith angle, θ.
Calculates the silhouette (projected) area of a cylinder.
Equation from Fig. 11.6 in Campbell, G. S., & Norman, J. M.
(1998). Environmental Biophysics. Springer.
"""
silhouette_area(body::AbstractBody, θ) = silhouette_area(shape(body), insulation(body), body, θ)

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
    length = (volume / (shape.b * shape.c))^(1 / 3)
    width = shape.b * length
    height = shape.c * length
    total = surface_area(shape, length, width, height)
    characteristic_dimension = width * 2 #volume^(1 / 3) 
    return Geometry(volume, characteristic_dimension, (; length, width, height), (; total))
end

# area functions
get_total_area(shape::Plate, insulation::Naked, body) = body.geometry.area.total
get_skin_area(shape::Plate, insulation::Naked, body) = body.geometry.area.total
get_evaporation_area(shape::Plate, insulation::Naked, body) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::Plate, insulation::Naked, body) = body.geometry.length.width / 2
get_r_insulation(shape::Plate, insulation::Naked, body) = body.geometry.length.width / 2
get_r_flesh(shape::Plate, insulation::Naked, body) = body.geometry.length.width / 2

function geometry(shape::Plate, fur::Fur)
    volume = shape.mass / shape.density
    length = (volume / (shape.b * shape.c))^(1 / 3)
    width = shape.b * length
    height = shape.c * length
    length_fur = length + fur.thickness * 2
    width_fur = width + fur.thickness * 2
    height_fur = height + fur.thickness * 2
    total = surface_area(shape, length_fur, width_fur, height_fur)
    skin = surface_area(shape, length, width, height)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    fat = 0.0u"m"
    characteristic_dimension = width_fur * 2 #volume^(1 / 3) 
    return Geometry(volume, characteristic_dimension, (; length, width, height, length_fur, width_fur, height_fur, fat), (; total, skin, convection))
end

# area functions
get_total_area(shape::Plate, insulation::Fur, body) = body.geometry.area.total
get_skin_area(shape::Plate, insulation::Fur, body) = body.geometry.area.skin
get_evaporation_area(shape::Plate, insulation::Fur, body) = body.geometry.area.convection

# "radii" functions
get_r_skin(shape::Plate, insulation::Fur, body) = body.geometry.length.width / 2
get_r_insulation(shape::Plate, insulation::Fur, body) = body.geometry.length.width_fur / 2
get_r_flesh(shape::Plate, insulation::Fur, body) = body.geometry.length.width / 2

function geometry(shape::Plate, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    volume = shape.mass / shape.density
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    r_flesh = (flesh_volume / (shape.b * shape.c))^(1 / 3) / 2
    length = (volume / (shape.b * shape.c))^(1 / 3)
    width = shape.b * length
    height = shape.c * length
    fat = length - r_flesh
    total = surface_area(shape, length, width, height)
    characteristic_dimension = width * 2 #volume^(1 / 3) 
    return Geometry(volume, characteristic_dimension, (; length, width, height, fat), (; total))
end

# area functions
get_total_area(shape::Plate, insulation::Fat, body) = body.geometry.area.total
get_skin_area(shape::Plate, insulation::Fat, body) = body.geometry.area.total
get_evaporation_area(shape::Plate, insulation::Fat, body) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::Plate, insulation::Fat, body) = body.geometry.length.width / 2
get_r_insulation(shape::Plate, insulation::Fat, body) = body.geometry.length.width / 2
get_r_flesh(shape::Plate, insulation::Fat, body) = body.geometry.length.width / 2 - body.geometry.length.fat

function geometry(shape::Plate, fur::Fur, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    volume = shape.mass / shape.density
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    r_flesh = (flesh_volume / (shape.b * shape.c))^(1 / 3) / 2
    length = (volume / (shape.b * shape.c))^(1 / 3)
    width = shape.b * length
    height = shape.c * length
    length_fur = length + fur.thickness * 2
    width_fur = width + fur.thickness * 2
    height_fur = height + fur.thickness * 2
    fat = length - r_flesh
    total = surface_area(shape, length_fur, width_fur, height_fur)
    skin = surface_area(shape, length, width, height)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = width_fur * 2 #volume^(1 / 3) 
    return Geometry(volume, characteristic_dimension, (; length, width, height, length_fur, width_fur, height_fur, fat), (; total, skin, convection))
end

# surface area functions
get_total_area(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.area.total
get_skin_area(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.area.skin
get_evaporation_area(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.area.convection

# "radii" functions
get_r_skin(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width / 2
get_r_insulation(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_fur / 2
get_r_flesh(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width / 2 - body.geometry.length.fat

function surface_area(shape::Plate, body)
    length = body.geometry.length.length
    width = body.geometry.length.width
    height = body.geometry.length.height
    surface_area(shape, length, width, height)
end
surface_area(shape::Plate, length, width, height) = length * width * 2 + length * height * 2 + width * height * 2

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
    radius = (volume / (shape.b * π * 2))^(1 / 3)
    length = shape.b * radius * 2
    total = surface_area(shape, radius, length)
    characteristic_dimension = radius * 2 #volume^(1 / 3) 
    return Geometry(volume, characteristic_dimension, (; length, radius), (; total))
end

# area functions
get_total_area(shape::Cylinder, insulation::Naked, body) = body.geometry.area.total
get_skin_area(shape::Cylinder, insulation::Naked, body) = body.geometry.area.total
get_evaporation_area(shape::Cylinder, insulation::Naked, body) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::Cylinder, insulation::Naked, body) = body.geometry.length.radius
get_r_insulation(shape::Cylinder, insulation::Naked, body) = body.geometry.length.radius
get_r_flesh(shape::Cylinder, insulation::Naked, body) = body.geometry.length.radius

function geometry(shape::Cylinder, fur::Fur)
    volume = shape.mass / shape.density
    radius_skin = (volume / (shape.b * π * 2))^(1 / 3)
    radius_fur = radius_skin + fur.thickness
    length_skin = shape.b * radius_skin * 2
    length_fur = 2 * radius_fur
    total = surface_area(shape, radius_fur, length_fur)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = radius_fur * 2 #volume^(1 / 3) 
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur, length_skin, length_fur), (; total, skin, convection))
end

# surface area functions
get_total_area(shape::Cylinder, insulation::Fur, body) = body.geometry.area.total
get_skin_area(shape::Cylinder, insulation::Fur, body) = body.geometry.area.skin
get_evaporation_area(shape::Cylinder, insulation::Fur, body) = body.geometry.area.convection

# "radii" functions
get_r_skin(shape::Cylinder, insulation::Fur, body) = body.geometry.length.radius_skin
get_r_insulation(shape::Cylinder, insulation::Fur, body) = body.geometry.length.radius_fur
get_r_flesh(shape::Cylinder, insulation::Fur, body) = body.geometry.length.radius_skin

function geometry(shape::Cylinder, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    radius = (volume / (shape.b * π * 2))^(1 / 3)
    length = shape.b * radius * 2
    radius_flesh = (flesh_volume / (π * length))^(1 / 2)
    fat = radius - radius_flesh
    total = surface_area(shape, radius, length)
    characteristic_dimension = radius * 2 #volume^(1 / 3) 
    return Geometry(volume, characteristic_dimension, (; radius, length, fat), (; total))
end

# surface area functions
get_total_area(shape::Cylinder, insulation::Fat, body) = body.geometry.area.total
get_skin_area(shape::Cylinder, insulation::Fat, body) = body.geometry.area.total
get_evaporation_area(shape::Cylinder, insulation::Fat, body) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::Cylinder, insulation::Fat, body) = body.geometry.length.radius
get_r_insulation(shape::Cylinder, insulation::Fat, body) = body.geometry.length.radius
get_r_flesh(shape::Cylinder, insulation::Fat, body) = body.geometry.length.radius - body.geometry.length.fat

function geometry(shape::Cylinder, fur::Fur, fat::Fat)
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    volume = shape.mass / shape.density
    flesh_volume = volume - fat_volume
    radius_skin = (volume / (shape.b * π * 2))^(1 / 3)
    radius_fur = radius_skin + fur.thickness
    length_skin = shape.b * radius_skin * 2
    length_fur = 2 * radius_fur
    radius_flesh = (flesh_volume / (π * length_skin))^(1 / 2)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_fur, length_fur)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = radius_fur * 2 #volume^(1 / 3) 
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur, length_skin, length_fur, fat), (; total, skin, convection))
end

# surface area functions
get_total_area(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.area.total
get_skin_area(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.area.skin
get_evaporation_area(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.area.convection

# "radii" functions
get_r_skin(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin
get_r_insulation(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.length.radius_fur
get_r_flesh(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat

# function surface_area(shape::Cylinder, body::AbstractBody)
#     r = body.geometry.length.radius
#     l = body.geometry.length.length
#     surface_area(shape, r, l)
# end
surface_area(shape::Cylinder, r, l) = 2 * π * r * l + 2 * π * r^2

function silhouette_area(shape::Cylinder, ::Naked, body::AbstractBody, θ)
    r = body.geometry.length.radius
    l = body.geometry.length.length
    return silhouette_area(shape, r, l, θ)
end
function silhouette_area(shape::Cylinder, insulation::Fat, body::AbstractBody, θ)
    r = body.geometry.length.radius
    l = body.geometry.length.length
    return silhouette_area(shape, r, l, θ)
end
function silhouette_area(shape::Cylinder, insulation::Fur, body::AbstractBody, θ)
    r = body.geometry.length.radius_fur
    l = body.geometry.length.length_fur
    return silhouette_area(shape, r, l, θ)
end
function silhouette_area(shape::Cylinder, insulation::CompositeInsulation, body::AbstractBody, θ)
    r = body.geometry.length.radius_fur
    l = body.geometry.length.length_fur
    return silhouette_area(shape, r, l, θ)
end
silhouette_area(shape::Cylinder, r, l, θ) = 2 * r * l * sin(θ) + π * r^2 * cos(θ)

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
    radius = ((3 / 4) * volume / π) ^ (1 / 3)
    total = surface_area(shape, radius)
    characteristic_dimension = radius * 2 #volume^(1 / 3)
    return Geometry(volume, characteristic_dimension, (; radius), (; total))
end

# area functions
get_total_area(shape::Sphere, insulation::Naked, body) = body.geometry.area.total
get_skin_area(shape::Sphere, insulation::Naked, body) = body.geometry.area.total
get_evaporation_area(shape::Sphere, insulation::Naked, body) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::Sphere, insulation::Naked, body) = body.geometry.length.radius
get_r_insulation(shape::Sphere, insulation::Naked, body) = body.geometry.length.radius
get_r_flesh(shape::Sphere, insulation::Naked, body) = body.geometry.length.radius

function geometry(shape::Sphere, fur::Fur)
    volume = shape.mass / shape.density
    radius_skin = ((3 / 4)* volume / π) ^ (1 / 3)
    radius_fur = radius_skin + fur.thickness
    total = surface_area(shape, radius_fur)u"m^2"
    skin = surface_area(shape, radius_skin)u"m^2"
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = radius_fur * 2 #volume^(1 / 3)
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur), (; total, skin, convection))
end

# surface area functions
get_total_area(shape::Sphere, insulation::Fur, body) = body.geometry.area.total
get_skin_area(shape::Sphere, insulation::Fur, body) = body.geometry.area.skin
get_evaporation_area(shape::Sphere, insulation::Fur, body) = body.geometry.area.convection

# "radii" functions
get_r_skin(shape::Sphere, insulation::Fur, body) = body.geometry.length.radius_skin
get_r_insulation(shape::Sphere, insulation::Fur, body) = body.geometry.length.radius_fur
get_r_flesh(shape::Sphere, insulation::Fur, body) = body.geometry.length.radius_skin

function geometry(shape::Sphere, fat::Fat)
    volume = shape.mass / shape.density
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume
    radius = ((3 / 4) * volume / π) ^ (1 / 3)
    radius_flesh = ((3 * flesh_volume) / (4 * π)) ^ (1 / 3)
    fat = radius - radius_flesh
    total = surface_area(shape, radius)
    characteristic_dimension = radius * 2 #volume^(1 / 3)
    return Geometry(volume, characteristic_dimension, (; radius, fat), (; total))
end

# surface area functions
get_total_area(shape::Sphere, insulation::Fat, body) = body.geometry.area.total
get_skin_area(shape::Sphere, insulation::Fat, body) = body.geometry.area.total
get_evaporation_area(shape::Sphere, insulation::Fat, body) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::Sphere, insulation::Fat, body) = body.geometry.length.radius
get_r_insulation(shape::Sphere, insulation::Fat, body) = body.geometry.length.radius
get_r_flesh(shape::Sphere, insulation::Fat, body) = body.geometry.length.radius - body.geometry.length.fat

function geometry(shape::Sphere, fur::Fur, fat::Fat)
    volume = shape.mass / shape.density
    fat_mass = shape.mass * fat.fraction
    fat_volume = fat_mass / fat.density
    flesh_volume = volume - fat_volume    
    radius_skin = ((3 / 4) * volume / π) ^ (1 / 3)
    radius_fur = radius_skin + fur.thickness
    radius_flesh = ((3 * flesh_volume) / (4 * π)) ^ (1 / 3)
    fat = radius_skin - radius_flesh
    total = surface_area(shape, radius_fur)u"m^2"
    skin = surface_area(shape, radius_skin)u"m^2"
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = radius_fur * 2 #volume^(1 / 3)
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur, fat), (; total, skin, convection))
end

# surface area functions
get_total_area(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.area.total
get_skin_area(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.area.skin
get_evaporation_area(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.area.convection

# "radii" functions
get_r_skin(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin
get_r_insulation(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_fur
get_r_flesh(shape::Sphere, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat

function calc_area_hair(fibre_diameter, fibre_density, skin)
    π * (fibre_diameter / 2) ^ 2 * (fibre_density * skin)
end
function surface_area(shape::Sphere, body::AbstractBody)
    r = body.geometry.length[1] / 2
    return surface_area(shape, r)
end
function surface_area(shape::Sphere, r)
    4 * π * r ^ 2
end

function silhouette_area(shape::Sphere, body::AbstractBody)
    r = body.geometry.length[1] / 2
    return silhouette_area(shape, r)
end
silhouette_area(shape::Sphere, r) = π * r ^ 2

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
    b_semi_minor = ((3 / 4) * volume / (π * shape.b)) ^ (1 / 3)
    c_semi_minor = b_semi_minor
    a_semi_major = b_semi_minor * shape.b
    e = ((ustrip(u"m", a_semi_major) ^ 2 - ustrip(u"m", c_semi_minor) ^ 2) ^ (1 / 2)) / ustrip(u"m", a_semi_major)
    total = surface_area(shape, ustrip(u"m", a_semi_major), ustrip(u"m", b_semi_minor), ustrip(u"m", b_semi_minor), e)
    characteristic_dimension = b_semi_minor * 2 #volume^(1 / 3)   
    return Geometry(volume, characteristic_dimension, (; a_semi_major, b_semi_minor, c_semi_minor), (; total))
end

# area functions
get_total_area(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.area.total
get_skin_area(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.area.total
get_evaporation_area(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.length.b_semi_minor
get_r_insulation(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.length.b_semi_minor
get_r_flesh(shape::Ellipsoid, insulation::Naked, body::AbstractBody) = body.geometry.length.b_semi_minor

function geometry(shape::Ellipsoid, fur::Fur)
    volume = shape.mass / shape.density
    b_semi_minor = ((3 / 4) * volume / (π * shape.b)) ^ (1 / 3)
    c_semi_minor = b_semi_minor
    a_semi_major = b_semi_minor * shape.b
    e = ((ustrip(u"m", a_semi_major) ^ 2 - ustrip(u"m", c_semi_minor) ^ 2) ^ (1 / 2)) / ustrip(u"m", a_semi_major)
    a_semi_major_fur = a_semi_major + fur.thickness
    b_semi_minor_fur = b_semi_minor + fur.thickness
    c_semi_minor_fur = c_semi_minor + fur.thickness
    e = ((a_semi_major ^ 2 - c_semi_minor ^ 2) ^ (1 / 2)) / a_semi_major
    e_fur = ((a_semi_major_fur ^ 2 - c_semi_minor_fur ^ 2) ^ (1 / 2)) / a_semi_major_fur
    total = surface_area(shape, ustrip(u"m", a_semi_major_fur),  ustrip(u"m", b_semi_minor_fur),  ustrip(u"m", c_semi_minor_fur),  e_fur)
    skin = surface_area(shape, ustrip(u"m", a_semi_major), ustrip(u"m", b_semi_minor), ustrip(u"m", c_semi_minor), e)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair 
    characteristic_dimension = b_semi_minor_fur * 2 #volume^(1 / 3)
    return Geometry(volume, characteristic_dimension, (; a_semi_major, b_semi_minor, c_semi_minor, a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur), (; total, skin, convection))
end

# surface area functions
get_total_area(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.area.total
get_skin_area(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.area.skin
get_evaporation_area(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.area.convection

# "radii" functions
get_r_skin(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.length.b_semi_minor
get_r_insulation(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.length.b_semi_minor_fur
get_r_flesh(shape::Ellipsoid, insulation::Fur, body::AbstractBody) = body.geometry.length.b_semi_minor

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
    characteristic_dimension = volume^(1 / 3)
    a_semi_major = a_flesh + fat
    b_semi_minor = b_flesh + fat
    c_semi_minor = c_flesh + fat
    e = ((a_semi_major ^ 2 - c_semi_minor ^ 2) ^ (1 / 2)) / a_semi_major
    total = surface_area(shape, ustrip(u"m", a_semi_major), ustrip(u"m", b_semi_minor), ustrip(u"m", c_semi_minor), e)
    characteristic_dimension = b_semi_minor * 2 #volume^(1 / 3)
    return Geometry(volume, characteristic_dimension, (; a_semi_major, b_semi_minor, c_semi_minor, fat), (; total))
end

# surface area functions
get_total_area(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.area.total
get_skin_area(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.area.total
get_evaporation_area(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.length.b_semi_minor
get_r_insulation(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.length.b_semi_minor
get_r_flesh(shape::Ellipsoid, insulation::Fat, body::AbstractBody) = body.geometry.length.b_semi_minor - body.geometry.length.fat

function geometry(shape::Ellipsoid, fur::Fur, fat::Fat)
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
    a_semi_major = a_flesh + fat
    b_semi_minor = b_flesh + fat
    c_semi_minor = c_flesh + fat
    e = ((a_semi_major ^ 2 - c_semi_minor ^ 2) ^ (1 / 2)) / a_semi_major
    a_semi_major_fur = a_semi_major + fur.thickness
    b_semi_minor_fur = b_semi_minor + fur.thickness
    c_semi_minor_fur = c_semi_minor + fur.thickness
    e = ((a_semi_major ^ 2 - c_semi_minor ^ 2) ^ (1 / 2)) / a_semi_major
    e_fur = ((a_semi_major_fur ^ 2 - c_semi_minor_fur ^ 2) ^ (1 / 2)) / a_semi_major_fur
    total = surface_area(shape, ustrip(u"m", a_semi_major_fur),  ustrip(u"m", b_semi_minor_fur),  ustrip(u"m", c_semi_minor_fur),  e_fur)
    skin = surface_area(shape, ustrip(u"m", a_semi_major), ustrip(u"m", b_semi_minor), ustrip(u"m", c_semi_minor), e)
    area_hair = calc_area_hair(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair 
    characteristic_dimension = b_semi_minor_fur * 2 #volume^(1 / 3)
    return Geometry(volume, characteristic_dimension, (; a_semi_major, b_semi_minor, c_semi_minor, a_semi_major_fur, b_semi_minor_fur, c_semi_minor_fur, fat), (; total, skin, convection))
end

# surface area functions
get_total_area(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.area.total
get_skin_area(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.area.skin
get_evaporation_area(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.area.convection

# "radii" functions
get_r_skin(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.length.b_semi_minor
get_r_insulation(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.length.b_semi_minor_fur
get_r_flesh(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody) = body.geometry.length.b_semi_minor - body.geometry.length.fat

function surface_area(shape::Ellipsoid, body::AbstractBody)
    a = body.geometry.length.a_semi_major
    b = body.geometry.length.b_semi_minor
    c = body.geometry.length.c_semi_minor
    return surface_area(shape, a, b, c)
end
function surface_area(shape::Ellipsoid, a, b, c)
    #e = ((a ^ 2 - c ^ 2) ^ 0.5 ) / a # eccentricity
    #2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)
    p =  1.6075
    return(4 * π * (((a ^ p * b ^ p + a ^ p * c ^ p + b ^ p * c ^ p)) / 3) ^ (1 / p))
end

function surface_area(shape::Ellipsoid, a, b, c, e)
    (2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)) * u"m^2"
end

# TODO: should we use this rather than the one below that ignores theta?
"""
    silhouette_area

Calculates the silhouette (projected) area of a prolate spheroid.
"""
function silhouette_area(shape::Ellipsoid, a, b, c, θ)
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

function silhouette_area(shape::Ellipsoid, ::Naked, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major
    b = body.geometry.length.b_semi_minor
    c = body.geometry.length.c_semi_minor
    return silhouette_area(shape, a, b, c, θ)
end

function silhouette_area(shape::Ellipsoid, insulation::Fur, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major_fur
    b = body.geometry.length.b_semi_minor_fur
    c = body.geometry.length.c_semi_minor_fur
    return silhouette_area(shape, a, b, c, θ)
end

function silhouette_area(shape::Ellipsoid, insulation::Fat, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major
    b = body.geometry.length.b_semi_minor
    c = body.geometry.length.c_semi_minor
    return silhouette_area(shape, a, b, c, θ)
end

function silhouette_area(shape::Ellipsoid, insulation::CompositeInsulation, body::AbstractBody, θ)
    a = body.geometry.length.a_semi_major_fur
    b = body.geometry.length.b_semi_minor_fur
    c = body.geometry.length.c_semi_minor_fur
    return silhouette_area(shape, a, b, c, θ)
end
#silhouette_area(shape::Ellipsoid, a, b, c, θ) = max(π * a * c, π * b * c)

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
    r_skin = (volume / (4 * π)) ^ (1 / 3)
    total = surface_area(shape)
    return Geometry(volume, characteristic_dimension, (; r_skin), (; total))
end

# surface area functions
get_total_area(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.area.total
get_skin_area(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.area.total
get_evaporation_area(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin
get_r_insulation(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin
get_r_flesh(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin

function surface_area(shape::LeopardFrog)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    Unitful.uconvert(u"m^2", (12.79 * Unitful.ustrip(mass_g) ^ 0.606)u"cm^2") # eq in Fig. 5
end

function silhouette_area(shape::LeopardFrog, θ)
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
    r_skin = (volume / (4 * π)) ^ (1 / 3)
    total, ventral = surface_area(shape)
    return Geometry(volume, characteristic_dimension, (; r_skin), (; total, ventral))
end

# surface area functions
get_total_area(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.area.total
get_skin_area(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.area.total
get_evaporation_area(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin
get_r_insulation(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin
get_r_flesh(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin

function surface_area(shape::DesertIguana)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    total = Unitful.uconvert(u"m^2", (10.4713 * Unitful.ustrip(mass_g) ^ 0.688)u"cm^2")
    ventral = Unitful.uconvert(u"m^2", (0.425 * Unitful.ustrip(mass_g) ^ 0.85)u"cm^2")
    return (; total, ventral)
end

function silhouette_area(shape::DesertIguana)
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
function mammal_fur_area(body::AbstractBody)
    fur_multiplier = get_total_area(body) / get_skin_area(body)
    return Unitful.uconvert(u"m^2", (1110.0 * Unitful.ustrip(u"kg", body.shape.mass) ^ 0.65 * fur_multiplier)u"cm^2")
end
