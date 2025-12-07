"""
    Cylinder <: AbstractShape

A cylindrical organism shape.
"""
mutable struct Cylinder{M,D,B} <: AbstractShape
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

function geometry(shape::Cylinder, fur::Fur)
    volume = shape.mass / shape.density
    radius_skin = (volume / (shape.b * π * 2))^(1 / 3)
    radius_fur = radius_skin + fur.thickness
    length_skin = shape.b * radius_skin * 2
    length_fur = 2 * radius_fur
    total = surface_area(shape, radius_fur, length_fur)
    skin = surface_area(shape, radius_skin, length_skin)
    area_hair = hair_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = radius_fur * 2 #volume^(1 / 3) 
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur, length_skin, length_fur), (; total, skin, convection))
end

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
    area_hair = hair_area(fur.fibre_diameter, fur.fibre_density, skin)
    convection = skin - area_hair
    characteristic_dimension = radius_fur * 2 #volume^(1 / 3) 
    return Geometry(volume, characteristic_dimension, (; radius_skin, radius_fur, length_skin, length_fur, fat), (; total, skin, convection))
end

# function surface_area(shape::Cylinder, body::AbstractBody)
#     r = body.geometry.length.radius
#     l = body.geometry.length.length
#     surface_area(shape, r, l)
# end
surface_area(shape::Cylinder, r, l) = 2 * π * r * l + 2 * π * r^2

# silhouette area functions

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

function silhouette_area(shape::Cylinder, insulation::Union{Naked,Fat}, body::AbstractBody)
    r = body.geometry.length.radius
    l = body.geometry.length.length
    normal = 2 * r * l
    parallel = π *r^2
    return (; normal, parallel)
end

function silhouette_area(shape::Cylinder, insulation::Union{Fur,CompositeInsulation}, body::AbstractBody)
    r = body.geometry.length.radius_fur
    l = body.geometry.length.length_fur
    normal = 2 * r * l
    parallel = π *r^2
    return (; normal, parallel)
end

# area and radii functions

# naked
total_area(shape::Cylinder, insulation::Naked, body) = body.geometry.area.total
skin_area(shape::Cylinder, insulation::Naked, body) = body.geometry.area.total
evaporation_area(shape::Cylinder, insulation::Naked, body) = body.geometry.area.total

skin_radius(shape::Cylinder, insulation::Naked, body) = body.geometry.length.radius
insulation_radius(shape::Cylinder, insulation::Naked, body) = body.geometry.length.radius
flesh_radius(shape::Cylinder, insulation::Naked, body) = body.geometry.length.radius

# fur
total_area(shape::Cylinder, insulation::Fur, body) = body.geometry.area.total
skin_area(shape::Cylinder, insulation::Fur, body) = body.geometry.area.skin
evaporation_area(shape::Cylinder, insulation::Fur, body) = body.geometry.area.convection

skin_radius(shape::Cylinder, insulation::Fur, body) = body.geometry.length.radius_skin
insulation_radius(shape::Cylinder, insulation::Fur, body) = body.geometry.length.radius_fur
flesh_radius(shape::Cylinder, insulation::Fur, body) = body.geometry.length.radius_skin

# fat
total_area(shape::Cylinder, insulation::Fat, body) = body.geometry.area.total
skin_area(shape::Cylinder, insulation::Fat, body) = body.geometry.area.total
evaporation_area(shape::Cylinder, insulation::Fat, body) = body.geometry.area.total

skin_radius(shape::Cylinder, insulation::Fat, body) = body.geometry.length.radius
insulation_radius(shape::Cylinder, insulation::Fat, body) = body.geometry.length.radius
flesh_radius(shape::Cylinder, insulation::Fat, body) = body.geometry.length.radius - body.geometry.length.fat

# fur and fat
total_area(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.area.total
skin_area(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.area.skin
evaporation_area(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.area.convection

skin_radius(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin
insulation_radius(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.length.radius_fur
flesh_radius(shape::Cylinder, insulation::CompositeInsulation, body) = body.geometry.length.radius_skin - body.geometry.length.fat