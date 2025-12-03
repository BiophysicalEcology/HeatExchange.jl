"""
    Plate <: Shape

A flat plate-shaped organism shape.
"""
mutable struct Plate{M,D,B,C} <: Shape
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

function surface_area(shape::Plate, body)
    length = body.geometry.length.length
    width = body.geometry.length.width
    height = body.geometry.length.height
    surface_area(shape, length, width, height)
end
surface_area(shape::Plate, length, width, height) = length * width * 2 + length * height * 2 + width * height * 2

# silhouette area functions

function silhouette_area(shape::Plate, insulation::Union{Naked,Fat}, body::AbstractBody)
    length = body.geometry.length.length
    width = body.geometry.length.width
    height = body.geometry.length.height
    normal = max(length * width, length * height, height * width)
    parallel = min(length * width, length * height, height * width)
    return (; normal, parallel)
end

function silhouette_area(shape::Plate, insulation::Union{Fur,CompositeInsulation}, body::AbstractBody)
    length = body.geometry.length.length_fur
    width = body.geometry.length.width_fur
    height = body.geometry.length.height_fur
    normal = max(length * width, length * height, height * width)
    parallel = min(length * width, length * height, height * width)
    return (; normal, parallel)
end

# area and radii functions

# naked
get_total_area(shape::Plate, insulation::Naked, body) = body.geometry.area.total
get_skin_area(shape::Plate, insulation::Naked, body) = body.geometry.area.total
get_evaporation_area(shape::Plate, insulation::Naked, body) = body.geometry.area.total

get_r_skin(shape::Plate, insulation::Naked, body) = body.geometry.length.width / 2
get_r_insulation(shape::Plate, insulation::Naked, body) = body.geometry.length.width / 2
get_r_flesh(shape::Plate, insulation::Naked, body) = body.geometry.length.width / 2

# fur
get_total_area(shape::Plate, insulation::Fur, body) = body.geometry.area.total
get_skin_area(shape::Plate, insulation::Fur, body) = body.geometry.area.skin
get_evaporation_area(shape::Plate, insulation::Fur, body) = body.geometry.area.convection

get_r_skin(shape::Plate, insulation::Fur, body) = body.geometry.length.width / 2
get_r_insulation(shape::Plate, insulation::Fur, body) = body.geometry.length.width_fur / 2
get_r_flesh(shape::Plate, insulation::Fur, body) = body.geometry.length.width / 2

# fat
get_total_area(shape::Plate, insulation::Fat, body) = body.geometry.area.total
get_skin_area(shape::Plate, insulation::Fat, body) = body.geometry.area.total
get_evaporation_area(shape::Plate, insulation::Fat, body) = body.geometry.area.total

get_r_skin(shape::Plate, insulation::Fat, body) = body.geometry.length.width / 2
get_r_insulation(shape::Plate, insulation::Fat, body) = body.geometry.length.width / 2
get_r_flesh(shape::Plate, insulation::Fat, body) = body.geometry.length.width / 2 - body.geometry.length.fat

# fur plus fat
get_total_area(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.area.total
get_skin_area(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.area.skin
get_evaporation_area(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.area.convection

get_r_skin(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width / 2
get_r_insulation(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width_fur / 2
get_r_flesh(shape::Plate, insulation::CompositeInsulation, body) = body.geometry.length.width / 2 - body.geometry.length.fat