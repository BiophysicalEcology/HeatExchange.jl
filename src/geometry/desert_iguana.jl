"""
    DesertIguana <: AbstractShape

A lizard-shaped organism. Based on the desert iguana (Porter and Tracy 1984).
"""
struct DesertIguana{M,D} <: AbstractShape
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

# area functions

function surface_area(shape::DesertIguana)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    total = Unitful.uconvert(u"m^2", (10.4713 * Unitful.ustrip(mass_g) ^ 0.688)u"cm^2")
    ventral = Unitful.uconvert(u"m^2", (0.425 * Unitful.ustrip(mass_g) ^ 0.85)u"cm^2")
    return (; total, ventral)
end

# silhouette area functions

function silhouette_area(shape::DesertIguana, ::NormalToSun)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    return Unitful.uconvert(u"m^2", (3.798 * Unitful.ustrip(mass_g) ^ 0.683)u"cm^2")
end

function silhouette_area(shape::DesertIguana, ::ParallelToSun)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    return Unitful.uconvert(u"m^2", (0.694 * Unitful.ustrip(mass_g) ^ 0.743)u"cm^2")
end

function silhouette_area(shape::DesertIguana, ::Intermediate)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    normal = Unitful.uconvert(u"m^2", (3.798 * Unitful.ustrip(mass_g) ^ 0.683)u"cm^2")
    parallel = Unitful.uconvert(u"m^2", (0.694 * Unitful.ustrip(mass_g) ^ 0.743)u"cm^2")
    return (normal + parallel) * 0.5
end

# surface area and radii functions
total_area(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.area.total
skin_area(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.area.total
evaporation_area(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.area.total

skin_radius(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin
insulation_radius(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin
flesh_radius(shape::DesertIguana, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin
