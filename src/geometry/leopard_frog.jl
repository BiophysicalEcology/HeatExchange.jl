"""
    LeopardFrog <: AbstractShape

A frog-shaped organism. Based on the leopard frog (Tracy 1976 Ecol. Monog.).
"""
struct LeopardFrog{M,D} <: AbstractShape
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

# area functions

function surface_area(shape::LeopardFrog)
    mass_g = Unitful.uconvert(u"g", shape.mass)
    Unitful.uconvert(u"m^2", (12.79 * Unitful.ustrip(mass_g) ^ 0.606)u"cm^2") # eq in Fig. 5
end

# silhouette area functions

function silhouette_area(shape::LeopardFrog, θ)
    area = surface_area(shape)
    pct = 1.38171e-6 * θ ^ 4 - 1.93335e-4 * θ ^ 3 + 4.75761e-3 * θ ^ 2 - 0.167912 * θ + 45.8228
    return pct * area / 100
end

# surface area and radii functions
get_total_area(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.area.total
get_skin_area(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.area.total
get_evaporation_area(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.area.total

# "radii" functions
get_r_skin(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin
get_r_insulation(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin
get_r_flesh(shape::LeopardFrog, insulation::Naked, body::AbstractBody) = body.geometry.length.r_skin