"""Ellipsoid shape factor `a²b²c²/(a²b²+a²c²+b²c²)` for semi-axes `a`, `b`, `c`."""
function ellipsoid_shape_factor(a, b, c)
    a2, b2, c2 = a^2, b^2, c^2
    (a2 * b2 * c2) / (a2 * b2 + a2 * c2 + b2 * c2)
end

"""Internal-temperature-gradient shape factor for `body` (shape-dispatched)."""
internal_gradient_shape_factor(body::AbstractBody) = internal_gradient_shape_factor(shape(body), body)

function internal_gradient_shape_factor(::Ellipsoid, body::AbstractBody)
    a, b, c = body.geometry.length[1], body.geometry.length[2], body.geometry.length[3]
    ellipsoid_shape_factor(a, b, c)
end

internal_gradient_shape_factor(::Cylinder, body::AbstractBody) = body.geometry.length[2]^2

internal_gradient_shape_factor(::Union{DesertIguana,LeopardFrog}, body::AbstractBody) = body.geometry.length[1]^2

"""
    surface_and_lung_temperature(; body, flesh_conductivity, specific_metabolic_heat_production, core_temperature)

Calculate surface and lung temperatures from core temperature and metabolic heat generation.

Uses shape-specific heat conduction equations (cylinder, ellipsoid, etc.) to compute
temperature gradients from core to surface.

# Keywords
- `body::AbstractBody`: Body geometry
- `flesh_conductivity`: Thermal conductivity of flesh
- `specific_metabolic_heat_production`: Specific metabolic heat generation (per unit volume)
- `core_temperature`: Core body temperature

# Returns
NamedTuple with:
- `surface_temperature`: Surface temperature
- `lung_temperature`: Lung temperature (intermediate between core and surface)
"""
function surface_and_lung_temperature(; body::AbstractBody, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
    return surface_and_lung_temperature(body, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
end

function surface_and_lung_temperature(body::AbstractBody, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
    surface_and_lung_temperature(shape(body), body, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
end

function surface_and_lung_temperature(::Union{Cylinder,DesertIguana,LeopardFrog}, body::AbstractBody, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    shape_factor = internal_gradient_shape_factor(body)
    surface_temperature = core_temperature - specific_metabolic_heat_production * shape_factor / (4 * flesh_conductivity)
    lung_temperature = (specific_metabolic_heat_production * shape_factor) / (8 * flesh_conductivity) + surface_temperature

    return (; surface_temperature, lung_temperature)
end

function surface_and_lung_temperature(::Plate, body::AbstractBody, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
    # flat slab: half-thickness h = height/2 (shortest dimension)
    # from plane-wall solution (Bird, Stewart & Lightfoot, Transport Phenomena)
    h = body.geometry.length.height_skin / 2
    surface_temperature = core_temperature - specific_metabolic_heat_production * h ^ 2 / (2 * flesh_conductivity)
    lung_temperature = (specific_metabolic_heat_production * h ^ 2) / (4 * flesh_conductivity) + surface_temperature

    return (; surface_temperature, lung_temperature)
end

function surface_and_lung_temperature(::Ellipsoid, body::AbstractBody, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
    shape_factor = internal_gradient_shape_factor(body)
    surface_temperature = core_temperature - (specific_metabolic_heat_production / (2 * flesh_conductivity)) * shape_factor
    lung_temperature = (specific_metabolic_heat_production / (4 * flesh_conductivity)) * shape_factor + surface_temperature

    return (; surface_temperature, lung_temperature)
end
