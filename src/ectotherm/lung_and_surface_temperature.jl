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

function surface_and_lung_temperature(shape::Cylinder, body::AbstractBody, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    flesh_radius = body.geometry.length[2]
    surface_temperature = core_temperature - specific_metabolic_heat_production * flesh_radius ^ 2 / (4 * flesh_conductivity)
    lung_temperature = (specific_metabolic_heat_production * flesh_radius ^ 2) / (8 * flesh_conductivity) + surface_temperature

    return (; surface_temperature, lung_temperature)
end

function surface_and_lung_temperature(shape::DesertIguana, body::AbstractBody, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    flesh_radius = body.geometry.length[1]
    surface_temperature = core_temperature - specific_metabolic_heat_production * flesh_radius ^ 2 / (4 * flesh_conductivity)
    lung_temperature = (specific_metabolic_heat_production * flesh_radius ^ 2) / (8 * flesh_conductivity) + surface_temperature

    return (; surface_temperature, lung_temperature)
end

function surface_and_lung_temperature(shape::LeopardFrog, body::AbstractBody, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    flesh_radius = body.geometry.length[1]
    surface_temperature = core_temperature - specific_metabolic_heat_production * flesh_radius ^ 2 / (4 * flesh_conductivity)
    lung_temperature = (specific_metabolic_heat_production * flesh_radius ^ 2) / (8 * flesh_conductivity) + surface_temperature

    return (; surface_temperature, lung_temperature)
end

function surface_and_lung_temperature(shape::Ellipsoid, body::AbstractBody, flesh_conductivity, specific_metabolic_heat_production, core_temperature)
    a = body.geometry.length[1] ^ 2
    b = body.geometry.length[2] ^ 2
    c = body.geometry.length[3] ^ 2
    x = ((a * b * c) / (a * b + a * c + b * c))
    surface_temperature = core_temperature - (specific_metabolic_heat_production / (2 * flesh_conductivity)) * x
    lung_temperature = (specific_metabolic_heat_production / (4 * flesh_conductivity)) * x + surface_temperature

    return (; surface_temperature, lung_temperature)
end
