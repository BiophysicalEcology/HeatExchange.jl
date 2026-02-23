"""
    Tsurf_and_Tlung(; body, flesh_conductivity, generated_specific_flux, core_temperature)

Calculate surface and lung temperatures from core temperature and metabolic heat generation.

Uses shape-specific heat conduction equations (cylinder, ellipsoid, etc.) to compute
temperature gradients from core to surface.

# Keywords
- `body::AbstractBody`: Body geometry
- `flesh_conductivity`: Thermal conductivity of flesh
- `generated_specific_flux`: Specific metabolic heat generation (per unit volume)
- `core_temperature`: Core body temperature

# Returns
NamedTuple with:
- `surface_temperature`: Surface temperature
- `lung_temperature`: Lung temperature (intermediate between core and surface)
"""
function Tsurf_and_Tlung(; body::AbstractBody, flesh_conductivity, generated_specific_flux, core_temperature)
    return Tsurf_and_Tlung(body, flesh_conductivity, generated_specific_flux, core_temperature)
end

function Tsurf_and_Tlung(body::AbstractBody, flesh_conductivity, generated_specific_flux, core_temperature)
    Tsurf_and_Tlung(shape(body), body, flesh_conductivity, generated_specific_flux, core_temperature)
end

function Tsurf_and_Tlung(shape::Cylinder, body::AbstractBody, flesh_conductivity, generated_specific_flux, core_temperature)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    flesh_radius = body.geometry.length[2]
    surface_temperature = core_temperature - generated_specific_flux * flesh_radius ^ 2 / (4 * flesh_conductivity)
    lung_temperature = (generated_specific_flux * flesh_radius ^ 2) / (8 * flesh_conductivity) + surface_temperature

    return (; surface_temperature, lung_temperature)
end

function Tsurf_and_Tlung(shape::DesertIguana, body::AbstractBody, flesh_conductivity, generated_specific_flux, core_temperature)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    flesh_radius = body.geometry.length[1]
    surface_temperature = core_temperature - generated_specific_flux * flesh_radius ^ 2 / (4 * flesh_conductivity)
    lung_temperature = (generated_specific_flux * flesh_radius ^ 2) / (8 * flesh_conductivity) + surface_temperature

    return (; surface_temperature, lung_temperature)
end

function Tsurf_and_Tlung(shape::LeopardFrog, body::AbstractBody, flesh_conductivity, generated_specific_flux, core_temperature)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    flesh_radius = body.geometry.length[1]
    surface_temperature = core_temperature - generated_specific_flux * flesh_radius ^ 2 / (4 * flesh_conductivity)
    lung_temperature = (generated_specific_flux * flesh_radius ^ 2) / (8 * flesh_conductivity) + surface_temperature

    return (; surface_temperature, lung_temperature)
end

function Tsurf_and_Tlung(shape::Ellipsoid, body::AbstractBody, flesh_conductivity, generated_specific_flux, core_temperature)
    a = body.geometry.length[1] ^ 2
    b = body.geometry.length[2] ^ 2
    c = body.geometry.length[3] ^ 2
    x = ((a * b * c) / (a * b + a * c + b * c))
    surface_temperature = core_temperature - (generated_specific_flux / (2 * flesh_conductivity)) * x
    lung_temperature = (generated_specific_flux / (4 * flesh_conductivity)) * x + surface_temperature

    return (; surface_temperature, lung_temperature)
end
