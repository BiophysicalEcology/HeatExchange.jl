"""
    Tsurf_and_Tlung(; body, k_flesh, Q_gen_spec, T_core)

Calculate surface and lung temperatures from core temperature and metabolic heat generation.

Uses shape-specific heat conduction equations (cylinder, ellipsoid, etc.) to compute
temperature gradients from core to surface.

# Keywords
- `body::AbstractBody`: Body geometry
- `k_flesh`: Thermal conductivity of flesh
- `Q_gen_spec`: Specific metabolic heat generation (per unit volume)
- `T_core`: Core body temperature

# Returns
NamedTuple with:
- `T_surface`: Surface temperature
- `T_lung`: Lung temperature (intermediate between core and surface)
"""
function Tsurf_and_Tlung(; body::AbstractBody, k_flesh, Q_gen_spec, T_core)
    return Tsurf_and_Tlung(body, k_flesh, Q_gen_spec, T_core)
end

function Tsurf_and_Tlung(body::AbstractBody, k_flesh, Q_gen_spec, T_core)
    Tsurf_and_Tlung(shape(body), body, k_flesh, Q_gen_spec, T_core)
end

function Tsurf_and_Tlung(shape::Cylinder, body::AbstractBody, k_flesh, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.length[2]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_flesh)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_flesh) + T_surface

    return (; T_surface, T_lung)
end

function Tsurf_and_Tlung(shape::DesertIguana, body::AbstractBody, k_flesh, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.length[1]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_flesh)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_flesh) + T_surface

    return (; T_surface, T_lung)
end

function Tsurf_and_Tlung(shape::LeopardFrog, body::AbstractBody, k_flesh, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.length[1]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_flesh)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_flesh) + T_surface

    return (; T_surface, T_lung)
end

function Tsurf_and_Tlung(shape::Ellipsoid, body::AbstractBody, k_flesh, Q_gen_spec, T_core)
    a = body.geometry.length[1] ^ 2
    b = body.geometry.length[2] ^ 2
    c = body.geometry.length[3] ^ 2
    x = ((a * b * c) / (a * b + a * c + b * c))
    T_surface = T_core - (Q_gen_spec / (2 * k_flesh)) * x
    T_lung = (Q_gen_spec / (4 * k_flesh)) * x + T_surface

    return (; T_surface, T_lung)
end
