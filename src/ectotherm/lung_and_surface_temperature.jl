"""
    Tsurf_and_Tlung(; kw...)
    Tsurf_and_Tlung(..) 

Computes ...

# Keywords
- `x`: x

"""
function Tsurf_and_Tlung(; body, k_flesh, Q_gen_spec, T_core)
    return Tsurf_and_Tlung(body, k_flesh, Q_gen_spec, T_core)
end
function Tsurf_and_Tlung(body::AbstractBody, k_flesh, Q_gen_spec, T_core)
    Tsurf_and_Tlung(shape(body), body, k_flesh, Q_gen_spec, T_core)
end
function Tsurf_and_Tlung(shape::Cylinder, body, k_flesh, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.length[2]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_flesh)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_flesh) + T_surface

    return (; T_surface, T_lung)
end
function Tsurf_and_Tlung(shape::DesertIguana, body, k_flesh, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.length[1]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_flesh)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_flesh) + T_surface

    return (; T_surface, T_lung)
end
function Tsurf_and_Tlung(shape::LeopardFrog, body, k_flesh, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.length[1]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_flesh)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_flesh) + T_surface

    return (; T_surface, T_lung)
end
function Tsurf_and_Tlung(shape::Ellipsoid, body, k_flesh, Q_gen_spec, T_core)
    a = body.geometry.length[1] ^ 2
    b = body.geometry.length[2] ^ 2
    c = body.geometry.length[3] ^ 2
    x = ((a * b * c) / (a * b + a * c + b * c))
    T_surface = T_core - (Q_gen_spec / (2 * k_flesh)) * x
    T_lung = (Q_gen_spec / (4 * k_flesh)) * x + T_surface

    return (; T_surface, T_lung)
end
