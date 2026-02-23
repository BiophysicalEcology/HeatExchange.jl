"""
    net_metabolic_heat(; body, conductivities, core_temperature, skin_temperature)

Calculate net metabolic heat generation required to maintain core-to-skin temperature gradient.

Uses shape-specific heat conduction equations through flesh and fat layers to compute
the metabolic heat production needed.

# Keywords
- `body::AbstractBody`: Body geometry
- `conductivities::ThermalConductivities`: Thermal conductivities (flesh, fat, insulation)
- `core_temperature`: Core body temperature
- `skin_temperature`: Skin temperature

# Returns
- `net_generated_flux`: Net metabolic heat generation (W)
"""
function net_metabolic_heat(; body::AbstractBody, conductivities::ThermalConductivities, core_temperature, skin_temperature)
    net_metabolic_heat(shape(body), body, conductivities, core_temperature, skin_temperature)
end
function net_metabolic_heat(
    shape::Union{Cylinder,Plate}, body::AbstractBody, conductivities::ThermalConductivities, core_temperature, skin_temperature
)
    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    net_generated_flux = (core_temperature - skin_temperature) / (
            (r_flesh ^ 2 / (4 * conductivities.flesh * volume)) +
            ((r_flesh^2 / (2 * conductivities.fat * volume)) * log(r_skin / r_flesh))
        )
    return net_generated_flux
end
function net_metabolic_heat(
    shape::Sphere, body::AbstractBody, conductivities::ThermalConductivities, core_temperature, skin_temperature
)
    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    net_generated_flux = (core_temperature - skin_temperature) / (
            (r_flesh ^ 2 / (6 * conductivities.flesh * volume)) +
            ((r_flesh^3 / (3 * conductivities.fat * volume)) * ((r_skin - r_flesh)/(r_flesh * r_skin)))
        )
    return net_generated_flux
end
function net_metabolic_heat(
    shape::Ellipsoid, body::AbstractBody, conductivities::ThermalConductivities, core_temperature, skin_temperature
)
    volume = flesh_volume(body)
    a_semi_major = body.geometry.length.a_semi_major_skin
    b_semi_minor = body.geometry.length.b_semi_minor_skin
    c_semi_minor = body.geometry.length.c_semi_minor_skin
    fat = body.geometry.length.fat
    a_semi_major_flesh = a_semi_major - fat
    b_semi_minor_flesh = b_semi_minor - fat
    c_semi_minor_flesh = c_semi_minor - fat

    a_square = min(a_semi_major_flesh^2, a_semi_major^2)
    b_square = min(b_semi_minor_flesh^2, b_semi_minor^2)
    c_square = min(c_semi_minor_flesh^2, c_semi_minor^2)

    ssqg =
        (a_square * b_square * c_square) /
        (a_square * b_square + a_square * c_square + b_square * c_square)

    bs = b_semi_minor
    bg = min(b_semi_minor, b_semi_minor_flesh)

    net_generated_flux = (core_temperature - skin_temperature) / (
            (ssqg / (2 * conductivities.flesh * volume)) +
            (((((3 * ssqg)^0.5)^3) / (3 * conductivities.fat * volume)) * ((bs - bg) / (bg * bs)))
        )
    return net_generated_flux
end
