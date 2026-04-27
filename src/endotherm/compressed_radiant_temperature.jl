"""
    compressed_radiant_temperature(; body, insulation, insulation_pars, conductivities, side, conductance_coefficient, core_temperature, substrate_temperature)

Calculate temperature at compressed insulation interface when body is in contact with substrate.

Used when conduction_fraction = 1 (full contact with ground), computing the temperature
at the compressed insulation-substrate interface.

# Keywords
- `body::AbstractBody`: Body geometry
- `insulation::InsulationProperties`: Computed insulation properties
- `insulation_pars::InsulationParameters`: Insulation parameters
- `conductivities::ThermalConductivities`: Thermal conductivities (flesh, fat)
- `side`: Body side (`:dorsal` or `:ventral`)
- `conductance_coefficient`: Substrate conductance coefficient
- `core_temperature`: Core body temperature
- `substrate_temperature`: Substrate temperature

# Returns
NamedTuple with:
- `compression_fraction`: Compression factor coefficient
- `compressed_insulation_temperature`: Temperature at compressed insulation interface
"""
function compressed_radiant_temperature(;
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    conductivities::ThermalConductivities,
    side,
    conductance_coefficient,
    core_temperature,
    substrate_temperature,
)
    compressed_radiant_temperature(
        shape(body), body, insulation, insulation_pars, conductivities, side, conductance_coefficient, core_temperature, substrate_temperature
    )
end

function compressed_radiant_temperature(
    shape::Union{Cylinder,Plate},
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    conductivities::ThermalConductivities,
    side,
    conductance_coefficient,
    core_temperature,
    substrate_temperature,
)
    volume = flesh_volume(body)
    length = body.geometry.length.length_skin
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    r_compressed = r_skin + insulation_pars.depth_compressed

    compression_fraction = (2 * π * insulation.conductivity_compressed * length) / (log(r_compressed / r_skin))
    geometric_divisor =
        1 +
        ((compression_fraction * r_flesh^2) / (4 * conductivities.flesh * volume)) +
        ((compression_fraction * r_flesh^2) / (2 * conductivities.fat * volume)) * log(r_skin / r_flesh)
    compressed_insulation_temperature_calc1 = (compression_fraction / geometric_divisor) * core_temperature + conductance_coefficient * substrate_temperature
    compressed_insulation_temperature_calc2 = conductance_coefficient + compression_fraction / geometric_divisor
    compressed_insulation_temperature = compressed_insulation_temperature_calc1 / compressed_insulation_temperature_calc2
    return (; compression_fraction, compressed_insulation_temperature)
end

function compressed_radiant_temperature(
    shape::Sphere,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    conductivities::ThermalConductivities,
    side,
    conductance_coefficient,
    core_temperature,
    substrate_temperature,
)
    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    r_compressed = r_skin + insulation_pars.depth_compressed

    compression_fraction = (4 * π * insulation.conductivity_compressed * r_compressed) / (r_compressed - r_skin)
    geometric_divisor =
        1 +
        ((compression_fraction * r_flesh^2.0) / (6 * conductivities.flesh * volume)) +
        ((compression_fraction * r_flesh^3) / (3 * conductivities.fat * volume)) *
        ((r_skin - r_flesh) / (r_flesh * r_skin))
    compressed_insulation_temperature_calc1 = (compression_fraction / geometric_divisor) * core_temperature + conductance_coefficient * substrate_temperature
    compressed_insulation_temperature_calc2 = conductance_coefficient + compression_fraction / geometric_divisor
    compressed_insulation_temperature = compressed_insulation_temperature_calc1 / compressed_insulation_temperature_calc2
    return (; compression_fraction, compressed_insulation_temperature)
end

function compressed_radiant_temperature(
    shape::Ellipsoid,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    conductivities::ThermalConductivities,
    side,
    conductance_coefficient,
    core_temperature,
    substrate_temperature,
)
    volume = flesh_volume(body)
    insulation_depth = if side == :dorsal
        insulation_pars.dorsal.depth
    else
        insulation_pars.ventral.depth
    end

    a_semi_major = body.geometry.length.a_semi_major
    b_semi_minor = body.geometry.length.b_semi_minor
    c_semi_minor = body.geometry.length.c_semi_minor

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
    bl = b_semi_minor + insulation_depth
    bl_compressed = b_semi_minor + insulation_pars.depth_compressed
    bg = min(b_semi_minor, b_semi_minor_flesh)

    compression_fraction =
        (3 * insulation.conductivity_compressed * volume * bl_compressed * bs) /
        ((((3 * ssqg)^0.5)^3) * (bl - bs))
    geometric_divisor =
        1 +
        ((compression_fraction * ssqg) / (2 * conductivities.flesh * volume)) +
        ((compression_fraction * (((3 * ssqg)^0.5)^3)) / (3 * conductivities.fat * volume)) * ((bs - bg) / (bs * bg))
    compressed_insulation_temperature_calc1 = (compression_fraction / geometric_divisor) * core_temperature + conductance_coefficient * substrate_temperature
    compressed_insulation_temperature_calc2 = conductance_coefficient + compression_fraction / geometric_divisor
    compressed_insulation_temperature = compressed_insulation_temperature_calc1 / compressed_insulation_temperature_calc2

    return (; compression_fraction, compressed_insulation_temperature)
end
