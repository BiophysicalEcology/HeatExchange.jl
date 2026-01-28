"""
    compressed_radiant_temperature(; body, insulation, insulation_pars, ks, side, cd, T_core, T_substrate)

Calculate temperature at compressed insulation interface when body is in contact with substrate.

Used when conduction_fraction = 1 (full contact with ground), computing the temperature
at the compressed insulation-substrate interface.

# Keywords
- `body::AbstractBody`: Body geometry
- `insulation::InsulationProperties`: Computed insulation properties
- `insulation_pars::InsulationParameters`: Insulation parameters
- `ks::ThermalConductivities`: Thermal conductivities (flesh, fat)
- `side`: Body side (`:dorsal` or `:ventral`)
- `cd`: Substrate conductance coefficient
- `T_core`: Core body temperature
- `T_substrate`: Substrate temperature

# Returns
NamedTuple with:
- `cf1`: Compression factor coefficient
- `T_ins_compressed`: Temperature at compressed insulation interface
"""
function compressed_radiant_temperature(;
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    ks::ThermalConductivities,
    side,
    cd,
    T_core,
    T_substrate,
)
    compressed_radiant_temperature(
        shape(body), body, insulation, insulation_pars, ks, side, cd, T_core, T_substrate
    )
end

function compressed_radiant_temperature(
    shape::Union{Cylinder,Plate},
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    ks::ThermalConductivities,
    side,
    cd,
    T_core,
    T_substrate,
)
    (; k_flesh, k_fat) = ks
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)
    r_compressed = r_skin + insulation_pars.depth_compressed
    k_compressed = insulation.conductivity_compressed

    cf1 = (2 * π * k_compressed * length) / (log(r_compressed / r_skin))
    dv5 =
        1 +
        ((cf1 * r_flesh^2) / (4 * k_flesh * volume)) +
        ((cf1 * r_flesh^2) / (2 * k_fat * volume)) * log(r_skin / RFLESH)
    T_ins_compressed_calc1 = (cf1 / dv5) * T_core + cd * T_substrate
    T_ins_compressed_calc2 = cd + cf1 / dv5
    T_ins_compressed = T_ins_compressed_calc1 / T_ins_compressed_calc2
    return (; cf1, T_ins_compressed)
end

function compressed_radiant_temperature(
    shape::Sphere,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    ks::ThermalConductivities,
    side,
    cd,
    T_core,
    T_substrate,
)
    (; k_flesh, k_fat) = ks
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)
    r_compressed = r_skin + insulation_pars.depth_compressed
    k_compressed = insulation.conductivity_compressed

    cf1 = (4 * π * k_compressed * r_compressed) / (r_compressed - r_skin)
    dv5 =
        1 +
        ((cf1 * r_flesh^2.0) / (6 * k_flesh * volume)) +
        ((cf1 * r_flesh^3) / (3 * k_fat * volume)) *
        ((r_skin - r_flesh) / (r_skin - r_flesh))
    T_ins_compressed_calc1 = (cf1 / dv5) * T_core + cd * T_substrate
    T_ins_compressed_calc2 = cd + cf1 / dv5
    T_ins_compressed = T_ins_compressed_calc1 / T_ins_compressed_calc2
    return (; cf1, T_ins_compressed)
end

function compressed_radiant_temperature(
    shape::Ellipsoid,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    ks::ThermalConductivities,
    side,
    cd,
    T_core,
    T_substrate,
)
    (; k_flesh, k_fat) = ks
    volume = flesh_volume(body)
    insulation_depth = if side == :dorsal
        insulation_pars.dorsal.depth
    else
        insulation_pars.ventral.depth
    end
    k_compressed = insulation.conductivity_compressed

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

    cf1 =
        (3 * k_compressed * volume * bl_compressed * bs) /
        ((((3 * ssqg)^0.5)^3) * (bl - bs))
    dv5 =
        1 +
        ((cf1 * ssqg) / (2 * k_flesh * volume)) +
        ((cf1 * (((3 * ssqg)^0.5)^3)) / (3 * k_fat * volume)) * ((bs - bg) / (bs * bg))
    T_ins_compressed_calc1 = (cf1 / dv5) * T_core + cd * T_substrate
    T_ins_compressed_calc2 = cd + cf1 / dv5
    T_ins_compressed = T_ins_compressed_calc1 / T_ins_compressed_calc2
    return (; cf1, T_ins_compressed)
end
