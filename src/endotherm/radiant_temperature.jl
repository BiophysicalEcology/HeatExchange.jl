function radiant_temperature(;
    body::AbstractBody,
    insulation::InsulationOutput,
    insulation_pars::InsulationParameters,
    org_temps::OrganismTemperatures,
    ks::ThermalConductivities,
    side,
    cd,
    longwave_depth_fraction,
    conduction_fraction,
    Q_evap,
    T_substrate,
)
    radiant_temperature(
        shape(body),
        body,
        insulation,
        insulation_pars,
        org_temps,
        ks,
        side,
        cd,
        longwave_depth_fraction,
        conduction_fraction,
        Q_evap,
        T_substrate,
    )
end

function radiant_temperature(
    shape::Union{Cylinder,Plate},
    body::AbstractBody,
    insulation::InsulationOutput,
    insulation_pars::InsulationParameters,
    org_temps::OrganismTemperatures,
    ks::ThermalConductivities,
    side,
    cd,
    longwave_depth_fraction,
    conduction_fraction,
    Q_evap,
    T_substrate,
)
    (; T_core, T_skin, T_insulation) = org_temps
    (; k_flesh, k_fat, k_insulation) = ks

    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    r_insulation = insulation_radius(body)
    insulation_depth = insulation.insulation_depths[side + 1]
    r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation_depth
    k_compressed = insulation.insulation_conductivity_compressed
    r_compressed = r_skin + insulation_pars.insulation_depth_compressed
    length = body.geometry.length.length_skin

    compression_fraction =
        (conduction_fraction * 2 * π * k_compressed * length) / log(r_compressed / r_skin)
    T_ins_compressed = if conduction_fraction > 0
        (compression_fraction * T_skin + cd * T_substrate) / (cd + compression_fraction)
    else
        0.0u"K"
    end

    cd1 =
        (k_compressed / log(r_compressed / r_skin)) * conduction_fraction +
        (k_insulation / log(r_insulation / r_skin)) * (1 - conduction_fraction)
    cd2 = (k_compressed / log(r_compressed / r_skin)) * conduction_fraction
    cd3 = (k_insulation / log(r_insulation / r_skin)) * (1 - conduction_fraction)

    dv1 =
        1 +
        ((2 * π * length * r_flesh^2 * cd1) / (4 * k_flesh * volume)) +
        ((2 * π * length * r_flesh^2 * cd1) / (2 * k_fat * volume)) * log(r_skin / r_flesh)

    dv2 =
        Q_evap * ((r_flesh^2 * cd1) / (4 * k_flesh * volume)) +
        Q_evap * ((r_flesh^2 * cd1) / (2 * k_fat * volume)) * log(r_skin / r_flesh)

    dv3 =
        ((2 * π * length) / dv1) *
        (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) *
        r_flesh^2 / (2 * volume)

    dv4 = if longwave_depth_fraction < 1
        cd2 + (k_insulation / log(r_insulation / r_radiation)) * (1 - conduction_fraction)
    else
        1.0
    end

    T_radiant = if longwave_depth_fraction < 1
        dv3 / dv4 +
        (T_ins_compressed * cd2) / dv4 +
        (
            T_insulation * (
                (k_insulation / log(r_insulation / r_radiation)) *
                (1 - conduction_fraction)
            )
        ) / dv4
    else
        T_insulation
    end
    cds = ConductanceCoeffs(cd1, cd2, cd3)
    dvs = DivisorCoeffs(dv1, dv2, dv3, dv4)
    return (; T_radiant, T_ins_compressed, cds, dvs)
end

function radiant_temperature(
    shape::Sphere,
    body::AbstractBody,
    insulation::InsulationOutput,
    insulation_pars::InsulationParameters,
    org_temps::OrganismTemperatures,
    ks::ThermalConductivities,
    side,
    cd,
    longwave_depth_fraction,
    conduction_fraction,
    Q_evap,
    T_substrate,
)
    (; T_core, T_skin, T_insulation) = org_temps
    (; k_flesh, k_fat, k_insulation) = ks

    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    r_insulation = insulation_radius(body)
    insulation_depth = insulation.insulation_depths[side + 1]
    r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation_depth
    k_compressed = insulation.insulation_conductivity_compressed
    r_compressed = r_skin + insulation_pars.insulation_depth_compressed

    compression_fraction =
        (conduction_fraction * 4 * π * k_compressed * r_compressed * r_skin) /
        (r_compressed - r_skin)

    T_ins_compressed = if conduction_fraction > 0
        (compression_fraction * T_skin + cd * T_substrate) / (cd + compression_fraction)
    else
        0.0u"K"
    end

    cd1 =
        ((k_compressed * r_compressed) / (r_compressed - r_skin)) * conduction_fraction +
        ((k_insulation * r_insulation) / (r_insulation - r_skin)) *
        (1 - conduction_fraction)

    cd2 = ((k_compressed * r_compressed) / (r_compressed - r_skin)) * conduction_fraction
    cd3 =
        ((k_insulation * r_insulation) / (r_insulation - r_skin)) *
        (1 - conduction_fraction)

    dv1 =
        1 +
        ((4 * π * r_skin * r_flesh^2 * cd1) / (6 * k_flesh * volume)) +
        ((4 * π * r_skin * r_flesh^3 * cd1) / (3 * k_fat * volume)) *
        ((r_skin - r_flesh) / (r_flesh * r_skin))

    dv2 =
        Q_evap * ((r_flesh^2 * cd1) / (6 * k_flesh * volume)) +
        Q_evap *
        ((r_flesh^3 * cd1) / (3 * k_fat * volume)) *
        ((r_skin - r_flesh) / (r_flesh * r_skin))

    dv3 =
        ((4 * π * r_skin) / dv1) *
        (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) *
        r_flesh^3 / (3 * volume * r_radiation)

    dv4 = if longwave_depth_fraction < 1
        cd2 +
        ((k_insulation * r_insulation) / (r_insulation - r_radiation)) *
        (1 - conduction_fraction)
    else
        1.0
    end

    T_radiant = if longwave_depth_fraction < 1
        dv3 / dv4 +
        (T_ins_compressed * cd2) / dv4 +
        (
            T_insulation * (
                (k_insulation * r_insulation) / (r_insulation - r_radiation) *
                (1 - conduction_fraction)
            )
        ) / dv4
    else
        T_insulation
    end
    cds = ConductanceCoeffs(cd1, cd2, cd3)
    dvs = DivisorCoeffs(dv1, dv2, dv3, dv4)
    return (; T_radiant, T_ins_compressed, cds, dvs)
end

function radiant_temperature(
    shape::Ellipsoid,
    body::AbstractBody,
    insulation::InsulationOutput,
    insulation_pars::InsulationParameters,
    org_temps::OrganismTemperatures,
    ks::ThermalConductivities,
    side,
    cd,
    longwave_depth_fraction,
    conduction_fraction,
    Q_evap,
    T_substrate,
)
    (; T_core, T_skin, T_insulation) = org_temps
    (; k_flesh, k_fat, k_insulation) = ks

    volume = flesh_volume(body)

    a_semi_major = body.geometry.length.a_semi_major_skin
    b_semi_minor = body.geometry.length.b_semi_minor_skin
    c_semi_minor = body.geometry.length.c_semi_minor_skin
    fat = body.geometry.length.fat
    a_semi_major_flesh = a_semi_major - fat
    b_semi_minor_flesh = b_semi_minor - fat
    c_semi_minor_flesh = c_semi_minor - fat

    insulation_depth = insulation.insulation_depths[side + 1]
    k_compressed = insulation.insulation_conductivity_compressed
    bl_compressed = b_semi_minor + insulation_pars.insulation_depth_compressed

    a_square = min(a_semi_major_flesh^2, a_semi_major^2)
    b_square = min(b_semi_minor_flesh^2, b_semi_minor^2)
    c_square = min(c_semi_minor_flesh^2, c_semi_minor^2)

    ssqg =
        (a_square * b_square * c_square) /
        (a_square * b_square + a_square * c_square + b_square * c_square)

    bg = min(b_semi_minor, b_semi_minor_flesh)
    bs = b_semi_minor
    bl = b_semi_minor + insulation_depth
    br = bs + longwave_depth_fraction * insulation_depth

    compression_fraction =
        (conduction_fraction * 3 * k_compressed * volume * bl_compressed * bs) /
        ((sqrt(3 * ssqg))^3 * (bl_compressed - bs))

    T_ins_compressed = if conduction_fraction > 0.0
        (compression_fraction * T_skin + cd * T_substrate) / (cd + compression_fraction)
    else
        0.0u"K"
    end

    cd1 =
        ((k_compressed * bl_compressed) / (bl_compressed - bs)) * conduction_fraction +
        ((k_insulation * bl) / (bl - bs)) * (1 - conduction_fraction)
    cd2 = ((k_compressed * bl_compressed) / (bl_compressed - bs)) * conduction_fraction
    cd3 = ((k_insulation * bl) / (bl - bs)) * (1 - conduction_fraction)
    dv1 =
        1 +
        (3 * bs * ssqg * cd1) / (2 * k_flesh * (sqrt(3 * ssqg)^3)) +
        (bs * cd1) / k_fat * ((bs - bg) / (bs * bg))

    dv2 =
        Q_evap * ((ssqg * cd1) / (2 * k_flesh * volume)) +
        Q_evap * ((sqrt(3 * ssqg)^3 * cd1) / (3 * k_fat * volume)) * ((bs - bg) / (bs * bg))

    dv3 =
        (bs / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) / br

    dv4 = if longwave_depth_fraction < 1
        cd2 + ((k_insulation * bl) / (bl - br)) * (1 - conduction_fraction)
    else
        1.0
    end

    T_radiant = if longwave_depth_fraction < 1
        dv3 / dv4 +
        (T_ins_compressed * cd2) / dv4 +
        (T_insulation * ((k_insulation * bl) / (bl - br) * (1 - conduction_fraction))) / dv4
    else
        T_insulation
    end
    cds = ConductanceCoeffs(cd1, cd2, cd3)
    dvs = DivisorCoeffs(dv1, dv2, dv3, dv4)
    return (; T_radiant, T_ins_compressed, cds, dvs)
end
