mean_skin_temperature(; body::AbstractBody, insulation, insulation_pars, Q_env, Q_evap_skin,
    ks::ThermalConductivities, T_core, T_insulation_calc, T_ins_compressed, cds::ConductanceCoeffs, conduction_fraction) =
        mean_skin_temperature(shape(body), body, insulation, insulation_pars, Q_env, Q_evap_skin,
            ks, T_core, T_insulation_calc, T_ins_compressed, cds, conduction_fraction)

function mean_skin_temperature(shape::Union{Cylinder,Plate}, body, insulation, insulation_pars, Q_env, Q_evap_skin, ks::ThermalConductivities, T_core, T_insulation_calc, T_ins_compressed, cds::ConductanceCoeffs, conduction_fraction)
    (; k_flesh, k_fat) = ks
    (; cd1, cd2, cd3) = cds
    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    k_compressed = insulation.insulation_conductivity_compressed
    if conduction_fraction < 1
        T_skin_calc1 = T_core - (((Q_env + Q_evap_skin) * r_flesh^2) / (4 * k_flesh * volume)) - 
            (((Q_env + Q_evap_skin) * r_flesh^2) / (2 * k_fat * volume)) * log(r_skin / r_flesh)
        T_skin_calc2 = ((Q_env * r_flesh^2) / (2 * cd1 * volume)) + ((T_ins_compressed * cd2) / cd1) + 
            ((T_insulation_calc * cd3) / cd1)
    else
        T_skin_calc1 = T_core - ((Q_env * r_flesh^2) / (4 * k_flesh * volume)) - 
            ((Q_env * r_flesh^2.) / (2 * k_fat * volume)) * log(r_skin / r_flesh)
        T_skin_calc2 = (((Q_env * r_flesh^2) / (2 * k_compressed * volume)) * log(r_compressed / r_skin)) + 
            T_ins_compressed
    end
    T_skin_mean = (T_skin_calc1 + T_skin_calc2) / 2
    return (; T_skin_mean, T_skin_calc1)
end

function mean_skin_temperature(shape::Sphere, body, insulation, insulation_pars, Q_env, Q_evap_skin,
        ks::ThermalConductivities, T_core, T_insulation_calc, T_ins_compressed, cds::ConductanceCoeffs, conduction_fraction)
    (; k_flesh, k_fat) = ks
    (; cd1, cd2, cd3) = cds
    k_compressed = insulation.insulation_conductivity_compressed
    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    if conduction_fraction < 1
        T_skin_calc1 = T_core - (((Q_env + Q_evap_skin) * r_flesh^2) / (6 * k_flesh * volume)) - 
            (((Q_env + Q_evap_skin) * r_flesh^3.) / (3 * k_fat * volume)) * ((r_skin - r_flesh) / (r_skin * r_flesh))
        T_skin_calc2 = ((Q_env * r_flesh^3) / (3 * cd1 * volume * r_skin)) + ((T_ins_compressed * cd2) / cd1) + 
            ((T_insulation_calc * cd3) / cd1)
    else
        T_skin_calc1 = T_core - ((Q_env * r_flesh^2) / (6 * k_flesh * volume)) - 
            ((Q_env * r_flesh^3) / (3 * k_fat * volume)) * ((r_skin - r_flesh) / (r_skin * r_flesh))
        T_skin_calc2 = (((Q_env * r_flesh^3) / (3 * k_compressed * volume)) * ((r_compressed - r_skin) / 
            (r_compressed * r_skin))) + T_ins_compressed
    end
    T_skin_mean = (T_skin_calc1 + T_skin_calc2) / 2
    return (; T_skin_mean, T_skin_calc1)
end

function mean_skin_temperature(shape::Ellipsoid, body, insulation, insulation_pars, Q_env, Q_evap_skin,
        ks::ThermalConductivities, T_core, T_insulation_calc, T_ins_compressed, cds::ConductanceCoeffs, conduction_fraction)
    (; k_flesh, k_fat) = ks
    (; cd1, cd2, cd3) = cds
    volume = flesh_volume(body)
    k_compressed = insulation.insulation_conductivity_compressed
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

    ssqg = (a_square * b_square * c_square) / 
        (a_square * b_square + a_square * c_square + b_square * c_square)

    bs = b_semi_minor
    bl_compressed = b_semi_minor + insulation_pars.insulation_depth_compressed
    bg = min(b_semi_minor, b_semi_minor_flesh)

    if conduction_fraction < 1
        T_skin_calc1 = T_core - (((Q_env + Q_evap_skin) * ssqg) / (2 * k_flesh * volume)) - 
            (((Q_env + Q_evap_skin) * (((3 * ssqg)^0.5)^3)) / (3 * k_fat * volume)) * ((bs - bg) / (bs * bg))
        T_skin_calc2 = ((Q_env * (((3 * ssqg)^0.5)^3)) / (3 * cd1 * volume * bs)) + 
            ((T_ins_compressed * cd2) / cd1) + ((T_insulation_calc * cd3) / cd1)
    else
        T_skin_calc1 = T_core - ((Q_env * ssqg) / (2 * k_flesh * volume)) - 
            ((Q_env * (((3 * ssqg)^0.5)^3)) / (3 * k_fat * volume)) * ((bs - bg) / (bs * bg))
        T_skin_calc2 = (((Q_env * (((3 * ssqg)^0.5)^3)) / 
            (3 * k_compressed * volume)) * ((bl_compressed - bs) / (bl_compressed * bs))) + 
            T_ins_compressed
    end
        T_skin_mean = (T_skin_calc1 + T_skin_calc2) / 2
    return (; T_skin_mean, T_skin_calc1)
end
