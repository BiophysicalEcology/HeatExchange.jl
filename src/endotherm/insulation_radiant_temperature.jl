insulation_radiant_temperature(; body::AbstractBody, insulation, insulation_pars, T_core,
    T_ins_compressed, env_temps::EnvironmentTemperatures, area_convection, hc, cd,
    k_insulation, Q_solar, Q_evap_insulation, Q_rads::RadiationCoeffs,
    cds::ConductanceCoeffs, dvs::DivisorCoeffs, longwave_depth_fraction, conduction_fraction, side) =
    insulation_radiant_temperature(shape(body), body, insulation, insulation_pars, T_core,
        T_ins_compressed, env_temps, area_convection, hc, cd, k_insulation, Q_solar, Q_evap_insulation,
        Q_rads, cds, dvs, longwave_depth_fraction, conduction_fraction, side)
        
function insulation_radiant_temperature(shape::Union{Cylinder,Plate}, body, insulation, insulation_pars,
    T_core, T_ins_compressed, env_temps::EnvironmentTemperatures, area_convection, hc,
    cd, k_insulation, Q_solar, Q_evap_insulation, Q_rads::RadiationCoeffs,
    cds::ConductanceCoeffs, dvs::DivisorCoeffs, longwave_depth_fraction, conduction_fraction, side)

    (; T_air, T_sky, T_ground, T_vegetation, T_bush, T_substrate) = env_temps
    (; cd1, cd2, cd3) = cds
    (; dv1, dv2, dv3, dv4) = dvs
    (; Q_rad1, Q_rad2, Q_rad3, Q_rad4) = Q_rads

    r_skin = skin_radius(body)
    r_insulation = insulation_radius(body)
    length = body.geometry.length.length_skin
    
    if side == 1
        r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation_pars.insulation_depth_dorsal
    else
        r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation_pars.insulation_depth_ventral
    end    

    if longwave_depth_fraction < 1
        T_ins1 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground - 
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins2 = ((2 * π * LEN) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_substrate - 
            Q_evap_insulation + Q_solar
        T_ins4 = (2 * π * length * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * 
            (((k_insulation / log(r_insulation / r_radiation)) * (1 - conduction_fraction)) / dv4) 
            + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(dv3 / dv4 + ((T_ins_compressed * cd2) / dv4) + (T_insulation_calc *
             ((k_insulation / log(r_insulation / r_radiation)) * (1 - conduction_fraction))) / dv4)
    else
        T_ins1 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + 
            Q_rad4 * T_ground
        T_ins2 = ((2 * π * length) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_substrate -
         Q_evap_insulation + Q_solar
        T_ins4 = (2 * π * length * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) +
            hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_insulation_calc
    end

    return(; T_insulation_calc, T_radiant2)
end

function insulation_radiant_temperature(shape::Sphere, body, insulation, insulation_pars,
    T_core, T_ins_compressed, env_temps::EnvironmentTemperatures, area_convection, hc,
    cd, k_insulation, Q_solar, Q_evap_insulation, Q_rads::RadiationCoeffs,
    cds::ConductanceCoeffs, dvs::DivisorCoeffs, longwave_depth_fraction, conduction_fraction, side)

    (; T_air, T_sky, T_ground, T_vegetation, T_bush, T_substrate) = env_temps
    (; cd1, cd2, cd3) = cds
    (; dv1, dv2, dv3, dv4) = dvs
    (; Q_rad1, Q_rad2, Q_rad3, Q_rad4) = Q_rads

    r_skin = skin_radius(body)
    r_insulation = insulation_radius(body)

    if side == 1
        r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation_pars.insulation_depth_dorsal
    else
        r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation_pars.insulation_depth_ventral
    end 

    if longwave_depth_fraction < 1
        T_ins1 = ((4 * π * r_skin) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground -
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_substrate - 
         Q_evap_insulation + Q_solar
        T_ins4 = (4 * π * r_skin * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) *
         ((((k_insulation * r_insulation) / (r_insulation - r_radiation)) * 
         (1 - conduction_fraction)) / dv4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(dv3 / dv4 + ((T_ins_compressed * cd2) / dv4) +
         (T_insulation_calc * (((k_insulation * r_insulation) / (r_insulation - r_radiation)) * 
         (1 - conduction_fraction))) / dv4)
    else
        T_ins1 = ((4 * π * r_skin) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_substrate - 
         Q_evap_insulation + Q_solar
        T_ins4 = (4 * π * r_skin * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_insulation_calc
    end
    return(; T_insulation_calc, T_radiant2)    
end

function insulation_radiant_temperature(shape::Ellipsoid, body, insulation, insulation_pars,
    T_core, T_ins_compressed, env_temps::EnvironmentTemperatures, area_convection, hc,
    cd, k_insulation, Q_solar, Q_evap_insulation, Q_rads::RadiationCoeffs,
    cds::ConductanceCoeffs, dvs::DivisorCoeffs, longwave_depth_fraction, conduction_fraction, side)

    (; T_air, T_sky, T_ground, T_vegetation, T_bush, T_substrate) = env_temps
    (; cd1, cd2, cd3) = cds
    (; dv1, dv2, dv3, dv4) = dvs
    (; Q_rad1, Q_rad2, Q_rad3, Q_rad4) = Q_rads

    volume = flesh_volume(body)

    a_semi_major = body.geometry.length.a_semi_major_skin
    b_semi_minor = body.geometry.length.b_semi_minor_skin
    c_semi_minor = body.geometry.length.c_semi_minor_skin

    fat = body.geometry.length.fat
    a_semi_major_flesh = a_semi_major - fat
    b_semi_minor_flesh = b_semi_minor - fat
    c_semi_minor_flesh = c_semi_minor - fat

    if side == 1
        insulation_depth = insulation_pars.insulation_depth_dorsal
    else
        insulation_depth = insulation_pars.insulation_depth_ventral
    end 

    a_square = min(a_semi_major_flesh^2, a_semi_major^2)
    b_square = min(b_semi_minor_flesh^2, b_semi_minor^2)
    c_square = min(c_semi_minor_flesh^2, c_semi_minor^2)

    ssqg = (a_square * b_square * c_square) / (a_square * b_square + a_square * c_square + b_square * c_square)

    bs = b_semi_minor
    bl = b_semi_minor + insulation_depth
    br = bs + longwave_depth_fraction * insulation_depth

    if longwave_depth_fraction < 1
        T_ins1 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground - 
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins2 = ((3 * volume * bs) / ((((3 * ssqg)^0.5)^3) * dv1)) * 
         (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + 
         cd * T_substrate - Q_evap_insulation + Q_solar
        T_ins4 = (3 * volume * bs * cd3) / ((((3 * ssqg)^0.5)^3) * dv1) + 
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((((k_insulation * bl) / (bl - br)) * 
         (1 - conduction_fraction)) / dv4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(dv3 / dv4 + ((T_ins_compressed * cd2) / dv4) + 
         (T_insulation_calc * (((k_insulation * bl) / (bl - br)) * (1 - conduction_fraction))) / dv4)
    else
        T_ins1 = ((3 * volume * bs) / ((((3 * ssqg)^0.5)^3) * dv1)) * 
         (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_substrate - 
         Q_evap_insulation + Q_solar
        T_ins4 = (3 * volume * bs * cd3) / ((((3 * ssqg)^0.5)^3) * dv1) + 
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_insulation_calc
    end 

    return (; T_insulation_calc, T_radiant2)
end
