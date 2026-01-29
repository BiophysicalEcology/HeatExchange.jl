"""
    insulation_radiant_temperature(; body, insulation, insulation_pars, env_temps, Q_rads, cds, dvs, side, area_convection, hc, cd, k_insulation, longwave_depth_fraction, conduction_fraction, Q_solar, Q_evap_insulation, T_core, T_ins_compressed)

Calculate insulation surface temperature and radiant temperature for heat exchange.

Solves the heat balance at the insulation-air interface accounting for convection,
longwave radiation, solar radiation, and evaporation from the insulation surface.

# Keywords
- `body::AbstractBody`: Body geometry
- `insulation::InsulationProperties`: Computed insulation properties
- `insulation_pars::InsulationParameters`: Insulation parameters
- `env_temps::EnvironmentTemperatures`: Environmental temperatures
- `Q_rads::RadiationCoeffs`: Radiation coefficients to sky, bush, vegetation, ground
- `cds::ConductanceCoeffs`: Conductance coefficients
- `dvs::DivisorCoeffs`: Divisor coefficients
- `side`: Body side (`:dorsal` or `:ventral`)
- `area_convection`: Surface area for convection
- `hc`: Convective heat transfer coefficient
- `cd`: Substrate conductance coefficient
- `k_insulation`: Effective insulation thermal conductivity
- `longwave_depth_fraction`: Fraction of insulation depth for longwave exchange
- `conduction_fraction`: Fraction of body in contact with substrate
- `Q_solar`: Solar heat gain
- `Q_evap_insulation`: Evaporative heat loss from insulation
- `T_core`: Core body temperature
- `T_ins_compressed`: Compressed insulation temperature

# Returns
NamedTuple with:
- `T_insulation_calc`: Calculated insulation surface temperature
- `T_radiant2`: Radiant temperature for longwave exchange
"""
function insulation_radiant_temperature(;
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    env_temps::EnvironmentTemperatures,
    Q_rads::RadiationCoeffs,
    cds::ConductanceCoeffs,
    dvs::DivisorCoeffs,
    side,
    area_convection,
    hc,
    cd,
    k_insulation,
    longwave_depth_fraction,
    conduction_fraction,
    Q_solar,
    Q_evap_insulation,
    T_core,
    T_ins_compressed,
)
    insulation_radiant_temperature(
        shape(body),
        body,
        insulation,
        insulation_pars,
        env_temps,
        Q_rads,
        cds,
        dvs,
        side,
        area_convection,
        hc,
        cd,
        k_insulation,
        longwave_depth_fraction,
        conduction_fraction,
        Q_solar,
        Q_evap_insulation,
        T_core,
        T_ins_compressed,
    )
end
function insulation_radiant_temperature(
    shape::Union{Cylinder,Plate},
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    env_temps::EnvironmentTemperatures,
    Q_rads::RadiationCoeffs,
    cds::ConductanceCoeffs,
    dvs::DivisorCoeffs,
    side,
    area_convection,
    hc,
    cd,
    k_insulation,
    longwave_depth_fraction,
    conduction_fraction,
    Q_solar,
    Q_evap_insulation,
    T_core,
    T_ins_compressed,
)
    T = env_temps
    (; cd1, cd2, cd3) = cds
    (; dv1, dv2, dv3, dv4) = dvs
    (; Q_rad1, Q_rad2, Q_rad3, Q_rad4) = Q_rads

    r_skin = skin_radius(body)
    r_insulation = insulation_radius(body)
    length = body.geometry.length.length_skin

    if side == :dorsal
        r_radiation =
            r_skin +
            insulation_pars.longwave_depth_fraction *
            insulation_pars.dorsal.depth
    else
        r_radiation =
            r_skin +
            insulation_pars.longwave_depth_fraction *
            insulation_pars.ventral.depth
    end

    if longwave_depth_fraction < 1
        T_ins1 =
            Q_rad1 * T.sky + Q_rad2 * T.bush + Q_rad3 * T.vegetation + Q_rad4 * T.ground -
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) *
            ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins2 = ((2 * π * LEN) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 =
            hc * area_convection * T.air - cd * T_ins_compressed + cd * T.substrate -
            Q_evap_insulation + Q_solar
        T_ins4 =
            (2 * π * length * cd3) / dv1 +
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * (
                (
                    (k_insulation / log(r_insulation / r_radiation)) *
                    (1 - conduction_fraction)
                ) / dv4
            )
        + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(
            dv3 / dv4 +
            ((T_ins_compressed * cd2) / dv4) +
            (
                T_insulation_calc * (
                    (k_insulation / log(r_insulation / r_radiation)) *
                    (1 - conduction_fraction)
                )
            ) / dv4,
        )
    else
        T_ins1 =
            Q_rad1 * T.sky + Q_rad2 * T.bush + Q_rad3 * T.vegetation + Q_rad4 * T.ground
        T_ins2 = ((2 * π * length) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 =
            hc * area_convection * T.air - cd * T_ins_compressed + cd * T.substrate -
            Q_evap_insulation + Q_solar
        T_ins4 =
            (2 * π * length * cd3) / dv1 +
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) +
            hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_insulation_calc
    end

    return (; T_insulation_calc, T_radiant2)
end
function insulation_radiant_temperature(
    shape::Sphere,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    env_temps::EnvironmentTemperatures,
    Q_rads::RadiationCoeffs,
    cds::ConductanceCoeffs,
    dvs::DivisorCoeffs,
    side,
    area_convection,
    hc,
    cd,
    k_insulation,
    longwave_depth_fraction,
    conduction_fraction,
    Q_solar,
    Q_evap_insulation,
    T_core,
    T_ins_compressed,
)
    T = env_temps
    (; cd1, cd2, cd3) = cds
    (; dv1, dv2, dv3, dv4) = dvs
    (; Q_rad1, Q_rad2, Q_rad3, Q_rad4) = Q_rads

    r_skin = skin_radius(body)
    r_insulation = insulation_radius(body)

    if side == :dorsal
        r_radiation =
            r_skin +
            insulation_pars.longwave_depth_fraction *
            insulation_pars.dorsal.depth
    else
        r_radiation =
            r_skin +
            insulation_pars.longwave_depth_fraction *
            insulation_pars.ventral.depth
    end

    if longwave_depth_fraction < 1
        T_ins1 = ((4 * π * r_skin) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 =
            Q_rad1 * T.sky + Q_rad2 * T.bush + Q_rad3 * T.vegetation + Q_rad4 * T.ground -
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) *
            ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins3 =
            hc * area_convection * T.air - cd * T_ins_compressed + cd * T.substrate -
            Q_evap_insulation + Q_solar
        T_ins4 =
            (4 * π * r_skin * cd3) / dv1 +
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * (
                (
                    ((k_insulation * r_insulation) / (r_insulation - r_radiation)) *
                    (1 - conduction_fraction)
                ) / dv4
            ) +
            hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(
            dv3 / dv4 +
            ((T_ins_compressed * cd2) / dv4) +
            (
                T_insulation_calc * (
                    ((k_insulation * r_insulation) / (r_insulation - r_radiation)) *
                    (1 - conduction_fraction)
                )
            ) / dv4,
        )
    else
        T_ins1 = ((4 * π * r_skin) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 =
            Q_rad1 * T.sky + Q_rad2 * T.bush + Q_rad3 * T.vegetation + Q_rad4 * T.ground
        T_ins3 =
            hc * area_convection * T.air - cd * T_ins_compressed + cd * T.substrate -
            Q_evap_insulation + Q_solar
        T_ins4 =
            (4 * π * r_skin * cd3) / dv1 +
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) +
            hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_insulation_calc
    end
    return (; T_insulation_calc, T_radiant2)
end

function insulation_radiant_temperature(
    shape::Ellipsoid,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    env_temps::EnvironmentTemperatures,
    Q_rads::RadiationCoeffs,
    cds::ConductanceCoeffs,
    dvs::DivisorCoeffs,
    side,
    area_convection,
    hc,
    cd,
    k_insulation,
    longwave_depth_fraction,
    conduction_fraction,
    Q_solar,
    Q_evap_insulation,
    T_core,
    T_ins_compressed,
)
    T = env_temps
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

    if side == :dorsal
        insulation_depth = insulation_pars.dorsal.depth
    else
        insulation_depth = insulation_pars.ventral.depth
    end

    a_square = min(a_semi_major_flesh^2, a_semi_major^2)
    b_square = min(b_semi_minor_flesh^2, b_semi_minor^2)
    c_square = min(c_semi_minor_flesh^2, c_semi_minor^2)

    ssqg =
        (a_square * b_square * c_square) /
        (a_square * b_square + a_square * c_square + b_square * c_square)

    bs = b_semi_minor
    bl = b_semi_minor + insulation_depth
    br = bs + longwave_depth_fraction * insulation_depth

    if longwave_depth_fraction < 1
        T_ins1 =
            Q_rad1 * T.sky + Q_rad2 * T.bush + Q_rad3 * T.vegetation + Q_rad4 * T.ground -
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) *
            ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins2 =
            ((3 * volume * bs) / ((((3 * ssqg)^0.5)^3) * dv1)) *
            (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 =
            hc * area_convection * T.air - cd * T_ins_compressed + cd * T.substrate -
            Q_evap_insulation + Q_solar
        T_ins4 =
            (3 * volume * bs * cd3) / ((((3 * ssqg)^0.5)^3) * dv1) +
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) *
            ((((k_insulation * bl) / (bl - br)) * (1 - conduction_fraction)) / dv4) +
            hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(
            dv3 / dv4 +
            ((T_ins_compressed * cd2) / dv4) +
            (
                T_insulation_calc *
                (((k_insulation * bl) / (bl - br)) * (1 - conduction_fraction))
            ) / dv4,
        )
    else
        T_ins1 =
            ((3 * volume * bs) / ((((3 * ssqg)^0.5)^3) * dv1)) *
            (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 =
            Q_rad1 * T.sky + Q_rad2 * T.bush + Q_rad3 * T.vegetation + Q_rad4 * T.ground
        T_ins3 =
            hc * area_convection * T.air - cd * T_ins_compressed + cd * T.substrate -
            Q_evap_insulation + Q_solar
        T_ins4 =
            (3 * volume * bs * cd3) / ((((3 * ssqg)^0.5)^3) * dv1) +
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) +
            hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_insulation_calc
    end

    return (; T_insulation_calc, T_radiant2)
end
