"""
    insulation_radiant_temperature(; body, insulation, insulation_pars, env_temps, coeffs, conductances, divisors, side, area_convection, heat_transfer_coefficient, substrate_conductance, insulation_conductivity, longwave_depth_fraction, conduction_fraction, solar_flux, insulation_evaporation_flux, core_temperature, compressed_insulation_temperature)

Calculate insulation surface temperature and radiant temperature for heat exchange.

Solves the heat balance at the insulation-air interface accounting for convection, longwave radiation, solar radiation, and evaporation from the insulation surface.

# Keywords
- `body::AbstractBody`: Body geometry
- `insulation::InsulationProperties`: Computed insulation properties
- `insulation_pars::InsulationParameters`: Insulation parameters
- `env_temps::EnvironmentTemperatures`: Environmental temperatures
- `coeffs::RadiationCoeffs`: Radiation coefficients to sky, bush, vegetation, ground
- `conductances::ConductanceCoeffs`: Conductance coefficients
- `divisors::DivisorCoeffs`: Divisor coefficients
- `side`: Body side (`:dorsal` or `:ventral`)
- `area_convection`: Surface area for convection
- `heat_transfer_coefficient`: Convective heat transfer coefficient
- `substrate_conductance`: Substrate conductance coefficient
- `insulation_conductivity`: Effective insulation thermal conductivity
- `longwave_depth_fraction`: Fraction of insulation depth for longwave exchange
- `conduction_fraction`: Fraction of body in contact with substrate
- `solar_flux`: Solar heat gain
- `insulation_evaporation_flux`: Evaporative heat loss from insulation
- `core_temperature`: Core body temperature
- `compressed_insulation_temperature`: Compressed insulation temperature

# Returns
NamedTuple with:
- `calculated_insulation_temperature`: Calculated insulation surface temperature
- `radiant_temperature2`: Radiant temperature for longwave exchange
"""
function insulation_radiant_temperature(;
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    env_temps::EnvironmentTemperatures,
    coeffs::RadiationCoeffs,
    conductances::ConductanceCoeffs,
    divisors::DivisorCoeffs,
    side,
    area_convection,
    heat_transfer_coefficient,
    substrate_conductance,
    insulation_conductivity,
    longwave_depth_fraction,
    conduction_fraction,
    solar_flux,
    insulation_evaporation_flux,
    core_temperature,
    compressed_insulation_temperature,
)
    insulation_radiant_temperature(
        shape(body),
        body,
        insulation,
        insulation_pars,
        env_temps,
        coeffs,
        conductances,
        divisors,
        side,
        area_convection,
        heat_transfer_coefficient,
        substrate_conductance,
        insulation_conductivity,
        longwave_depth_fraction,
        conduction_fraction,
        solar_flux,
        insulation_evaporation_flux,
        core_temperature,
        compressed_insulation_temperature,
    )
end
function insulation_radiant_temperature(
    shape::Union{Cylinder,Plate},
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    env_temps::EnvironmentTemperatures,
    coeffs::RadiationCoeffs,
    conductances::ConductanceCoeffs,
    divisors::DivisorCoeffs,
    side,
    area_convection,
    heat_transfer_coefficient,
    substrate_conductance,
    insulation_conductivity,
    longwave_depth_fraction,
    conduction_fraction,
    solar_flux,
    insulation_evaporation_flux,
    core_temperature,
    compressed_insulation_temperature,
)
    T = env_temps
    total_conductance = conductances.total
    compressed_conductance = conductances.compressed
    uncompressed_conductance = conductances.uncompressed
    geometric_divisor = divisors.geometric
    evaporative_divisor = divisors.evaporative
    numerator_divisor = divisors.numerator
    radiative_divisor = divisors.radiative
    
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
        ins_calc1 =
            coeffs.sky * env_temps.sky + coeffs.bush * env_temps.bush + coeffs.vegetation * env_temps.vegetation + coeffs.ground * env_temps.ground -
            (coeffs.sky + coeffs.bush + coeffs.vegetation + coeffs.ground) *
            ((numerator_divisor / radiative_divisor) + ((compressed_insulation_temperature * compressed_conductance) / radiative_divisor))
        ins_calc2 = ((2 * π * LEN) / geometric_divisor) * (core_temperature * total_conductance - evaporative_divisor - compressed_insulation_temperature * compressed_conductance)
        ins_calc3 =
            heat_transfer_coefficient * area_convection * env_temps.air - substrate_conductance * compressed_insulation_temperature + substrate_conductance * env_temps.substrate -
            insulation_evaporation_flux + solar_flux
        ins_calc4 =
            (2 * π * length * uncompressed_conductance) / geometric_divisor +
            (coeffs.sky + coeffs.bush + coeffs.vegetation + coeffs.ground) * (
                (
                    (insulation_conductivity / log(r_insulation / r_radiation)) *
                    (1 - conduction_fraction)
                ) / radiative_divisor
            )
        + heat_transfer_coefficient * area_convection
        calculated_insulation_temperature = u"K"((ins_calc1 + ins_calc2 + ins_calc3) / ins_calc4)
        radiant_temperature2 = u"K"(
            numerator_divisor / radiative_divisor +
            ((compressed_insulation_temperature * compressed_conductance) / radiative_divisor) +
            (
                calculated_insulation_temperature * (
                    (insulation_conductivity / log(r_insulation / r_radiation)) *
                    (1 - conduction_fraction)
                )
            ) / radiative_divisor,
        )
    else
        ins_calc1 =
            coeffs.sky * env_temps.sky + coeffs.bush * env_temps.bush + coeffs.vegetation * env_temps.vegetation + coeffs.ground * env_temps.ground
        ins_calc2 = ((2 * π * length) / geometric_divisor) * (core_temperature * total_conductance - evaporative_divisor - compressed_insulation_temperature * compressed_conductance)
        ins_calc3 =
            heat_transfer_coefficient * area_convection * env_temps.air - substrate_conductance * compressed_insulation_temperature + substrate_conductance * env_temps.substrate -
            insulation_evaporation_flux + solar_flux
        ins_calc4 =
            (2 * π * length * uncompressed_conductance) / geometric_divisor +
            (coeffs.sky + coeffs.bush + coeffs.vegetation + coeffs.ground) +
            heat_transfer_coefficient * area_convection
        calculated_insulation_temperature = u"K"((ins_calc1 + ins_calc2 + ins_calc3) / ins_calc4)
        radiant_temperature2 = calculated_insulation_temperature
    end

    return (; calculated_insulation_temperature, radiant_temperature2)
end
function insulation_radiant_temperature(
    shape::Sphere,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    env_temps::EnvironmentTemperatures,
    coeffs::RadiationCoeffs,
    conductances::ConductanceCoeffs,
    divisors::DivisorCoeffs,
    side,
    area_convection,
    heat_transfer_coefficient,
    substrate_conductance,
    insulation_conductivity,
    longwave_depth_fraction,
    conduction_fraction,
    solar_flux,
    insulation_evaporation_flux,
    core_temperature,
    compressed_insulation_temperature,
)
    T = env_temps
    total_conductance = conductances.total
    compressed_conductance = conductances.compressed
    uncompressed_conductance = conductances.uncompressed
    geometric_divisor = divisors.geometric
    evaporative_divisor = divisors.evaporative
    numerator_divisor = divisors.numerator
    radiative_divisor = divisors.radiative
    
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
        ins_calc1 = ((4 * π * r_skin) / geometric_divisor) * (core_temperature * total_conductance - evaporative_divisor - compressed_insulation_temperature * compressed_conductance)
        ins_calc2 =
            coeffs.sky * env_temps.sky + coeffs.bush * env_temps.bush + coeffs.vegetation * env_temps.vegetation + coeffs.ground * env_temps.ground -
            (coeffs.sky + coeffs.bush + coeffs.vegetation + coeffs.ground) *
            ((numerator_divisor / radiative_divisor) + ((compressed_insulation_temperature * compressed_conductance) / radiative_divisor))
        ins_calc3 =
            heat_transfer_coefficient * area_convection * env_temps.air - substrate_conductance * compressed_insulation_temperature + substrate_conductance * env_temps.substrate -
            insulation_evaporation_flux + solar_flux
        ins_calc4 =
            (4 * π * r_skin * uncompressed_conductance) / geometric_divisor +
            (coeffs.sky + coeffs.bush + coeffs.vegetation + coeffs.ground) * (
                (
                    ((insulation_conductivity * r_insulation) / (r_insulation - r_radiation)) *
                    (1 - conduction_fraction)
                ) / radiative_divisor
            ) +
            heat_transfer_coefficient * area_convection
        calculated_insulation_temperature = u"K"((ins_calc1 + ins_calc2 + ins_calc3) / ins_calc4)
        radiant_temperature2 = u"K"(
            numerator_divisor / radiative_divisor +
            ((compressed_insulation_temperature * compressed_conductance) / radiative_divisor) +
            (
                calculated_insulation_temperature * (
                    ((insulation_conductivity * r_insulation) / (r_insulation - r_radiation)) *
                    (1 - conduction_fraction)
                )
            ) / radiative_divisor,
        )
    else
        ins_calc1 = ((4 * π * r_skin) / geometric_divisor) * (core_temperature * total_conductance - evaporative_divisor - compressed_insulation_temperature * compressed_conductance)
        ins_calc2 =
            coeffs.sky * env_temps.sky + coeffs.bush * env_temps.bush + coeffs.vegetation * env_temps.vegetation + coeffs.ground * env_temps.ground
        ins_calc3 =
            heat_transfer_coefficient * area_convection * env_temps.air - substrate_conductance * compressed_insulation_temperature + substrate_conductance * env_temps.substrate -
            insulation_evaporation_flux + solar_flux
        ins_calc4 =
            (4 * π * r_skin * uncompressed_conductance) / geometric_divisor +
            (coeffs.sky + coeffs.bush + coeffs.vegetation + coeffs.ground) +
            heat_transfer_coefficient * area_convection
        calculated_insulation_temperature = u"K"((ins_calc1 + ins_calc2 + ins_calc3) / ins_calc4)
        radiant_temperature2 = calculated_insulation_temperature
    end
    return (; calculated_insulation_temperature, radiant_temperature2)
end

function insulation_radiant_temperature(
    shape::Ellipsoid,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    env_temps::EnvironmentTemperatures,
    coeffs::RadiationCoeffs,
    conductances::ConductanceCoeffs,
    divisors::DivisorCoeffs,
    side,
    area_convection,
    heat_transfer_coefficient,
    substrate_conductance,
    insulation_conductivity,
    longwave_depth_fraction,
    conduction_fraction,
    solar_flux,
    insulation_evaporation_flux,
    core_temperature,
    compressed_insulation_temperature,
)
    T = env_temps
    total_conductance = conductances.total
    compressed_conductance = conductances.compressed
    uncompressed_conductance = conductances.uncompressed
    geometric_divisor = divisors.geometric
    evaporative_divisor = divisors.evaporative
    numerator_divisor = divisors.numerator
    radiative_divisor = divisors.radiative
    
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
        ins_calc1 =
            coeffs.sky * env_temps.sky + coeffs.bush * env_temps.bush + coeffs.vegetation * env_temps.vegetation + coeffs.ground * env_temps.ground -
            (coeffs.sky + coeffs.bush + coeffs.vegetation + coeffs.ground) *
            ((numerator_divisor / radiative_divisor) + ((compressed_insulation_temperature * compressed_conductance) / radiative_divisor))
        ins_calc2 =
            ((3 * volume * bs) / ((((3 * ssqg)^0.5)^3) * geometric_divisor)) *
            (core_temperature * total_conductance - evaporative_divisor - compressed_insulation_temperature * compressed_conductance)
        ins_calc3 =
            heat_transfer_coefficient * area_convection * env_temps.air - substrate_conductance * compressed_insulation_temperature + substrate_conductance * env_temps.substrate -
            insulation_evaporation_flux + solar_flux
        ins_calc4 =
            (3 * volume * bs * uncompressed_conductance) / ((((3 * ssqg)^0.5)^3) * geometric_divisor) +
            (coeffs.sky + coeffs.bush + coeffs.vegetation + coeffs.ground) *
            ((((insulation_conductivity * bl) / (bl - br)) * (1 - conduction_fraction)) / radiative_divisor) +
            heat_transfer_coefficient * area_convection
        calculated_insulation_temperature = u"K"((ins_calc1 + ins_calc2 + ins_calc3) / ins_calc4)
        radiant_temperature2 = u"K"(
            numerator_divisor / radiative_divisor +
            ((compressed_insulation_temperature * compressed_conductance) / radiative_divisor) +
            (
                calculated_insulation_temperature *
                (((insulation_conductivity * bl) / (bl - br)) * (1 - conduction_fraction))
            ) / radiative_divisor,
        )
    else
        ins_calc1 =
            ((3 * volume * bs) / ((((3 * ssqg)^0.5)^3) * geometric_divisor)) *
            (core_temperature * total_conductance - evaporative_divisor - compressed_insulation_temperature * compressed_conductance)
        ins_calc2 =
            coeffs.sky * env_temps.sky + coeffs.bush * env_temps.bush + coeffs.vegetation * env_temps.vegetation + coeffs.ground * env_temps.ground
        ins_calc3 =
            heat_transfer_coefficient * area_convection * env_temps.air - substrate_conductance * compressed_insulation_temperature + substrate_conductance * env_temps.substrate -
            insulation_evaporation_flux + solar_flux
        ins_calc4 =
            (3 * volume * bs * uncompressed_conductance) / ((((3 * ssqg)^0.5)^3) * geometric_divisor) +
            (coeffs.sky + coeffs.bush + coeffs.vegetation + coeffs.ground) +
            heat_transfer_coefficient * area_convection
        calculated_insulation_temperature = u"K"((ins_calc1 + ins_calc2 + ins_calc3) / ins_calc4)
        radiant_temperature2 = calculated_insulation_temperature
    end

    return (; calculated_insulation_temperature, radiant_temperature2)
end
