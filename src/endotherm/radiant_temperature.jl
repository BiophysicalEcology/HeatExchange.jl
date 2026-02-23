"""
    radiant_temperature(; body, insulation, insulation_pars, org_temps, conductivities, side, substrate_conductance, longwave_depth_fraction, conduction_fraction, evaporation_flux, substrate_temperature)

Calculate the radiant temperature at the insulation surface for longwave radiation exchange.

Dispatches on body shape (Cylinder, Plate, Sphere, Ellipsoid) to use shape-specific
heat conduction equations through the insulation layer.

# Keywords
- `body::AbstractBody`: Body geometry
- `insulation::InsulationProperties`: Computed insulation properties
- `insulation_pars::InsulationParameters`: Insulation parameters
- `org_temps::OrganismTemperatures`: Organism temperatures (core, skin, insulation)
- `conductivities::ThermalConductivities`: Thermal conductivities (flesh, fat, insulation)
- `side`: Body side (`:dorsal` or `:ventral`)
- `substrate_conductance`: Substrate conductance coefficient
- `longwave_depth_fraction`: Fraction of insulation depth for longwave exchange
- `conduction_fraction`: Fraction of body in contact with substrate
- `evaporation_flux`: Evaporative heat loss
- `substrate_temperature`: Substrate temperature

# Returns
NamedTuple with:
- `radiant_temperature`: Radiant temperature for longwave exchange
- `compressed_insulation_temperature`: Compressed insulation temperature
- `conductances::ConductanceCoeffs`: Conductance coefficients
- `divisors::DivisorCoeffs`: Divisor coefficients for calculations
"""
function radiant_temperature(;
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    org_temps::OrganismTemperatures,
    conductivities::ThermalConductivities,
    side,
    substrate_conductance,
    longwave_depth_fraction,
    conduction_fraction,
    evaporation_flux,
    substrate_temperature,
)
    radiant_temperature(
        shape(body),
        body,
        insulation,
        insulation_pars,
        org_temps,
        conductivities,
        side,
        substrate_conductance,
        longwave_depth_fraction,
        conduction_fraction,
        evaporation_flux,
        substrate_temperature,
    )
end
function radiant_temperature(
    shape::Union{Cylinder,Plate},
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    org_temps::OrganismTemperatures,
    conductivities::ThermalConductivities,
    side,
    substrate_conductance,
    longwave_depth_fraction,
    conduction_fraction,
    evaporation_flux,
    substrate_temperature,
)
    (; core_temperature, skin_temperature, insulation_temperature) = org_temps

    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    r_insulation = insulation_radius(body)
    insulation_depth = getproperty(insulation.fibres, side).depth
    r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation_depth
    compressed_conductivity = insulation.conductivity_compressed
    r_compressed = r_skin + insulation_pars.depth_compressed
    length = body.geometry.length.length_skin

    compression_fraction =
        (conduction_fraction * 2 * π * compressed_conductivity * length) / log(r_compressed / r_skin)
    compressed_insulation_temperature = if conduction_fraction > 0
        (compression_fraction * skin_temperature + substrate_conductance * substrate_temperature) / (substrate_conductance + compression_fraction)
    else
        0.0u"K"
    end

    total_conductance =
        (compressed_conductivity / log(r_compressed / r_skin)) * conduction_fraction +
        (conductivities.insulation / log(r_insulation / r_skin)) * (1 - conduction_fraction)
    compressed_conductance = (compressed_conductivity / log(r_compressed / r_skin)) * conduction_fraction
    uncompressed_conductance = (conductivities.insulation / log(r_insulation / r_skin)) * (1 - conduction_fraction)

    geometric_divisor =
        1 +
        ((2 * π * length * r_flesh^2 * total_conductance) / (4 * conductivities.flesh * volume)) +
        ((2 * π * length * r_flesh^2 * total_conductance) / (2 * conductivities.fat * volume)) * log(r_skin / r_flesh)

    evaporative_divisor =
        evaporation_flux * ((r_flesh^2 * total_conductance) / (4 * conductivities.flesh * volume)) +
        evaporation_flux * ((r_flesh^2 * total_conductance) / (2 * conductivities.fat * volume)) * log(r_skin / r_flesh)

    numerator_divisor =
        ((2 * π * length) / geometric_divisor) *
        (core_temperature * total_conductance - evaporative_divisor - compressed_insulation_temperature * compressed_conductance - insulation_temperature * uncompressed_conductance) *
        r_flesh^2 / (2 * volume)

    radiative_divisor = if longwave_depth_fraction < 1
        compressed_conductance + (conductivities.insulation / log(r_insulation / r_radiation)) * (1 - conduction_fraction)
    else
        1.0
    end

    radiant_temperature = if longwave_depth_fraction < 1
        numerator_divisor / radiative_divisor +
        (compressed_insulation_temperature * compressed_conductance) / radiative_divisor +
        (
            insulation_temperature * (
                (conductivities.insulation / log(r_insulation / r_radiation)) *
                (1 - conduction_fraction)
            )
        ) / radiative_divisor
    else
        insulation_temperature
    end
    conductances = ConductanceCoeffs(total_conductance, compressed_conductance, uncompressed_conductance)
    divisors = DivisorCoeffs(geometric_divisor, evaporative_divisor, numerator_divisor, radiative_divisor)
    return (; radiant_temperature, compressed_insulation_temperature, conductances, divisors)
end
function radiant_temperature(
    shape::Sphere,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    org_temps::OrganismTemperatures,
    conductivities::ThermalConductivities,
    side,
    substrate_conductance,
    longwave_depth_fraction,
    conduction_fraction,
    evaporation_flux,
    substrate_temperature,
)
    (; core_temperature, skin_temperature, insulation_temperature) = org_temps

    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    r_insulation = insulation_radius(body)
    insulation_depth = getproperty(insulation.fibres, side).depth
    r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation_depth
    compressed_conductivity = insulation.conductivity_compressed
    r_compressed = r_skin + insulation_pars.depth_compressed

    compression_fraction =
        (conduction_fraction * 4 * π * compressed_conductivity * r_compressed * r_skin) /
        (r_compressed - r_skin)

    compressed_insulation_temperature = if conduction_fraction > 0
        (compression_fraction * skin_temperature + substrate_conductance * substrate_temperature) / (substrate_conductance + compression_fraction)
    else
        0.0u"K"
    end

    total_conductance =
        ((compressed_conductivity * r_compressed) / (r_compressed - r_skin)) * conduction_fraction +
        ((conductivities.insulation * r_insulation) / (r_insulation - r_skin)) *
        (1 - conduction_fraction)

    compressed_conductance = ((compressed_conductivity * r_compressed) / (r_compressed - r_skin)) * conduction_fraction
    uncompressed_conductance =
        ((conductivities.insulation * r_insulation) / (r_insulation - r_skin)) *
        (1 - conduction_fraction)

    geometric_divisor =
        1 +
        ((4 * π * r_skin * r_flesh^2 * total_conductance) / (6 * conductivities.flesh * volume)) +
        ((4 * π * r_skin * r_flesh^3 * total_conductance) / (3 * conductivities.fat * volume)) *
        ((r_skin - r_flesh) / (r_flesh * r_skin))

    evaporative_divisor =
        evaporation_flux * ((r_flesh^2 * total_conductance) / (6 * conductivities.flesh * volume)) +
        evaporation_flux *
        ((r_flesh^3 * total_conductance) / (3 * conductivities.fat * volume)) *
        ((r_skin - r_flesh) / (r_flesh * r_skin))

    numerator_divisor =
        ((4 * π * r_skin) / geometric_divisor) *
        (core_temperature * total_conductance - evaporative_divisor - compressed_insulation_temperature * compressed_conductance - insulation_temperature * uncompressed_conductance) *
        r_flesh^3 / (3 * volume * r_radiation)

    radiative_divisor = if longwave_depth_fraction < 1
        compressed_conductance +
        ((conductivities.insulation * r_insulation) / (r_insulation - r_radiation)) *
        (1 - conduction_fraction)
    else
        1.0
    end

    radiant_temperature = if longwave_depth_fraction < 1
        numerator_divisor / radiative_divisor +
        (compressed_insulation_temperature * compressed_conductance) / radiative_divisor +
        (
            insulation_temperature * (
                (conductivities.insulation * r_insulation) / (r_insulation - r_radiation) *
                (1 - conduction_fraction)
            )
        ) / radiative_divisor
    else
        insulation_temperature
    end
    conductances = ConductanceCoeffs(total_conductance, compressed_conductance, uncompressed_conductance)
    divisors = DivisorCoeffs(geometric_divisor, evaporative_divisor, numerator_divisor, radiative_divisor)
    return (; radiant_temperature, compressed_insulation_temperature, conductances, divisors)
end
function radiant_temperature(
    shape::Ellipsoid,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    org_temps::OrganismTemperatures,
    conductivities::ThermalConductivities,
    side,
    substrate_conductance,
    longwave_depth_fraction,
    conduction_fraction,
    evaporation_flux,
    substrate_temperature,
)
    (; core_temperature, skin_temperature, insulation_temperature) = org_temps

    volume = flesh_volume(body)

    a_semi_major = body.geometry.length.a_semi_major_skin
    b_semi_minor = body.geometry.length.b_semi_minor_skin
    c_semi_minor = body.geometry.length.c_semi_minor_skin
    fat = body.geometry.length.fat
    a_semi_major_flesh = a_semi_major - fat
    b_semi_minor_flesh = b_semi_minor - fat
    c_semi_minor_flesh = c_semi_minor - fat

    insulation_depth = getproperty(insulation.fibres, side).depth
    compressed_conductivity = insulation.conductivity_compressed
    bl_compressed = b_semi_minor + insulation_pars.depth_compressed

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
        (conduction_fraction * 3 * compressed_conductivity * volume * bl_compressed * bs) /
        ((sqrt(3 * ssqg))^3 * (bl_compressed - bs))

    compressed_insulation_temperature = if conduction_fraction > 0.0
        (compression_fraction * skin_temperature + substrate_conductance * substrate_temperature) / (substrate_conductance + compression_fraction)
    else
        0.0u"K"
    end

    total_conductance =
        ((compressed_conductivity * bl_compressed) / (bl_compressed - bs)) * conduction_fraction +
        ((conductivities.insulation * bl) / (bl - bs)) * (1 - conduction_fraction)
    compressed_conductance = ((compressed_conductivity * bl_compressed) / (bl_compressed - bs)) * conduction_fraction
    uncompressed_conductance = ((conductivities.insulation * bl) / (bl - bs)) * (1 - conduction_fraction)
    geometric_divisor =
        1 +
        (3 * bs * ssqg * total_conductance) / (2 * conductivities.flesh * (sqrt(3 * ssqg)^3)) +
        (bs * total_conductance) / conductivities.fat * ((bs - bg) / (bs * bg))

    evaporative_divisor =
        evaporation_flux * ((ssqg * total_conductance) / (2 * conductivities.flesh * volume)) +
        evaporation_flux * ((sqrt(3 * ssqg)^3 * total_conductance) / (3 * conductivities.fat * volume)) * ((bs - bg) / (bs * bg))

    numerator_divisor =
        (bs / geometric_divisor) * (core_temperature * total_conductance - evaporative_divisor - compressed_insulation_temperature * compressed_conductance - insulation_temperature * uncompressed_conductance) / br

    radiative_divisor = if longwave_depth_fraction < 1
        compressed_conductance + ((conductivities.insulation * bl) / (bl - br)) * (1 - conduction_fraction)
    else
        1.0
    end

    radiant_temperature = if longwave_depth_fraction < 1
        numerator_divisor / radiative_divisor +
        (compressed_insulation_temperature * compressed_conductance) / radiative_divisor +
        (insulation_temperature * ((conductivities.insulation * bl) / (bl - br) * (1 - conduction_fraction))) / radiative_divisor
    else
        insulation_temperature
    end
    conductances = ConductanceCoeffs(total_conductance, compressed_conductance, uncompressed_conductance)
    divisors = DivisorCoeffs(geometric_divisor, evaporative_divisor, numerator_divisor, radiative_divisor)
    return (; radiant_temperature, compressed_insulation_temperature, conductances, divisors)
end
