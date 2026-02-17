"""
    mean_skin_temperature(; body, insulation, insulation_pars, ks, cds, conduction_fraction, environment_flux, skin_evaporation_flux, core_temperature, calculated_insulation_temperature, compressed_insulation_temperature)

Calculate the mean skin temperature for an endotherm given heat fluxes and body geometry.

Dispatches on body shape (Cylinder, Plate, Sphere, Ellipsoid) to use shape-specific
heat conduction equations through flesh and fat layers.

# Keywords
- `body::AbstractBody`: Body geometry
- `insulation::InsulationProperties`: Computed insulation properties
- `insulation_pars::InsulationParameters`: Insulation parameters
- `ks::ThermalConductivities`: Thermal conductivities (flesh, fat, insulation)
- `cds::ConductanceCoeffs`: Conductance coefficients
- `conduction_fraction`: Fraction of body in contact with substrate
- `environment_flux`: Environmental heat flux
- `skin_evaporation_flux`: Evaporative heat loss from skin
- `core_temperature`: Core body temperature
- `calculated_insulation_temperature`: Calculated insulation surface temperature
- `compressed_insulation_temperature`: Compressed insulation temperature

# Returns
NamedTuple with:
- `mean_skin_temperature`: Mean skin temperature
- `skin_temperature_calc1`: Skin temperature from core-to-skin heat flow
"""
function mean_skin_temperature(;
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    ks::ThermalConductivities,
    cds::ConductanceCoeffs,
    conduction_fraction,
    environment_flux,
    skin_evaporation_flux,
    core_temperature,
    calculated_insulation_temperature,
    compressed_insulation_temperature,
)
    mean_skin_temperature(
        shape(body),
        body,
        insulation,
        insulation_pars,
        ks,
        cds,
        conduction_fraction,
        environment_flux,
        skin_evaporation_flux,
        core_temperature,
        calculated_insulation_temperature,
        compressed_insulation_temperature,
    )
end
function mean_skin_temperature(
    shape::Union{Cylinder,Plate},
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    ks::ThermalConductivities,
    cds::ConductanceCoeffs,
    conduction_fraction,
    environment_flux,
    skin_evaporation_flux,
    core_temperature,
    calculated_insulation_temperature,
    compressed_insulation_temperature,
)
    (; flesh=k_flesh, fat=k_fat) = ks
    (; cd1, cd2, cd3) = cds
    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    k_compressed = insulation.conductivity_compressed
    r_compressed = r_skin + insulation_pars.depth_compressed
    if conduction_fraction < 1
        skin_temperature_calc1 =
            core_temperature - (((environment_flux + skin_evaporation_flux) * r_flesh^2) / (4 * k_flesh * volume)) -
            (((environment_flux + skin_evaporation_flux) * r_flesh^2) / (2 * k_fat * volume)) *
            log(r_skin / r_flesh)
        skin_temperature_calc2 =
            ((environment_flux * r_flesh^2) / (2 * cd1 * volume)) +
            ((compressed_insulation_temperature * cd2) / cd1) +
            ((calculated_insulation_temperature * cd3) / cd1)
    else
        skin_temperature_calc1 =
            core_temperature - ((environment_flux * r_flesh^2) / (4 * k_flesh * volume)) -
            ((environment_flux * r_flesh^2.0) / (2 * k_fat * volume)) * log(r_skin / r_flesh)
        skin_temperature_calc2 =
            (
                ((environment_flux * r_flesh^2) / (2 * k_compressed * volume)) *
                log(r_compressed / r_skin)
            ) + compressed_insulation_temperature
    end
    mean_skin_temperature = (skin_temperature_calc1 + skin_temperature_calc2) / 2
    return (; mean_skin_temperature, skin_temperature_calc1)
end
function mean_skin_temperature(
    shape::Sphere,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    ks::ThermalConductivities,
    cds::ConductanceCoeffs,
    conduction_fraction,
    environment_flux,
    skin_evaporation_flux,
    core_temperature,
    calculated_insulation_temperature,
    compressed_insulation_temperature,
)
    (; flesh=k_flesh, fat=k_fat) = ks
    (; cd1, cd2, cd3) = cds
    k_compressed = insulation.conductivity_compressed
    volume = flesh_volume(body)
    r_skin = skin_radius(body)
    r_flesh = flesh_radius(body)
    r_compressed = r_skin + insulation_pars.depth_compressed
    if conduction_fraction < 1
        skin_temperature_calc1 =
            core_temperature - (((environment_flux + skin_evaporation_flux) * r_flesh^2) / (6 * k_flesh * volume)) -
            (((environment_flux + skin_evaporation_flux) * r_flesh^3.0) / (3 * k_fat * volume)) *
            ((r_skin - r_flesh) / (r_skin * r_flesh))
        skin_temperature_calc2 =
            ((environment_flux * r_flesh^3) / (3 * cd1 * volume * r_skin)) +
            ((compressed_insulation_temperature * cd2) / cd1) +
            ((calculated_insulation_temperature * cd3) / cd1)
    else
        skin_temperature_calc1 =
            core_temperature - ((environment_flux * r_flesh^2) / (6 * k_flesh * volume)) -
            ((environment_flux * r_flesh^3) / (3 * k_fat * volume)) *
            ((r_skin - r_flesh) / (r_skin * r_flesh))
        skin_temperature_calc2 =
            (
                ((environment_flux * r_flesh^3) / (3 * k_compressed * volume)) *
                ((r_compressed - r_skin) / (r_compressed * r_skin))
            ) + compressed_insulation_temperature
    end
    mean_skin_temperature = (skin_temperature_calc1 + skin_temperature_calc2) / 2
    return (; mean_skin_temperature, skin_temperature_calc1)
end
function mean_skin_temperature(
    shape::Ellipsoid,
    body::AbstractBody,
    insulation::InsulationProperties,
    insulation_pars::InsulationParameters,
    ks::ThermalConductivities,
    cds::ConductanceCoeffs,
    conduction_fraction,
    environment_flux,
    skin_evaporation_flux,
    core_temperature,
    calculated_insulation_temperature,
    compressed_insulation_temperature,
)
    (; flesh=k_flesh, fat=k_fat) = ks
    (; cd1, cd2, cd3) = cds
    volume = flesh_volume(body)
    k_compressed = insulation.conductivity_compressed
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
    bl_compressed = b_semi_minor + insulation_pars.depth_compressed
    bg = min(b_semi_minor, b_semi_minor_flesh)

    if conduction_fraction < 1
        skin_temperature_calc1 =
            core_temperature - (((environment_flux + skin_evaporation_flux) * ssqg) / (2 * k_flesh * volume)) -
            (((environment_flux + skin_evaporation_flux) * (((3 * ssqg)^0.5)^3)) / (3 * k_fat * volume)) *
            ((bs - bg) / (bs * bg))
        skin_temperature_calc2 =
            ((environment_flux * (((3 * ssqg)^0.5)^3)) / (3 * cd1 * volume * bs)) +
            ((compressed_insulation_temperature * cd2) / cd1) +
            ((calculated_insulation_temperature * cd3) / cd1)
    else
        skin_temperature_calc1 =
            core_temperature - ((environment_flux * ssqg) / (2 * k_flesh * volume)) -
            ((environment_flux * (((3 * ssqg)^0.5)^3)) / (3 * k_fat * volume)) *
            ((bs - bg) / (bs * bg))
        skin_temperature_calc2 =
            (
                ((environment_flux * (((3 * ssqg)^0.5)^3)) / (3 * k_compressed * volume)) *
                ((bl_compressed - bs) / (bl_compressed * bs))
            ) + compressed_insulation_temperature
    end
    mean_skin_temperature = (skin_temperature_calc1 + skin_temperature_calc2) / 2
    return (; mean_skin_temperature, skin_temperature_calc1)
end
