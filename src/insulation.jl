"""
    insulation_thermal_conductivity(fibres, air_conductivity; depth=nothing)

Compute effective thermal conductivity and radiation parameters for insulation fibres.

Implements the method of Conley & Porter (1986) and Kowalski (1983) to calculate
effective thermal conductivity, absorption coefficient, and optical thickness
of an insulation layer (fur or feathers).

# Arguments
- `fibre::FibreProperties`: Fibre properties (diameter, length, density, depth, reflectance, conductivity).
- `air_conductivity`: Thermal conductivity of air (W/m/K).

# Keywords
- `depth`: Override `fibres.depth` with this value (for compressed insulation calculations).

# Returns
NamedTuple with:
- `effective_conductivity`: Effective insulation thermal conductivity (W/m/K).
- `absorption_coefficient`: Average absorption coefficient (m⁻¹).
- `optical_thickness_factor`: Optical thickness (dimensionless).

# References
- Conley, K. E., & Porter, W. P. (1986). Modeling fur insulation. *Journal of Thermal Biology*.
- Kowalski, K. A. (1983). *A model of heat transfer through mammalian fur*. PhD Thesis,
  University of Wisconsin.
"""
function insulation_thermal_conductivity(fibre::FibreProperties, air_conductivity; depth=nothing)
    # Allow depth overrides when compressed
    insulation_depth = isnothing(depth) ? fibre.depth : depth
    density = fibre.density

    w = 1.0  # weighting factor placeholder

    # Effective density (Conley & Porter 1986, p.254)
    effective_density = density * (fibre.length / insulation_depth)

    # Distance between fibre centers (assuming uniform spacing)
    l_unit = 1.0 / sqrt(effective_density)
    fibre_spacing = l_unit - fibre.diameter

    # Initialize variables
    a_air = r_air = r_fibre = a_fibre = kx = ky = 0.0

    # Check for feasible fibre density/diameter
    if fibre_spacing > 0.0u"m"
        a_air = (w / 2.0) * (l_unit - fibre.diameter)
        r_air = u"K*m*W^-1"(l_unit / (air_conductivity * a_air))
        r_fibre = u"K*m*W^-1"(
            (
                (fibre.diameter * air_conductivity) +
                (l_unit - fibre.diameter) * fibre.conductivity
            ) / (w * fibre.diameter * fibre.conductivity * air_conductivity),
        )
        a_fibre = density * ((u"cm"(fibre.diameter) / 2.0)^2 * π)
        kx = a_fibre * fibre.conductivity + (1.0 - a_fibre) * air_conductivity
        ky = (2.0 / r_air) + (1.0 / r_fibre)
    else
        # No space between fibres — recalculate geometry
        if fibre_spacing < 0.0u"m"
            l_unit = fibre.diameter
            effective_density = 1.0 / l_unit^2
            density = (effective_density * insulation_depth) / fibre.length
            fibre_spacing = l_unit - fibre.diameter
        end

        a_air = (w / 2.0) * (l_unit - fibre.diameter)
        r_air =
            2.0 /
            (sqrt(effective_density) * air_conductivity * (l_unit - fibre.diameter) * w)
        r_fibre =
            (
                fibre.diameter * air_conductivity +
                (l_unit - fibre.diameter) * fibre.conductivity
            ) / (w * fibre.diameter * fibre.conductivity * air_conductivity)
        a_fibre = effective_density * ((fibre.diameter / 2.0)^2 * π)
        kx = a_fibre * fibre.conductivity + (1.0 - a_fibre) * air_conductivity
        ky = (2.0 / r_air) + (1.0 / r_fibre)
    end

    # Effective thermal conductivity (eq. 3-28, Kowalski 1983)
    effective_conductivity = (ky + kx) / 2.0

    # Ensure air_conductivity < effective_conductivity < fibre.conductivity
    if effective_conductivity > fibre.conductivity
        effective_conductivity = fibre.conductivity
    elseif effective_conductivity < air_conductivity
        effective_conductivity = air_conductivity
    end

    # Absorption and optical parameters
    absorption_coefficient = (0.67 / π) * u"m^-2"(effective_density) * u"m"(fibre.diameter)
    optical_thickness_factor = absorption_coefficient * u"m"(insulation_depth)

    return (; effective_conductivity, absorption_coefficient, optical_thickness_factor)
end

"""
    insulation_properties(insulation, insulation_temperature, ventral_fraction) -> InsulationProperties

Compute parameters for heat conduction and infrared radiation through insulation (fur or plumage).

# Arguments
- `insulation::InsulationParameters`: Insulation parameters with dorsal/ventral `FibreProperties`.
- `insulation_temperature`: Temperature of insulation for computing air conductivity (K).
- `ventral_fraction`: Fraction of body surface that is ventral (0-1).

# Returns
`InsulationProperties` containing:
- `fibres`: `BodyRegionValues{FibreProperties}` for average/dorsal/ventral regions.
- `conductivities`: Effective thermal conductivities (W/m/K).
- `absorption_coefficients`: Absorption coefficients for radiation (m⁻¹).
- `optical_thickness`: Optical thickness factors (dimensionless).
- `insulation_test`: Bare-skin test parameter (zero if no insulation).
- `conductivity_compressed`: Conductivity of compressed ventral insulation (W/m/K).
"""
function insulation_properties(insulation::InsulationParameters, insulation_temperature, ventral_fraction)
    (; dorsal, ventral, depth_compressed) = insulation

    # Physical constants
    air_conductivity = dry_air_properties(insulation_temperature).thermal_conductivity

    # Bare-skin test (dorsal only)
    insulation_test = dorsal.density * dorsal.diameter * dorsal.length * dorsal.depth

    # Weighted average fibre properties
    avg_fibres = interpolate(dorsal, ventral, ventral_fraction)

    # Adjusted ventral values (accounting for partial insulation coverage)
    ventral_adj = interpolate(dorsal, ventral, min(ventral_fraction * 2.0, 1.0))

    fibres = BodyRegionValues(avg_fibres, dorsal, ventral_adj)

    # Compute insulation thermal parameters for each region
    conductivity_compressed = 0.0u"W/m/K"
    if insulation_test <= 0.0u"m"
        conductivities = BodyRegionValues(0.0u"W/m/K", 0.0u"W/m/K", 0.0u"W/m/K")
        absorption_coefficients = BodyRegionValues(0.0u"m^-1", 0.0u"m^-1", 0.0u"m^-1")
        optical_thickness = BodyRegionValues(0.0, 0.0, 0.0)
    else
        avg = insulation_thermal_conductivity(fibres.average, air_conductivity)
        dors = insulation_thermal_conductivity(fibres.dorsal, air_conductivity)
        vent = insulation_thermal_conductivity(fibres.ventral, air_conductivity)
        conductivities = BodyRegionValues(
            avg.effective_conductivity, dors.effective_conductivity, vent.effective_conductivity
        )
        absorption_coefficients = BodyRegionValues(
            avg.absorption_coefficient, dors.absorption_coefficient, vent.absorption_coefficient
        )
        optical_thickness = BodyRegionValues(
            avg.optical_thickness_factor, dors.optical_thickness_factor, vent.optical_thickness_factor
        )
        # Compressed ventral insulation conductivity
        vent_compressed = insulation_thermal_conductivity(fibres.ventral, air_conductivity; depth=depth_compressed)
        conductivity_compressed = vent_compressed.effective_conductivity
    end

    return InsulationProperties(
        fibres,
        conductivities,
        absorption_coefficients,
        optical_thickness,
        insulation_test,
        conductivity_compressed,
    )
end
