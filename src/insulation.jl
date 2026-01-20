"""
    insulation_thermal_conductivity(fibre_density, fibre_length, insulation_depth, 
    fibre_diameter, air_conductivity, fibre_conductivity)

Compute feasible insulation fibre (hair/feather) spacing and parameters needed for 
conduction and infrared radiation through insulation.

This function implements the method described by Conley & Porter (1986) and Kowalski (1983).
It calculates the effective thermal conductivity, absorption coefficient, and optical thickness
of an insulation layer (fur or feathers) based on geometric and physical fibre/fiber parameters.

# Arguments
- `fibre_density::Quantity`: fibre density (fibres per m²)
- `fibre_length::Quantity`: hair (or feather) length (m)
- `insulation_depth::Quantity`: insulation depth (m)
- `fibre_diameter::Quantity`: fibre diameter (m)
- `air_conductivity::Quantity`: thermal conductivity of air (W m⁻¹ K⁻¹)
- `fibre_conductivity::Quantity`: thermal conductivity of fibre (W m⁻¹ K⁻¹)

# Returns
`Vector{Quantity}` of length 3:
1. `effective_conductivity` — effective insulation thermal conductivity (W m⁻¹ K⁻¹)
2. `absorption_coefficient` — average absorption coefficient (m⁻¹)
3. `optical_thickness_factor` — optical thickness (dimensionless)

# References
- Conley, K. E., & Porter, W. P. (1986). Modeling fur insulation. *Journal of Thermal Biology*.
- Kowalski, K. A. (1983). *A model of heat transfer through mammalian fur*. PhD Thesis,
     University of Wisconsin.
"""
function insulation_thermal_conductivity(;
    fibre_density::Quantity,
    fibre_length::Quantity,
    insulation_depth::Quantity,
    fibre_diameter::Quantity,
    air_conductivity::Quantity,
    fibre_conductivity::Quantity,
)
    w = 1.0  # unused weighting factor

    # Effective density (Conley & Porter 1986, p.254)
    effective_density = fibre_density * (fibre_length / insulation_depth)

    # Distance between fibre centers (assuming uniform spacing)
    l_unit = 1.0 / sqrt(effective_density)
    fibre_spacing = l_unit - fibre_diameter

    # Initialize variables
    a_air = r_air = r_fibre = a_fibre = kx = ky = 0.0

    # Check for feasible fibre density/diameter
    if fibre_spacing > 0.0u"m"
        a_air = (w / 2.0) * (l_unit - fibre_diameter)
        r_air = u"K*m*W^-1"(l_unit / (air_conductivity * a_air))
        r_fibre = u"K*m*W^-1"(
            (
                (fibre_diameter * air_conductivity) +
                (l_unit - fibre_diameter) * fibre_conductivity
            ) / (w * fibre_diameter * fibre_conductivity * air_conductivity),
        )
        a_fibre = fibre_density * ((u"cm"(fibre_diameter) / 2.0)^2 * π)
        kx = a_fibre * fibre_conductivity + (1.0 - a_fibre) * air_conductivity
        ky = (2.0 / r_air) + (1.0 / r_fibre)
    else
        # No space between fibres — recalculate density
        if fibre_spacing < 0.0u"m"
            l_unit = fibre_diameter
            effective_density = 1.0 / l_unit^2
            fibre_density = (effective_density * insulation_depth) / fibre_length
        end

        a_air = (w / 2.0) * (l_unit - fibre_diameter)
        r_air = u"K*m*W^-1"(l_unit / (air_conductivity * a_air))
        r_fibre = u"K*m*W^-1"(
            (
                (u"cm"(fibre_diameter) * air_conductivity) +
                (l_unit - fibre_diameter) * fibre_conductivity
            ) / (w * fibre_diameter * fibre_conductivity * air_conductivity),
        )
        a_fibre = effective_density * ((fibre_diameter / 2.0)^2 * π)
        kx = a_fibre * fibre_conductivity + (1.0 - a_fibre) * air_conductivity
        ky = (2.0 / r_air) + (1.0 / r_fibre)
        fibre_spacing = l_unit - fibre_diameter
        r_air =
            2.0 /
            (sqrt(effective_density) * air_conductivity * (l_unit - fibre_diameter) * w)
        r_fibre =
            (
                fibre_diameter * air_conductivity +
                (l_unit - fibre_diameter) * fibre_conductivity
            ) / (w * fibre_diameter * fibre_conductivity * air_conductivity)
        a_air = (w / 2.0) * (l_unit - fibre_diameter)
        ky = (2.0 / r_air) + (1.0 / r_fibre)
    end

    # Effective thermal conductivity (eq. 3-28, Kowalski 1983)
    effective_conductivity = (ky + kx) / 2.0

    # Ensure air_conductivity < effective_conductivity < fibre_conductivity
    if effective_conductivity > fibre_conductivity
        effective_conductivity = fibre_conductivity
    elseif effective_conductivity < air_conductivity
        effective_conductivity = air_conductivity
    end

    # Absorption and optical parameters
    absorption_coefficient = (0.67 / π) * u"m^-2"(effective_density) * u"m"(fibre_diameter)
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
    (; dorsal, ventral, depth_compressed, fibre_conductivity) = insulation

    # Physical constants
    air_conductivity = dry_air_properties(insulation_temperature).thermal_conductivity

    # Bare-skin test (dorsal only)
    insulation_test = dorsal.density * dorsal.diameter * dorsal.length * dorsal.depth

    # Weighted average fibre properties
    avg_fibres = FibreProperties(;
        diameter=dorsal.diameter * (1 - ventral_fraction) + ventral.diameter * ventral_fraction,
        length=dorsal.length * (1 - ventral_fraction) + ventral.length * ventral_fraction,
        density=dorsal.density * (1 - ventral_fraction) + ventral.density * ventral_fraction,
        depth=dorsal.depth * (1 - ventral_fraction) + ventral.depth * ventral_fraction,
        reflectance=dorsal.reflectance * (1 - ventral_fraction) + ventral.reflectance * ventral_fraction,
    )

    # Adjusted ventral values (accounting for partial insulation coverage)
    pven_v = min(ventral_fraction * 2.0, 1.0)
    ventral_adj = FibreProperties(;
        diameter=dorsal.diameter * (1 - pven_v) + ventral.diameter * pven_v,
        length=dorsal.length * (1 - pven_v) + ventral.length * pven_v,
        density=dorsal.density * (1 - pven_v) + ventral.density * pven_v,
        depth=dorsal.depth * (1 - pven_v) + ventral.depth * pven_v,
        reflectance=dorsal.reflectance * (1 - pven_v) + ventral.reflectance * pven_v,
    )

    fibres = BodyRegionValues(avg_fibres, dorsal, ventral_adj)

    # Compute insulation thermal parameters for each region
    conductivity_compressed = 0.0u"W/m/K"
    if insulation_test <= 0.0u"m"
        conductivities = BodyRegionValues(0.0u"W/m/K", 0.0u"W/m/K", 0.0u"W/m/K")
        absorption_coefficients = BodyRegionValues(0.0u"m^-1", 0.0u"m^-1", 0.0u"m^-1")
        optical_thickness = BodyRegionValues(0.0, 0.0, 0.0)
    else
        avg = insulation_thermal_conductivity(;
            fibre_density=fibres.average.density,
            fibre_length=fibres.average.length,
            insulation_depth=fibres.average.depth,
            fibre_diameter=fibres.average.diameter,
            air_conductivity,
            fibre_conductivity,
        )
        dors = insulation_thermal_conductivity(;
            fibre_density=fibres.dorsal.density,
            fibre_length=fibres.dorsal.length,
            insulation_depth=fibres.dorsal.depth,
            fibre_diameter=fibres.dorsal.diameter,
            air_conductivity,
            fibre_conductivity,
        )
        vent = insulation_thermal_conductivity(;
            fibre_density=fibres.ventral.density,
            fibre_length=fibres.ventral.length,
            insulation_depth=fibres.ventral.depth,
            fibre_diameter=fibres.ventral.diameter,
            air_conductivity,
            fibre_conductivity,
        )
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
        vent_compressed = insulation_thermal_conductivity(;
            fibre_density=fibres.ventral.density,
            fibre_length=fibres.ventral.length,
            insulation_depth=depth_compressed,
            fibre_diameter=fibres.ventral.diameter,
            air_conductivity,
            fibre_conductivity,
        )
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
