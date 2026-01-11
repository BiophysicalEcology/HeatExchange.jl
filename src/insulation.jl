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
function insulation_thermal_conductivity(; fibre_density::Quantity, fibre_length::Quantity, 
    insulation_depth::Quantity, fibre_diameter::Quantity, air_conductivity::Quantity, 
    fibre_conductivity::Quantity)

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
        r_fibre = u"K*m*W^-1"(((fibre_diameter * air_conductivity) + 
            (l_unit - fibre_diameter) * fibre_conductivity) / 
                (w * fibre_diameter * fibre_conductivity * air_conductivity))
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
        r_fibre = u"K*m*W^-1"(((u"cm"(fibre_diameter) * air_conductivity) 
            + (l_unit - fibre_diameter) * fibre_conductivity) / 
            (w * fibre_diameter * fibre_conductivity * air_conductivity))
        a_fibre = effective_density * ((fibre_diameter / 2.0)^2 * π)
        kx = a_fibre * fibre_conductivity + (1.0 - a_fibre) * air_conductivity
        ky = (2.0 / r_air) + (1.0 / r_fibre)
        fibre_spacing = l_unit - fibre_diameter
        r_air = 2.0 / (sqrt(effective_density) * air_conductivity * 
            (l_unit - fibre_diameter) * w)
        r_fibre = (fibre_diameter * air_conductivity + 
            (l_unit - fibre_diameter) * fibre_conductivity) / 
            (w * fibre_diameter * fibre_conductivity * air_conductivity)
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
    insulation_properties(insulation)

Compute parameters for heat conduction and infrared radiation through insulation (fur or plumage).

# Arguments
`insulation`, an InsulationPars struct containing
- `fibre_diameter_dorsal`, `fibre_diameter_ventral` : fibre diameters
- `fibre_length_dorsal`, `fibre_length_ventral` : fibre lengths
- `fibre_density_dorsal`, `fibre_density_ventral` : fibre densities (fibres/area)
- `insulation_reflectance_dorsal`, `insulation_reflectance_ventral` : insulation solar reflectance
- `insulation_depth_compressed` : Compressed insulation depth for ventral insulation
- `ventral_fraction` : Fraction of ventral surface covered by insulation
- `fibre_conductivity` : Thermal conductivity of fibre fibre (energy per time per length per unit temperature)

- `insulation_temperature` : Temperature of insulation
- `ventral_fraction` : Fraction of body surface that is ventral insulation
- `insulation_depth_dorsal`, `insulation_depth_ventral` : Insulation depths (m)

# Returns
`effective_conductivity` (avg, dorsal, ventral)  
`absorption_coefficient` (avg, dorsal, ventral)  
`optical_thickness_factor` (avg, dorsal, ventral)  
`fibre_diameter` (avg, dorsal, ventral)  
`fibre_length` (avg, dorsal, ventral)  
`fibre_density` (avg, dorsal, ventral)  
`insulation_depths` (avg, dorsal, ventral)  
`insulation_reflectance` (avg, dorsal, ventral)
`insulation_test` : Bare-skin test parameter  
`insulation_conductivity_compressed` 
"""
function insulation_properties(; insulation, insulation_temperature, ventral_fraction)

    (; fibre_diameter_dorsal,
    fibre_diameter_ventral, 
    fibre_length_dorsal, 
    fibre_length_ventral,
    insulation_depth_dorsal,
    insulation_depth_ventral,
    fibre_density_dorsal,
    fibre_density_ventral,
    insulation_reflectance_dorsal,
    insulation_reflectance_ventral,
    insulation_depth_compressed,
    fibre_conductivity) = insulation 

    # Physical constants
    air_conductivity = dry_air_properties(insulation_temperature).k_air
    
    # Initialisation
    insulation_conductivity_compressed = 0.0u"W/m/K"
    insulation_test = fibre_density_dorsal * fibre_diameter_dorsal * 
        fibre_length_dorsal * insulation_depth_dorsal  # bare-skin test, TODO note this is dorsal only

    # Weighted averages
    fibre_density = fibre_density_dorsal * (1 - ventral_fraction) + 
        fibre_density_ventral * ventral_fraction
    fibre_diameter = fibre_diameter_dorsal * (1 - ventral_fraction) + 
        fibre_diameter_ventral * ventral_fraction
    fibre_length = fibre_length_dorsal * (1 - ventral_fraction) + 
        fibre_length_ventral * ventral_fraction
    insulation_depth  = insulation_depth_dorsal * (1 - ventral_fraction) + 
        insulation_depth_ventral * ventral_fraction
    insulation_reflectance  = insulation_reflectance_dorsal * (1 - ventral_fraction) + 
        insulation_reflectance_ventral * ventral_fraction

    # Arrays for body regions: 1 = average, 2 = dorsal, 3 = ventral
    effective_conductivities = fill(0.0u"W/m/K", 3)
    absorption_coefficients = fill(0.0u"m^-1", 3)
    optical_thickness_factors = fill(0.0, 3)
    fibre_diameters = fill(0.0u"m", 3)
    fibre_lengths = fill(0.0u"m", 3)
    fibre_densities = fill(0.0u"m^-2", 3)
    insulation_depths = fill(0.0u"m", 3)
    insulation_reflectances = fill(0.0, 3)

    # Average insulation values
    fibre_diameters[1] = fibre_diameter
    fibre_lengths[1] = fibre_length
    fibre_densities[1] = fibre_density
    insulation_depths[1] = insulation_depth
    insulation_reflectances[1] = insulation_reflectance

    # Dorsal values
    fibre_diameters[2] = fibre_diameter_dorsal
    fibre_lengths[2] = fibre_length_dorsal
    fibre_densities[2] = fibre_density_dorsal
    insulation_depths[2] = insulation_depth_dorsal
    insulation_reflectances[2] = insulation_reflectance_dorsal

    # Ventral values (accounting for partial insulation coverage)
    pven_v = min(ventral_fraction * 2.0, 1.0)
    fibre_diameters[3] = fibre_diameter_dorsal * (1 - pven_v) + fibre_diameter_ventral * pven_v
    fibre_lengths[3] = fibre_length_dorsal * (1 - pven_v) + fibre_length_ventral * pven_v
    fibre_densities[3] = fibre_density_dorsal * (1 - pven_v) + fibre_density_ventral * pven_v
    insulation_depths[3] = insulation_depth_dorsal * (1 - pven_v) + insulation_depth_ventral * pven_v
    insulation_reflectances[3] = insulation_reflectance_dorsal * (1 - pven_v) + 
        insulation_reflectance_ventral * pven_v
    # Compute insulation thermal parameters
    for i in 1:3
        if insulation_test <= 0.0u"m"
            effective_conductivities[i] = 0.0u"W/m/K"
            absorption_coefficients[i] = 0.0u"m^-1"
            optical_thickness_factors[i] = 0.0
        else
            (; effective_conductivity, absorption_coefficient, optical_thickness_factor) = 
                insulation_thermal_conductivity(fibre_density = fibre_densities[i], 
                    fibre_length = fibre_lengths[i], insulation_depth = insulation_depths[i],
                    fibre_diameter = fibre_diameters[i], air_conductivity = air_conductivity, 
                    fibre_conductivity = fibre_conductivity)
            effective_conductivities[i] = effective_conductivity
            absorption_coefficients[i] = absorption_coefficient
            optical_thickness_factors[i] = optical_thickness_factor

            # Compressed ventral insulation conductivity
            if i == 3
                (; effective_conductivity) = insulation_thermal_conductivity(fibre_density = 
                    fibre_densities[i], fibre_length = fibre_lengths[i], 
                    insulation_depth = insulation_depth_compressed,
                    fibre_diameter = fibre_diameters[i], 
                    air_conductivity = air_conductivity, fibre_conductivity = fibre_conductivity)
                insulation_conductivity_compressed = effective_conductivity
            end
        end
    end

    return (; effective_conductivities, absorption_coefficients, optical_thickness_factors,
                fibre_diameters, fibre_lengths, fibre_densities, insulation_depths, 
                insulation_reflectances, insulation_test, insulation_conductivity_compressed)
end