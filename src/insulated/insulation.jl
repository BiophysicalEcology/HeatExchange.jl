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
insulation_thermal_conductivity(fibre::FibreProperties, air_conductivity; smoothing::SmoothingStrategy=HardBound()) =
    insulation_thermal_conductivity(fibre, air_conductivity, fibre.depth; smoothing)

function insulation_thermal_conductivity(fibre::FibreProperties, air_conductivity, depth;
                                         smoothing::SmoothingStrategy=HardBound())
    # Canonicalise inputs to SI before any arithmetic. Otherwise
    # `density * (fibre.length / insulation_depth)` carries the literal
    # unit ratio of whichever inputs were used (e.g. mm/m gives a non-
    # canonical FreeUnits) and every downstream type becomes call-site-
    # dependent. Caller-side unit choices then make the whole function
    # type-unstable from a downstream view.
    insulation_depth = uconvert(u"m", depth)
    fibre_length     = uconvert(u"m", fibre.length)
    fibre_diameter   = uconvert(u"m", fibre.diameter)
    density          = uconvert(u"m^-2", fibre.density)

    w = 1.0  # weighting factor placeholder

    # Effective density (Conley & Porter 1986, p.254)
    effective_density = density * (fibre_length / insulation_depth)

    # Distance between fibre centers (assuming uniform spacing)
    l_unit = 1.0 / sqrt(effective_density)
    fibre_spacing = l_unit - fibre_diameter

    # Resolve crowded-fibre case (no space between fibres) by re-deriving
    # density from a hard-packed lattice. The reassignment must preserve
    # the element type (e.g. `Dual` under AD) — naively setting
    # `l_unit = fibre_diameter` would collapse the branch back to
    # `Float64` even when the open case carried a `Dual`, producing a
    # `Union` join in inference and ~kB / call of boxed allocations.
    # `+ zero(l_unit)` is a no-op value-wise but promotes to the wider
    # type.
    if fibre_spacing < 0.0u"m"
        l_unit = fibre_diameter + zero(l_unit)
        effective_density = uconvert(u"m^-2", 1.0 / l_unit^2)
        density = uconvert(u"m^-2", (effective_density * insulation_depth) / fibre_length)
        fibre_spacing = l_unit - fibre_diameter
    end

    a_air = (w / 2.0) * (l_unit - fibre_diameter)
    # Cast r_air, r_fibre to a fixed unit (K·m/W) in both branches —
    # without the cast the two arms have different FreeUnits orderings,
    # which makes ky / effective_conductivity / the whole return Union-typed.
    r_air = if fibre_spacing > 0.0u"m"
        u"K*m*W^-1"(l_unit / (air_conductivity * a_air))
    else
        u"K*m*W^-1"(2.0 / (sqrt(effective_density) * air_conductivity * (l_unit - fibre_diameter) * w))
    end
    r_fibre = u"K*m*W^-1"(
        ((fibre_diameter * air_conductivity) + (l_unit - fibre_diameter) * fibre.conductivity) /
        (w * fibre_diameter * fibre.conductivity * air_conductivity)
    )
    # `a_fibre` is the dimensionless fibre-cross-section fraction.
    # NoUnits forces the result to plain Float64; without it the two
    # branches return Quantity{NoDims, ...} with different FreeUnits
    # orderings (cm-derived vs m-derived), making the local Union-typed.
    a_fibre = if fibre_spacing > 0.0u"m"
        NoUnits(density * (u"cm"(fibre_diameter) / 2.0)^2 * π)
    else
        NoUnits(effective_density * (fibre_diameter / 2.0)^2 * π)
    end
    kx = a_fibre * fibre.conductivity + (1.0 - a_fibre) * air_conductivity
    ky = (2.0 / r_air) + (1.0 / r_fibre)

    # Effective thermal conductivity (eq. 3-28, Kowalski 1983).
    # Cast to canonical W/m/K. The clamp below would otherwise mix three
    # different FreeUnits orderings (computed `(ky+kx)/2`, fibre.conductivity,
    # air_conductivity) and produce a Union return type.
    effective_conductivity = u"W/m/K"((ky + kx) / 2.0)
    effective_conductivity = safe_clamp(smoothing, effective_conductivity,
        u"W/m/K"(air_conductivity), u"W/m/K"(fibre.conductivity); scale=u"W/m/K"(air_conductivity))

    # Absorption and optical parameters
    absorption_coefficient = (0.67 / π) * u"m^-2"(effective_density) * u"m"(fibre_diameter)
    optical_thickness_factor = absorption_coefficient * u"m"(insulation_depth)

    return (; effective_conductivity, absorption_coefficient, optical_thickness_factor)
end

"""
    insulation_properties(insulation, insulation_temperature, ventral_fraction) -> InsulationProperties

Compute parameters for heat conduction and longwave radiation through insulation (fur or plumage).

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
function insulation_properties(insulation::InsulationParameters, insulation_temperature, ventral_fraction;
                               smoothing::SmoothingStrategy=HardBound())
    (; dorsal, ventral, depth_compressed) = insulation

    # Physical constants
    air_conductivity = dry_air_properties(insulation_temperature).thermal_conductivity

    # Bare-skin test (dorsal only)
    insulation_test = (dorsal.density * dorsal.diameter) * (dorsal.length * dorsal.depth)

    # Weighted average fibre properties
    avg_fibres = interpolate(dorsal, ventral, ventral_fraction)

    # Adjusted ventral values (accounting for partial insulation coverage).
    # `ventral_fraction` is dimensionless in [0, 1]; the kink is at `2v = 1`.
    # `scale = 0.01` (vs the default `oneunit(...) = 1.0`) tightens the
    # SmoothBound transition to a fractional width of 0.01·ε around the kink —
    # well inside the [0, 1] range without overlapping the physical extremes.
    ventral_adj = interpolate(dorsal, ventral, safe_min(smoothing, ventral_fraction * 2.0, 1.0; scale=0.01))

    fibres = BodyRegionValues(avg_fibres, dorsal, ventral_adj)

    # Bare-skin (zero density / depth / length) makes `insulation_thermal_conductivity`
    # divide by zero, so we can't use a `mask * value` blend — IEEE keeps NaN through
    # `0 * NaN`. We need a Union-free return type under AD (Dual through
    # `air_conductivity` and `fibres`), so both branches must produce the same
    # concrete component types. The bare-skin branch can't use literal `0.0u"..."`
    # because the canonical-W/m/K branch's Unitful FreeUnits ordering differs from
    # a literal's. We anchor the types by calling `insulation_thermal_conductivity`
    # once on `fibres.average` (always defined), then `zero(...)` of its component
    # fields gives the bare-skin values. `zero(NaN) == 0.0`, so the NaN-producing
    # call on truly-bare-skin inputs is harmless: we consume the type, not the value.
    if insulation_test > zero(insulation_test)
        avg  = insulation_thermal_conductivity(fibres.average, air_conductivity; smoothing)
        dors = insulation_thermal_conductivity(fibres.dorsal,  air_conductivity; smoothing)
        vent = insulation_thermal_conductivity(fibres.ventral, air_conductivity; smoothing)
        vent_compressed = insulation_thermal_conductivity(fibres.ventral, air_conductivity, depth_compressed; smoothing)
        conductivities = BodyRegionValues(
            avg.effective_conductivity, dors.effective_conductivity, vent.effective_conductivity,
        )
        absorption_coefficients = BodyRegionValues(
            avg.absorption_coefficient, dors.absorption_coefficient, vent.absorption_coefficient,
        )
        optical_thickness = BodyRegionValues(
            avg.optical_thickness_factor, dors.optical_thickness_factor, vent.optical_thickness_factor,
        )
        conductivity_compressed = vent_compressed.effective_conductivity
    else
        # Single anchor call for type-only consumption.
        anchor = insulation_thermal_conductivity(fibres.average, air_conductivity; smoothing)
        zk = zero(anchor.effective_conductivity)
        za = zero(anchor.absorption_coefficient)
        zo = zero(anchor.optical_thickness_factor)
        conductivities          = BodyRegionValues(zk, zk, zk)
        absorption_coefficients = BodyRegionValues(za, za, za)
        optical_thickness       = BodyRegionValues(zo, zo, zo)
        conductivity_compressed = zk
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
