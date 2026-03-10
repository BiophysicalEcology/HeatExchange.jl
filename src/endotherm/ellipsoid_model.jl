"""
    ellipsoid_endotherm(; posture, mass, density, core_temperature,
               insulation_depth, insulation_conductivity, oxygen_extraction_efficiency, stress_factor,
               air_temperature, wind_speed, relative_humidity, q10,
               minimum_metabolic_rate=missing, atmospheric_pressure, metabolic_factor=1)

Ellipsoid endotherm model.

Implements the model from Porter & Kearney (2009) *PNAS* 106:19666–19672,
with rough water loss estimates. Valid for conditions without solar radiation and
where air, ground, and sky temperatures are equal.

# Arguments
- `air_temperature`: Air temperature.
- `wind_speed`: Wind speed (length/time).
- `relative_humidity`: Relative humidity (fractional).
- `atmospheric_pressure`: Atmospheric_pressue (pressure).
- `mass`: Body mass (mass).
- `density`: Body density (mass/volume).
- `posture`: Ratio of long to short axis of a prolate ellipsoid.
- `insulation_depth`: Insulation depth (length).
- `insulation_conductivity`: Insulation conductivity (energy/time/length/temperature).
- `emissivity`: Insulation emissivity (-)
- `core_temperature`: Core temperature.
- `minimum_metabolic_rate`: Optional user-specified minimum metabolic rate (energy/time).
- `metabolic_rate_equation`: Allometric function to compute metabolic rate (default: Kleiber())
- `metabolic_multiplier`: Scaling factor for the metabolic_rate_equation, 
    e.g. to adjust to a field metabolic rate.
- `q10`: Temperature dependence factor for metabolism.
- `oxygen_fraction`: Fractional oxygen in the atmosphere (-).
- `oxygen_extraction_efficiency`: Oxygen extraction efficiency (fraction).
- `stress_factor`: Fraction of basal metabolism at which evaporative water loss begins.

# Returns
A `NamedTuple` with:
- `skin_temperature`, `upper_critical_air_temperature`, `lower_critical_air_temperature`,
- `generated_flux`, `final_generated_flux`, `respiration_heat_flow`, `evaporation_heat_flow`,
- `oxygen_consumption_rate`, `respiratory_water_loss_rate`, `total_water_loss_rate`,
- `basal_metabolic_rate_fraction`, `fractional_mass_loss`.

# References
Porter, W. P., & Kearney, M. R. (2009). *Size, shape, and the thermal niche of endotherms*.  
PNAS, 106(46), 19666–19672.
"""
ellipsoid_endotherm(::Missing, ::Missing, ::Missing, ::Missing; kwargs...) = missing
function ellipsoid_endotherm(
    air_temperature::Quantity,
    wind_speed::Quantity,
    relative_humidity::Real,
    atmospheric_pressure::Quantity;
    mass::Quantity,
    density::Quantity,
    posture::Real,
    insulation_depth::Quantity,
    insulation_conductivity::Quantity,
    emissivity::Real,
    core_temperature::Quantity,
    minimum_metabolic_rate=missing,
    metabolic_rate_equation::MetabolicRateEquation=Kleiber(),
    metabolic_multiplier::Real=1.0,
    q10::Real,
    oxygen_fraction::Real,
    oxygen_extraction_efficiency::Real,
    stress_factor::Real,
)
    ellipsoid_endotherm(
        air_temperature,
        wind_speed,
        relative_humidity,
        atmospheric_pressure,
        mass,
        density,
        posture,
        insulation_depth,
        insulation_conductivity,
        emissivity,
        core_temperature,
        minimum_metabolic_rate,
        metabolic_rate_equation,
        metabolic_multiplier,
        q10,
        oxygen_fraction,
        oxygen_extraction_efficiency,
        stress_factor,
    )
end
function ellipsoid_endotherm(
    air_temperature::Quantity,
    wind_speed::Quantity,
    relative_humidity::Real,
    atmospheric_pressure::Quantity,
    mass::Quantity,
    density::Quantity,
    posture::Real,
    insulation_depth::Quantity,
    insulation_conductivity::Quantity,
    emissivity::Real,
    core_temperature::Quantity,
    minimum_metabolic_rate::Union{Missing,Quantity},
    metabolic_rate_equation::MetabolicRateEquation,
    metabolic_multiplier::Real,
    q10::Real,
    oxygen_fraction::Real,
    oxygen_extraction_efficiency::Real,
    stress_factor::Real,
)

    # local aliases to match Porter & Kearney 2009 formulae
    ϵ = emissivity
    v = wind_speed
    minimum_generated_flux = minimum_metabolic_rate

    # avoid divide-by-zero
    posture = posture == 1 ? 1.01 : posture

    # estimate basal metabolism if not provided
    if isnothing(minimum_metabolic_rate) || minimum_metabolic_rate === missing
        allometric_estimate =
            metabolic_rate(metabolic_rate_equation, mass, core_temperature) * metabolic_multiplier
        minimum_generated_flux = allometric_estimate * q10^((ustrip(u"°C", core_temperature) - 37) / 10)
    end

    # constants
    a_coef = 0.6
    b_coef = 0.5
    c_p_air = 1005.8u"J/kg/K"

    # body geometry
    V = mass / density # volume
    b = ((3 * V) / (4 * π * posture))^(1 / 3)
    c = b
    a = b * posture

    body_conductivity = (0.5 + (6.14 * ustrip(u"m", b)) + 0.439)u"W/m/K" # thermal conductivity of body
    numerator = a^2 * b^2 * c^2
    denominator = a^2 * b^2 + a^2 * c^2 + b^2 * c^2
    S2 = numerator / denominator # Eq. 5
    R_b = S2 / (2 * body_conductivity * V) # resistance of body

    a_o = b * posture + insulation_depth
    b_o = b + insulation_depth
    c_o = c + insulation_depth

    e_o = sqrt(a_o^2 - c_o^2) / a_o
    # outer area of ellipsoid # TODO use general function
    outer_area = 2 * π * b_o^2 + 2 * π * ((a_o * b_o) / e_o) * asin(e_o)
    R_ins = (b_o - b) / (insulation_conductivity * outer_area) # insulation resistance

    # air properties
    (; dynamic_viscosity, thermal_conductivity, density) = dry_air_properties(air_temperature, atmospheric_pressure)
    μ, air_conductivity, ρ_air = dynamic_viscosity, thermal_conductivity, density

    V = (4 / 3) * π * a * b * c # recompute volume
    L_c = V^(1 / 3) # characteristic dimension
    e = sqrt(a^2 - c^2) / a # eccentricity
    A = 2 * π * b^2 + 2 * π * ((a * b) / e) * asin(e) # area of skin

    reynolds_number = ρ_air * v * L_c / μ # Reynolds number
    prandtl_number = (μ * c_p_air) / air_conductivity # Prandtl number

    q′′′_numerator =
        2 * A * body_conductivity * air_conductivity * (2 + a_coef * (reynolds_number^b_coef) * prandtl_number^(1 / 3)) * (core_temperature - air_temperature)
    q′′′_denominator =
        2 * body_conductivity * L_c * V + A * S2 * air_conductivity * (2 + a_coef * reynolds_number^b_coef * prandtl_number^(1 / 3))
    q′′′ = q′′′_numerator / q′′′_denominator
    skin_temperature = core_temperature - (q′′′ * S2) / (2 * body_conductivity) # skin temperature, Eq. 6

    grashof_number = abs((ρ_air^2) * (1 / air_temperature) * Unitful.gn * (L_c^3) * (skin_temperature - air_temperature) / μ^2) # Grashof number
    free_nusselt_number = 2 + 0.6 * (grashof_number^0.25) * (prandtl_number^(1 / 3))
    forced_nusselt_number = 0.37 * reynolds_number^0.6
    combined_nusselt_number = (free_nusselt_number^3 + forced_nusselt_number^3)^(1 / 3) #  Eq. 24
    h_cv = combined_nusselt_number * air_conductivity / L_c # Eq. 25
    R_cv = 1 / (h_cv * outer_area)
    R_rad = u"K/W"(1 / (4 * outer_area * ϵ * Unitful.σ * air_temperature^3)) # Eq. 39
    R_total = R_b + R_ins + (R_cv * R_rad) / (R_cv + R_rad)

    upper_critical_air_temperature = core_temperature - (minimum_generated_flux * stress_factor * R_total)
    lower_critical_air_temperature = core_temperature - minimum_generated_flux * R_total

    required_generated_flux = (core_temperature - air_temperature) / R_total
    final_generated_flux = max(required_generated_flux, minimum_generated_flux)

    oxygen_consumption_rate = u"ml/hr"(Joules_to_O2(final_generated_flux))
    ρ_vap_f = wet_air_properties(air_temperature, relative_humidity, atmospheric_pressure).vapour_density # inhaled air
    ρ_vap_c = wet_air_properties(core_temperature, 1, atmospheric_pressure).vapour_density # exhaled air (saturated)

    respiratory_water_loss_rate = u"g/hr"(
        (oxygen_consumption_rate / oxygen_fraction / oxygen_extraction_efficiency) * (ρ_vap_c - ρ_vap_f)
    )
    latent_heat = enthalpy_of_vaporisation(air_temperature)
    respiration_heat_flow = u"W"(respiratory_water_loss_rate * latent_heat)

    basal_metabolic_rate_fraction = final_generated_flux / minimum_generated_flux

    evaporation_heat_flow = -required_generated_flux + minimum_generated_flux
    total_water_loss_rate = u"g/hr"(max((u"J/hr"(evaporation_heat_flow) / latent_heat), 0.0u"kg/hr"))
    fractional_mass_loss = u"kg/hr"(total_water_loss_rate) / mass

    return (;
        skin_temperature=u"°C"(skin_temperature),
        upper_critical_air_temperature=u"°C"(upper_critical_air_temperature),
        lower_critical_air_temperature=u"°C"(lower_critical_air_temperature),
        required_generated_flux,
        final_generated_flux,
        respiration_heat_flow,
        evaporation_heat_flow,
        oxygen_consumption_rate=u"mL/hr"(oxygen_consumption_rate),
        respiratory_water_loss_rate=u"g/hr"(respiratory_water_loss_rate),
        total_water_loss_rate=u"g/hr"(total_water_loss_rate),
        basal_metabolic_rate_fraction,
        fractional_mass_loss,
    )
end
