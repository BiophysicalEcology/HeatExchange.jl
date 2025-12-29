"""
    ellipsoid(; posture, mass, density, core_temperature,
               insulation_depth, k_insulation, oxygen_extraction_efficiency, stress_factor,
               air_temperature, wind_speed, relative_humidity, q10,
               minimum_metabolic_rate=missing, P_atmos, metabolic_factor=1)

Ellipsoid endotherm model.

Implements the model from Porter & Kearney (2009) *PNAS* 106:19666–19672,  
with rough water loss estimates. Valid for conditions without solar radiation and  
where air, ground, and sky temperatures are equal.

# Arguments
- `posture`: Ratio of long to short axis of a prolate ellipsoid.
- `mass`: Body mass (mass).
- `density`: Body density (mass/volume).
- `core_temperature`: Core temperature.
- `insulation_depth`: Insulation depth (length).
- `k_insulation`: Insulation conductivity (energy/time/length/temperature).
- `oxygen_extraction_efficiency`: Oxygen extraction efficiency (fraction).
- `stress_factor`: Fraction of basal metabolism at which evaporative water loss begins.
- `air_temperature`: Air temperature.
- `wind_speed`: Wind speed (length/time).
- `relative_humidity`: Relative humidity (%).
- `P_atmos`: Atmospheric_pressue (pressure).
- `q10`: Temperature dependence factor for metabolism.
- `minimum_metabolic_rate`: Optional user-specified minimum metabolic rate (energy/time).
- `f_O2`: Fractional oxygen in the atmosphere (-).
- `metabolic_multiplier`: Scaling factor for the mouse–elephant basal rate.

# Returns
A `NamedTuple` with:
- `skin_temperature`, `upper_critical_air_temperature`, `lower_critical_air_temperature`,
- `Q_gen`, `Q_gen_final`, `Q_respiration`, `Q_evap`,
- `O2_consumption_rate`, `respiratory_water_loss_rate`, `total_water_loss_rate`,
- `basal_metabolic_rate_fraction`, `fractional_mass_loss`.

# References
Porter, W. P., & Kearney, M. R. (2009). *Size, shape, and the thermal niche of endotherms*.  
PNAS, 106(46), 19666–19672.
"""
function ellipsoid_endotherm(;
    mrate_equation::MetabolicRateEquation = Kleiber(),
    posture,
    mass,
    density,
    core_temperature,
    insulation_depth,
    insulation_conductivity,
    emissivity,
    oxygen_extraction_efficiency,
    stress_factor,
    air_temperature,
    wind_speed,
    relative_humidity,
    P_atmos,
    q10,
    minimum_metabolic_rate=missing,
    metabolic_multiplier=1,
    lethal_desiccation,
    f_O2,
)

    # refinition of variables and parameters to match Porter & Kearney 2009 formulae
    T_f = air_temperature |> u"K" # fluid temperature
    T_c = core_temperature |> u"K" # core body temperature
    k_insulation = insulation_conductivity
    ϵ = emissivity
    v = wind_speed
    Q_gen_min = minimum_metabolic_rate

    # avoid divide-by-zero
    posture = posture == 1 ? 1.01 : posture

    # estimate basal metabolism if not provided
    if isnothing(minimum_metabolic_rate) || minimum_metabolic_rate === missing
        allometric_estimate = metabolic_rate(mrate_equation, mass) * metabolic_multiplier
        Q_gen_min = allometric_estimate * q10^((ustrip(u"°C", core_temperature) - 37) / 10)
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

    k_b = (0.5 + (6.14 * ustrip(u"m", b)) + 0.439)u"W/m/K" # thermal conductivity of body
    numerator = a^2 * b^2 * c^2
    denominator = a^2 * b^2 + a^2 * c^2 + b^2 * c^2
    S2 = numerator / denominator # Eq. 5
    R_b = S2 / (2 * k_b * V) # resistance of body

    a_o = b * posture + insulation_depth
    b_o = b + insulation_depth
    c_o = c + insulation_depth

    e_o = sqrt(a_o^2 - c_o^2) / a_o
    # outer area of ellipsoid # TODO use general function
    A_o = 2 * π * b_o^2 + 2 * π * ((a_o * b_o) / e_o) * asin(e_o)
    R_ins = (b_o - b) / (k_insulation * A_o) # insulation radius

    # air properties
    (; μ, k_air, ρ_air) = dry_air_properties(T_f, P_atmos)

    V = (4 / 3) * π * a * b * c # recompute volume
    L_c = V^(1 / 3) # characteristic dimension
    e = sqrt(a^2 - c^2) / a # eccentricity
    A = 2 * π * b^2 + 2 * π * ((a * b) / e) * asin(e) # area of skin

    Re = ρ_air * v * L_c / μ # Reynold's number
    Pr = (μ * c_p_air) / k_air # Prandtl number

    q′′′_numerator = 2 * A * k_b * k_air * (2 + a_coef * (Re^b_coef) * Pr^(1 / 3)) * (T_c - T_f)
    q′′′_denominator = 2 * k_b * L_c * V + A * S2 * k_air * (2 + a_coef * Re^b_coef * Pr^(1 / 3))
    q′′′ = q′′′_numerator / q′′′_denominator
    T_s = T_c - (q′′′ * S2) / (2 * k_b) # skin temperature, Eq. 6

    Gr = abs((ρ_air^2) * (1 / T_f) * Unitful.gn * (L_c^3) * (T_s - T_f) / μ^2) # Grashof number
    Nu_free = 2 + 0.6 * (Gr^0.25) * (Pr^(1 / 3))
    Nu_forced = 0.37 * Re^0.6
    Nu_total = (Nu_free^3 + Nu_forced^3)^(1 / 3) #  Eq. 24
    h_cv = Nu_total * k_air / L_c # Eq. 25
    R_cv = 1 / (h_cv * A_o)
    R_rad = u"K/W"(1 / (4 * A_o * ϵ * Unitful.σ * T_f^3)) # Eq. 39
    R_total = R_b + R_ins + (R_cv * R_rad) / (R_cv + R_rad)

    upper_critical_air_temperature = T_c - (Q_gen_min * stress_factor * R_total)
    lower_critical_air_temperature = T_c - Q_gen_min * R_total

    Q_gen_required = (T_c - T_f) / R_total
    Q_gen_final = max(Q_gen_required, Q_gen_min)

    O2_consumption_rate = u"ml/hr"(Joules_to_O2(Q_gen_final))
    ρ_vap_f = wet_air_properties(T_f, relative_humidity, P_atmos).ρ_vap # inhaled air
    ρ_vap_c = wet_air_properties(T_c, 1, P_atmos).ρ_vap # exhaled air (saturated)

    respiratory_water_loss_rate = u"g/hr"((O2_consumption_rate / f_O2 / oxygen_extraction_efficiency) *
                                  (ρ_vap_c - ρ_vap_f))
    latent_heat = enthalpy_of_vaporisation(T_f)
    Q_respiration = u"W"(respiratory_water_loss_rate * latent_heat)

    basal_metabolic_rate_fraction = Q_gen_final / Q_gen_min

    Q_evap = -Q_gen_required + Q_gen_min
    total_water_loss_rate = u"g/hr"(max((u"J/hr"(Q_evap) / latent_heat), 0.0u"kg/hr"))
    fractional_mass_loss = u"kg/hr"(total_water_loss_rate) / mass

    return (;
        skin_temperature=u"°C"(T_s),
        upper_critical_air_temperature=u"°C"(upper_critical_air_temperature),
        lower_critical_air_temperature=u"°C"(lower_critical_air_temperature),
        Q_gen_required,
        Q_gen_final,
        Q_respiration,
        Q_evap,
        O2_consumption_rate=u"mL/hr"(O2_consumption_rate),
        respiratory_water_loss_rate=u"g/hr"(respiratory_water_loss_rate),
        total_water_loss_rate=u"g/hr"(total_water_loss_rate),
        basal_metabolic_rate_fraction,
        fractional_mass_loss,
    )
end
