"""
    AndrewsPough2 <: MetabolicRateEquation

Empirical model of standard (fasted, inactive period) and resting (fasted, active period)
metabolic rate of squamates, Eq. 2 in Andrews & Pough 1985. Physiol. Zool. 58:214-231.

# Arguments

- `mass_normalisation`: mass normalisation constant.
- `mass_exponent`: mass allometric exponent.
- `thermal_sensitivity`: thermal senstivity constant.
- `metabolic_state`: metabolic state (0 = standard, 1 = resting)
- `O2conversion`: `OxygenJoulesConversion`, `Typical` by default.

The `metabolic_rate` equation reports metabolic rate as oxygen consumption rate given mass
(g) and body temperature (°C). Predictions are capped betwen 1 °C and 50 °C body temperature.
"""
@kwdef struct AndrewsPough2 <: MetabolicRateEquation
    mass_normalisation::Float64 = 0.013
    mass_exponent::Float64 = 0.8
    thermal_sensitivity::Float64 = 0.038
    metabolic_state::Float64 = 0.0
    O2conversion::OxygenJoulesConversion = Typical()
end

"""
    metabolic_rate(eq::MetabolicRateEquation, mass, core_temperature)

Calculate metabolic rate using an empirical allometric equation.

# Arguments
- `eq::MetabolicRateEquation`: The metabolic rate model (e.g., `Kleiber()`, `McKechnieWolf()`, `AndrewsPough2()`)
- `mass`: Body mass
- `core_temperature`: Core temperature (required for `AndrewsPough2`, ignored for others)

# Returns
- `metabolic_heat_flow`: Metabolic heat generation rate (W)
"""
function metabolic_rate(eq::AndrewsPough2, mass, core_temperature)
    (; mass_normalisation, mass_exponent, thermal_sensitivity, metabolic_state, O2conversion) = eq

    mass_g = ustrip(u"g", mass)
    temperature_celsius = clamp(ustrip(u"°C", core_temperature), 1.0, 50.0)
    oxygen_flow = mass_normalisation * mass_g^mass_exponent * 10^(thermal_sensitivity * temperature_celsius) * 10.0 ^ metabolic_state

    metabolic_heat_flow = u"W"(O2_to_Joules(O2conversion, (oxygen_flow)u"ml/hr", 0.8))

    return metabolic_heat_flow
end

"""
    Kleiber <: MetabolicRateEquation

Kleiber's empirical model of basal metabolic rate for mammals.

Kleiber, M. 1947. Body size and metabolic rate. Physiological Reviews 27:511–541.
"""
struct Kleiber <: MetabolicRateEquation end

function metabolic_rate(::Kleiber, mass, core_temperature=nothing)
    mass_kg = ustrip(u"kg", mass)
    m_rate_kcal_day = 70.0 * mass_kg^0.75 # general equation for mammals, p. 535, kcal/day
    metabolic_heat_flow = (m_rate_kcal_day * 4.185 / (24 * 3.6))u"W" # convert to Watts

    return (metabolic_heat_flow)
end

"""
    McKechnieWolf <: MetabolicRateEquation

McKechnie & Wolf's empirical model of basal metabolic rate for birds.

McKechnie, A. E., and B. O. Wolf. 2004. The Allometry of Avian Basal Metabolic Rate:
    Good Predictions Need Good Data. Physiological and Biochemical Zoology:
    Ecological and Evolutionary Approaches 77:502–521.
"""
struct McKechnieWolf <: MetabolicRateEquation end

function metabolic_rate(::McKechnieWolf, mass, core_temperature=nothing)
    mass_g = ustrip(u"g", mass)
    metabolic_heat_flow = (10.0^(-1.461 + 0.669 * log10(mass_g)))u"W"

    return (metabolic_heat_flow)
end

"""
    PlantDarkRespiration <: MetabolicRateEquation

Leaf/plant dark (mitochondrial) respiration using mass scaling and Arrhenius
temperature dependence (Reich et al. 2006; CLM/LPJ parameterisation).

    R = mass_normalisation × mass^mass_exponent
        × exp(-activation_energy / (k × T))
        / exp(-activation_energy / (k × reference_temperature))

# Parameters
- `mass_normalisation` — respiration rate per unit mass at reference temperature (W/kg)
- `mass_exponent` — allometric exponent; default 1.0 (near-isometric for plants; Reich 2006)
- `activation_energy` — Arrhenius activation energy (eV); default 0.65 eV
- `reference_temperature` — normalisation temperature (K); default 298.15 K (25 °C)

# References
- Reich, P. B., et al. (2006). Universal scaling of respiratory metabolism, size and nitrogen in plants.
  *Nature*, 439, 457–461.
"""
@kwdef struct PlantDarkRespiration <: MetabolicRateEquation
    mass_normalisation::Float64    = 4.6e-4    # W/kg at 25 °C; broad C3 leaf estimate
    mass_exponent::Float64         = 1.0
    activation_energy::Float64     = 0.65      # eV
    reference_temperature::Float64 = 298.15    # K
end

function metabolic_rate(eq::PlantDarkRespiration, mass, temperature)
    (; mass_normalisation, mass_exponent, activation_energy, reference_temperature) = eq
    k_eV = 8.617333e-5   # Boltzmann constant, eV/K
    T    = ustrip(u"K", temperature)
    M    = ustrip(u"kg", mass)
    rate = mass_normalisation * M^mass_exponent *
           exp(-activation_energy / (k_eV * T)) /
           exp(-activation_energy / (k_eV * reference_temperature))
    return rate * u"W"
end

"""
    metabolic_rate(::Nothing, mass, body_temperature)

Returns zero metabolic heat production. Use `model = nothing` in
`MetabolismParameters` to disable metabolic heat generation entirely.
"""
metabolic_rate(::Nothing, mass, body_temperature) = 0.0u"W"

# conversion from O2 consumption to Joules

"""
    Typical <: OxygenJoulesConversion

Standard conversion between O₂ consumption and energy use.

# Arguments (Joules_to_O2) / Outputs (O2_to_Joules)
- `metabolic_heat_flow` — metabolic rate (W).

# Outputs (Joules_to_O2) / Arguments (O2_to_Joules)
- `oxygen_flow_standard` — oxygen consumption rate at standard temperature and pressure (L/s)

TODO - add a ref?
"""
struct Typical <: OxygenJoulesConversion end

"""
    O2_to_Joules(conversion::OxygenJoulesConversion, oxygen_flow_standard, rq)
    O2_to_Joules(oxygen_flow_standard)

Convert oxygen consumption rate to metabolic heat rate.

# Arguments
- `conversion::OxygenJoulesConversion`: Conversion model (e.g., `Typical()`, `Kleiber1961()`)
- `oxygen_flow_standard`: Oxygen consumption rate at standard temperature and pressure
- `rq`: Respiratory quotient (CO₂ produced / O₂ consumed)

# Returns
- `metabolic_heat_flow`: Metabolic heat rate (W)
"""
function O2_to_Joules(::Typical, oxygen_flow_standard, rq)
    O2_energy_equivalent = 20.1u"J/ml"
    metabolic_heat_flow = oxygen_flow_standard * O2_energy_equivalent
    return (metabolic_heat_flow)
end

"""
    Joules_to_O2(conversion::OxygenJoulesConversion, metabolic_heat_flow, rq)
    Joules_to_O2(metabolic_heat_flow)

Convert metabolic heat rate to oxygen consumption rate.

# Arguments
- `conversion::OxygenJoulesConversion`: Conversion model (e.g., `Typical()`, `Kleiber1961()`)
- `metabolic_heat_flow`: Metabolic heat rate (W)
- `rq`: Respiratory quotient (CO₂ produced / O₂ consumed)

# Returns
- `oxygen_flow_standard`: Oxygen consumption rate at standard temperature and pressure
"""
function Joules_to_O2(::Typical, metabolic_heat_flow, rq)
    O2_energy_equivalent = 20.1u"J/ml"
    oxygen_flow_standard = metabolic_heat_flow / O2_energy_equivalent
    return (oxygen_flow_standard)
end

"""
    Kleiber1961 <: OxygenJoulesConversion

Kleiber's conversion between O₂ consumption and energy use, from
the bottom of Table 7.3.

# Arguments (Joules_to_O2) / Outputs (O2_to_Joules)
- `metabolic_heat_flow` — metabolic rate (W).

# Outputs (Joules_to_O2) / Arguments (O2_to_Joules)
- `oxygen_flow_standard` — oxygen consumption rate at standard temperature and pressure (L/s)

Kleiber, M. 1961. The Fire of Life. An Introduction to Animal Energetics.
"""
struct Kleiber1961 <: OxygenJoulesConversion end

function Joules_to_O2(::Kleiber1961, metabolic_heat_flow, rq)
    carbohydrate_energy_equivalent = u"J/ml"(5.0u"kcal"/1u"L")
    fat_energy_equivalent = u"J/ml"(4.7u"kcal"/1u"L")
    protein_energy_equivalent = u"J/ml"(4.5u"kcal"/1u"L")
    if rq ≥ 1.0
        oxygen_flow_standard = metabolic_heat_flow / carbohydrate_energy_equivalent
    elseif rq ≤ 0.7
        oxygen_flow_standard = metabolic_heat_flow / fat_energy_equivalent
    else
        oxygen_flow_standard = metabolic_heat_flow / protein_energy_equivalent
    end
    return (oxygen_flow_standard)
end

function O2_to_Joules(::Kleiber1961, oxygen_flow_standard, rq)
    carbohydrate_energy_equivalent = u"J/ml"(5.0u"kcal"/1u"L")
    fat_energy_equivalent = u"J/ml"(4.7u"kcal"/1u"L")
    protein_energy_equivalent = u"J/ml"(4.5u"kcal"/1u"L")
    if rq ≥ 1.0
        metabolic_heat_flow = oxygen_flow_standard * carbohydrate_energy_equivalent
    elseif rq ≤ 0.7
        metabolic_heat_flow = oxygen_flow_standard * fat_energy_equivalent
    else
        metabolic_heat_flow = oxygen_flow_standard * protein_energy_equivalent
    end
    return (metabolic_heat_flow)
end


# Fallbacks: why are these needed ?
O2_to_Joules(::Missing) = missing
Joules_to_O2(::Missing) = missing

Joules_to_O2(metabolic_heat_flow) = Joules_to_O2(Typical(), metabolic_heat_flow, one(metabolic_heat_flow))
O2_to_Joules(oxygen_flow_standard) = O2_to_Joules(Typical(), oxygen_flow_standard, one(oxygen_flow_standard))

# TODO
# """
#     ElliottDavison <: OxygenJoulesConversion

# Elliott and Davison conversion between O₂ consumption and energy use.

# # Arguments (Joules_to_O2) / Outputs (O2_to_Joules)
# - `metabolic_flux` — metabolic rate (W).

# # Outputs (Joules_to_O2) / Arguments (O2_to_Joules)
# - V_O2_STP, oxygen consumption rate, standard temperature and pressure, L/s

# Elliott, J. M., and W. Davison. 1975. Oecologia 19:195–201.
# """
# struct ElliottDavison <: OxygenJoulesConversion end

# function Joules_to_O2(::ElliottDavison, metabolic_flux, rq)
#     if rq ≥ 1.0
#         # carbohydrate metabolism
#         # carbohydrates ≈ 4193 cal/g
#         # L/s = (J/s) * (cal/J) * (kcal/cal) * (L O2 / kcal)
#         V_O2_STP = metabolic_flux * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/5.057u"kcal")
#     end
#     if rq ≤ 0.7
#         # fat metabolism; fats ≈ 9400 cal/g
#         V_O2_STP = metabolic_flux * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.7u"kcal")
#     else
#         # protein metabolism (RQ ≈ 0.8); ≈ 4300 cal/g
#         V_O2_STP = metabolic_flux * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.5u"kcal")
#     end
#     return (V_O2_STP)
# end

# function O2_to_Joules(::ElliottDavison, V_O2_STP, rq)
#     if rq ≥ 1.0
#         # carbohydrate metabolism
#         # carbohydrates ≈ 4193 cal/g
#         # L/s = (J/s) * (cal/J) * (kcal/cal) * (L O2 / kcal)
#         metabolic_flux  = V_O2_STP / ((1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/5.057u"kcal"))
#     end
#     if rq ≤ 0.7
#         # fat metabolism; fats ≈ 9400 cal/g
#         metabolic_flux  = V_O2_STP / ((1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.7u"kcal"))
#     else
#         # protein metabolism (RQ ≈ 0.8); ≈ 4300 cal/g
#         metabolic_flux  = V_O2_STP / ((1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.5u"kcal"))
#     end
#     return (metabolic_flux)
# end
