"""
    AndrewsPough2 <: MetabolicRateEquation

Empirical model of standard (fasted, inactive period) and resting (fasted, active period)
metabolic rate of squamates, Eq. 2 in Andrews & Pough 1985. Physiol. Zool. 58:214-231.

# Arguments

- `M1`: mass normalisation constant.
- `M2`: mass allometric exponent. 
- `M3`: thermal senstivity constant. 
- `M4`: metabolic state (0 = standard, 1 = resting) 
- `O2conversion`: `OxygenJoulesConversion`, `Typical` by default.

The `metabolic_rate` equation reports metabolic rate as oxygen consumption rate given mass
(g) and body temperature (°C). Predictions are capped betwen 1 °C and 50 °C body temperature.
"""
@kwdef struct AndrewsPough2 <: MetabolicRateEquation
    M1::Float64 = 0.013
    M2::Float64 = 0.8
    M3::Float64 = 0.038
    M4::Float64 = 0.0
    O2conversion::OxygenJoulesConversion = Typical()
end

"""
    metabolic_rate(eq::MetabolicRateEquation, mass, body_temperature)

Calculate metabolic rate using an empirical allometric equation.

# Arguments
- `eq::MetabolicRateEquation`: The metabolic rate model (e.g., `Kleiber()`, `McKechnieWolf()`, `AndrewsPough2()`)
- `mass`: Body mass
- `body_temperature`: Body temperature (required for `AndrewsPough2`, ignored for others)

# Returns
- `metabolic_flux`: Metabolic heat generation rate (W)
"""
function metabolic_rate(eq::AndrewsPough2, mass, body_temperature)
    (; M1, M2, M3, M4, O2conversion) = eq

    mass_g = ustrip(u"g", mass)
    T_c = clamp(ustrip(u"°C", body_temperature), 1.0, 50.0)
    V_O2 = M1 * mass_g^M2 * 10^(M3 * T_c) * 10.0 ^ M4

    metabolic_flux = u"W"(O2_to_Joules(O2conversion, (V_O2)u"ml/hr", 0.8))

    return metabolic_flux
end

"""
    Kleiber <: MetabolicRateEquation

Kleiber's empirical model of basal metabolic rate for mammals.

Kleiber, M. 1947. Body size and metabolic rate. Physiological Reviews 27:511–541.
"""
struct Kleiber <: MetabolicRateEquation end

function metabolic_rate(::Kleiber, mass, body_temperature=nothing)
    mass_kg = ustrip(u"kg", mass)
    m_rate_kcal_day = 70.0 * mass_kg^0.75 # general equation for mammals, p. 535, kcal/day
    metabolic_flux = (m_rate_kcal_day * 4.185 / (24 * 3.6))u"W" # convert to Watts

    return (metabolic_flux)
end

"""
    McKechnieWolf <: MetabolicRateEquation

McKechnie & Wolf's empirical model of basal metabolic rate for birds.

McKechnie, A. E., and B. O. Wolf. 2004. The Allometry of Avian Basal Metabolic Rate:
    Good Predictions Need Good Data. Physiological and Biochemical Zoology:
    Ecological and Evolutionary Approaches 77:502–521.
"""
struct McKechnieWolf <: MetabolicRateEquation end

function metabolic_rate(::McKechnieWolf, mass, body_temperature=nothing)
    mass_g = ustrip(u"g", mass)
    metabolic_flux = (10.0^(-1.461 + 0.669 * log10(mass_g)))u"W"

    return (metabolic_flux)
end

# conversion from O2 consumption to Joules

"""
    Typical <: OxygenJoulesConversion

Standard conversion between O₂ consumption and energy use.

# Arguments (Joules_to_O2) / Outputs (O2_to_Joules)
- `metabolic_flux` — metabolic rate (W).

# Outputs (Joules_to_O2) / Arguments (O2_to_Joules)
- V_O2_STP, oxygen consumption rate, standard temperature and pressure, L/s

TODO - add a ref?
"""
struct Typical <: OxygenJoulesConversion end

"""
    O2_to_Joules(conversion::OxygenJoulesConversion, V_O2_STP, rq)
    O2_to_Joules(V_O2_STP)

Convert oxygen consumption rate to metabolic heat rate.

# Arguments
- `conversion::OxygenJoulesConversion`: Conversion model (e.g., `Typical()`, `Kleiber1961()`)
- `V_O2_STP`: Oxygen consumption rate at STP
- `rq`: Respiratory quotient (CO₂ produced / O₂ consumed)

# Returns
- `metabolic_flux`: Metabolic heat rate (W)
"""
function O2_to_Joules(::Typical, V_O2_STP, rq)
    oxidation_flux = 20.1u"J/ml"
    metabolic_flux = V_O2_STP * oxidation_flux
    return (metabolic_flux)
end

"""
    Joules_to_O2(conversion::OxygenJoulesConversion, metabolic_flux, rq)
    Joules_to_O2(metabolic_flux)

Convert metabolic heat rate to oxygen consumption rate.

# Arguments
- `conversion::OxygenJoulesConversion`: Conversion model (e.g., `Typical()`, `Kleiber1961()`)
- `metabolic_flux`: Metabolic heat rate (W)
- `rq`: Respiratory quotient (CO₂ produced / O₂ consumed)

# Returns
- `V_O2_STP`: Oxygen consumption rate at STP
"""
function Joules_to_O2(::Typical, metabolic_flux, rq)
    oxidation_flux = 20.1u"J/ml"
    V_O2_STP = metabolic_flux / oxidation_flux
    return (V_O2_STP)
end

"""
    Kleiber1961 <: OxygenJoulesConversion

Kleiber's conversion between O₂ consumption and energy use, from
the bottom of Table 7.3.

# Arguments (Joules_to_O2) / Outputs (O2_to_Joules)
- `metabolic_flux` — metabolic rate (W).

# Outputs (Joules_to_O2) / Arguments (O2_to_Joules)
- V_O2_STP, oxygen consumption rate, standard temperature and pressure, L/s

Kleiber, M. 1961. The Fire of Life. An Introduction to Animal Energetics.
"""
struct Kleiber1961 <: OxygenJoulesConversion end

function Joules_to_O2(::Kleiber1961, metabolic_flux, rq)
    carbohydrate_oxidation_flux = u"J/ml"(5.0u"kcal"/1u"L")
    fat_oxidation_flux = u"J/ml"(4.7u"kcal"/1u"L")
    protein_oxidation_flux = u"J/ml"(4.5u"kcal"/1u"L")
    if rq ≥ 1.0
        V_O2_STP = metabolic_flux / carbohydrate_oxidation_flux
    end
    if rq ≤ 0.7
        V_O2_STP = metabolic_flux / fat_oxidation_flux
    else
        V_O2_STP = metabolic_flux / protein_oxidation_flux
    end
    return (V_O2_STP)
end

function O2_to_Joules(::Kleiber1961, V_O2_STP, rq)
    carbohydrate_oxidation_flux = u"J/ml"(5.0u"kcal"/1u"L")
    fat_oxidation_flux = u"J/ml"(4.7u"kcal"/1u"L")
    protein_oxidation_flux = u"J/ml"(4.5u"kcal"/1u"L")
    if rq ≥ 1.0
        metabolic_flux = V_O2_STP * carbohydrate_oxidation_flux
    end
    if rq ≤ 0.7
        metabolic_flux = V_O2_STP * fat_oxidation_flux
    else
        metabolic_flux = V_O2_STP * protein_oxidation_flux
    end
    return (metabolic_flux)
end


# Fallbacks: why are these needed ?
O2_to_Joules(::Missing) = missing
Joules_to_O2(::Missing) = missing

Joules_to_O2(metabolic_flux) = Joules_to_O2(Typical(), metabolic_flux, one(metabolic_flux))
O2_to_Joules(V_O2_STP) = O2_to_Joules(Typical(), V_O2_STP, one(V_O2_STP))

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
