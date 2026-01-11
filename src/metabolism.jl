"""
    AndrewsPough2 <: MetabolicRateEquation

Empirical model of standard (fasted, inactive period) and resting (fasted, active period)
 metabolic rate of squamates, Eq. 2 in Andrews & Pough 1985. Physiol. Zool. 58:214-231.

# Arguments
- `mass` — body mass.
- `T_body` — body temperature.
- `M1` — mass normalisation constant.
- `M2` — mass allometric exponent. 
- `M3` — thermal senstivity constant. 
- `M4` — metabolic state (0 = standard, 1 = resting) 

# Outputs
- Q_metab, metabolic rate (W)
- V_O2, oxygen consumption rate (ml/h)

The equation reports metabolic rate as oxygen consumption rate given mass (g) and body
temperature (°C). The result is internally converted to watts (0.0056 conversion factor)
and both results are reported. Predictions are capped betwen 1 °C and 50 °C body temperature.
"""
struct AndrewsPough2 <: MetabolicRateEquation end
#TODO make it so M1 etc. are in AndrewPough2 struct, make them all take mass or mass and temp
function metabolic_rate(::AndrewsPough2, mass, T_body; M1=0.013, M2=0.8, M3=0.038, M4=0.0, 
        O2conversion::OxygenJoulesConversion=Typical())
    mass_g = ustrip(u"g", mass)
    T_c = ustrip(u"°C", T_body)

    if T_c > 1.0
        if T_c > 50.0
            V_O2 = M1 * mass_g^M2 * 10^(M3 * 50.0) * 10.0 ^ M4
         else
            V_O2 = M1 * mass_g^M2 * 10^(M3 * T_c) * 10.0 ^ M4
        end
    else
        V_O2 = M1 * mass_g^M2 * 10^(M3 * 1.0) * 10.0 ^ M4
    end
    
    Q_metab = u"W"(O2_to_Joules(O2conversion, (V_O2)u"ml/hr", 0.8))

    return (Q_metab)
end

"""
    Kleiber <: MetabolicRateEquation

Kleiber's empirical model of basal metabolic rate for mammals.

# Arguments
- `mass` — body mass.

# Outputs
- Q_metab, metabolic rate (W)

Kleiber, M. 1947. Body size and metabolic rate. Physiological Reviews 27:511–541.
"""
struct Kleiber <: MetabolicRateEquation end

function metabolic_rate(::Kleiber, mass)
    mass_kg = ustrip(u"kg", mass)
    m_rate_kcal_day = 70.0 * mass_kg^0.75 # general equation for mammals, p. 535, kcal/day
    Q_metab = (m_rate_kcal_day * 4.185 / (24 * 3.6))u"W" # convert to Watts
    return (Q_metab)
end

"""
    McKechnieWolf <: MetabolicRateEquation

Kleiber's empirical model of basal metabolic rate for birds.

# Arguments
- `mass` — body mass.

# Outputs
- Q_metab, metabolic rate (W)

McKechnie, A. E., and B. O. Wolf. 2004. The Allometry of Avian Basal Metabolic Rate: 
    Good Predictions Need Good Data. Physiological and Biochemical Zoology: 
    Ecological and Evolutionary Approaches 77:502–521.
"""
struct McKechnieWolf <: MetabolicRateEquation end

function metabolic_rate(::McKechnieWolf, mass)
    mass_g = ustrip(u"g", mass)
    Q_metab = (10.0^(-1.461 + 0.669 * log10(mass_g)))u"W"
    return (Q_metab)
end

# conversion from O2 consumption to Joules

"""
    Typical <: OxygenJoulesConversion

Standard conversion between O₂ consumption and energy use.

# Arguments (Joules_to_O2) / Outputs (O2_to_Joules)
- `Q_metab` — metabolic rate (W).

# Outputs (Joules_to_O2) / Arguments (O2_to_Joules)
- V_O2_STP, oxygen consumption rate, standard temperature and pressure, L/s

TODO - add a ref?
"""
struct Typical <: OxygenJoulesConversion end

function O2_to_Joules(::Typical, V_O2_STP, rq)
    Q_ox = 20.1u"J/ml"
    Q_metab = V_O2_STP * Q_ox
    return (Q_metab)
end

function Joules_to_O2(::Typical, Q_metab, rq)
    Q_ox = 20.1u"J/ml"
    V_O2_STP = Q_metab / Q_ox
    return (V_O2_STP)
end

"""
    Kleiber1961 <: OxygenJoulesConversion

Kleiber's conversion between O₂ consumption and energy use, from
the bottom of Table 7.3.

# Arguments (Joules_to_O2) / Outputs (O2_to_Joules)
- `Q_metab` — metabolic rate (W).

# Outputs (Joules_to_O2) / Arguments (O2_to_Joules)
- V_O2_STP, oxygen consumption rate, standard temperature and pressure, L/s

Kleiber, M. 1961. The Fire of Life. An Introduction to Animal Energetics.
"""
struct Kleiber1961 <: OxygenJoulesConversion end

function Joules_to_O2(::Kleiber1961, Q_metab, rq)
    Q_ox_carbohydrate = u"J/ml"(5.0u"kcal"/1u"L")
    Q_ox_fat = u"J/ml"(4.7u"kcal"/1u"L")
    Q_ox_protein = u"J/ml"(4.5u"kcal"/1u"L")
    if rq ≥ 1.0
        V_O2_STP = Q_metab / Q_ox_carbohydrate
    end
    if rq ≤ 0.7
        V_O2_STP = Q_metab / Q_ox_fat
    else
        V_O2_STP = Q_metab / Q_ox_protein
    end
    return (V_O2_STP)
end

function O2_to_Joules(::Kleiber1961, V_O2_STP, rq)
    Q_ox_carbohydrate = u"J/ml"(5.0u"kcal"/1u"L")
    Q_ox_fat = u"J/ml"(4.7u"kcal"/1u"L")
    Q_ox_protein = u"J/ml"(4.5u"kcal"/1u"L")
    if rq ≥ 1.0
        Q_metab = V_O2_STP * Q_ox_carbohydrate
    end
    if rq ≤ 0.7
        Q_metab = V_O2_STP * Q_ox_fat
    else
        Q_metab = V_O2_STP * Q_ox_protein
    end
    return (Q_metab)
end

# TODO
# """
#     ElliottDavison <: OxygenJoulesConversion

# Elliott and Davison conversion between O₂ consumption and energy use.

# # Arguments (Joules_to_O2) / Outputs (O2_to_Joules)
# - `Q_metab` — metabolic rate (W).

# # Outputs (Joules_to_O2) / Arguments (O2_to_Joules)
# - V_O2_STP, oxygen consumption rate, standard temperature and pressure, L/s

# Elliott, J. M., and W. Davison. 1975. Oecologia 19:195–201.
# """
# struct ElliottDavison <: OxygenJoulesConversion end

# function Joules_to_O2(::ElliottDavison, Q_metab, rq)
#     if rq ≥ 1.0
#         # carbohydrate metabolism
#         # carbohydrates ≈ 4193 cal/g
#         # L/s = (J/s) * (cal/J) * (kcal/cal) * (L O2 / kcal)
#         V_O2_STP = Q_metab * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/5.057u"kcal")
#     end
#     if rq ≤ 0.7
#         # fat metabolism; fats ≈ 9400 cal/g
#         V_O2_STP = Q_metab * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.7u"kcal")
#     else
#         # protein metabolism (RQ ≈ 0.8); ≈ 4300 cal/g
#         V_O2_STP = Q_metab * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.5u"kcal")
#     end
#     return (V_O2_STP)
# end

# function O2_to_Joules(::ElliottDavison, V_O2_STP, rq)
#     if rq ≥ 1.0
#         # carbohydrate metabolism
#         # carbohydrates ≈ 4193 cal/g
#         # L/s = (J/s) * (cal/J) * (kcal/cal) * (L O2 / kcal)
#         Q_metab  = V_O2_STP / ((1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/5.057u"kcal"))
#     end
#     if rq ≤ 0.7
#         # fat metabolism; fats ≈ 9400 cal/g
#         Q_metab  = V_O2_STP / ((1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.7u"kcal"))
#     else
#         # protein metabolism (RQ ≈ 0.8); ≈ 4300 cal/g
#         Q_metab  = V_O2_STP / ((1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.5u"kcal"))
#     end
#     return (Q_metab)
# end

"""
    vapour_pressure(T)
    vapour_pressure(formulation, T)

Calculates saturation vapour pressure (Pa) for a given air temperature.

# Arguments
- `T`: air temperature in K.

The `Typical` formulation is used by default.
"""
O2_to_Joules(::Missing) = missing
Joules_to_O2(::Missing) = missing

Joules_to_O2(Q_metab) = Joules_to_O2(Typical(), Q_metab, one(Q_metab))
O2_to_Joules(V_O2_STP) = O2_to_Joules(Typical(), V_O2_STP, one(V_O2_STP))