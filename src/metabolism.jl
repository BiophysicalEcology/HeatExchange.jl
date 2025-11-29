abstract type MetabolicRateEquation end

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

function metabolic_rate(::AndrewsPough2, mass, T_body; M1=0.013, M2=0.8, M3=0.038, M4=0.0)
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

    Q_metab = (0.0056 * V_O2)u"W"
    V_O2 = (V_O2)u"ml/hr"

    return (; Q_metab, V_O2)
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
    return (; Q_metab)
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
    mass_kg = ustrip(u"g", mass)
    Q_metab = 10^(-1.461 + 0.669 * log10(mass_kg))u"W"
    return (; Q_metab)
end
