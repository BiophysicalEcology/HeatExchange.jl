"""
    metabolic_rate(model, mass, temperature)

Metabolic heat production. `model` is a `(variable, taxon)` pair for
`BiologicalScaling.allometric`, e.g. `(StandardMetabolicRate(), Squamate())` or
`(BasalMetabolicRate(), EutherianMammal())`, or `nothing` for none.
"""
metabolic_rate((variable, taxon)::Tuple{AbstractMetabolicRate,AbstractTaxon}, mass, temperature) =
    allometric(variable, taxon, mass, temperature)
metabolic_rate(::Nothing, mass, temperature) = 0.0u"W"

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
