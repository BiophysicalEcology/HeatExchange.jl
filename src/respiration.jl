"""
    respiration(rates, resp_pars, atmos, mass, lung_temperature, air_temperature; kw...)

Computes respiratory heat and water loss via mass flow through the lungs
given gas concentrations, pressure, respiration rate and humidity.

# Arguments
- `rates::MetabolicRates`: Metabolic rates (metabolic, sum, minimum)
- `resp_pars::RespirationParameters`: Respiration parameters (oxygen_extraction_efficiency, pant, respiratory_quotient, etc.)
- `atmos::AtmosphericConditions`: Atmospheric conditions (relative_humidity, atmospheric_pressure)
- `mass`: Body mass
- `lung_temperature`: Lung/core temperature
- `air_temperature`: Air temperature

# Keywords
- `gas_fractions::GasFractions`: Gas fractions for O2, CO2, N2 concentrations
- `O2conversion::OxygenJoulesConversion`: Model to convert O2 to Watts

# Returns
NamedTuple with balance, respiration_flux, m_resp, generated_flux, V_air, V_O2_STP, molar_fluxes
"""
function respiration(
    rates::MetabolicRates,
    resp_pars::RespirationParameters,
    atmos::AtmosphericConditions,
    mass,
    lung_temperature,
    air_temperature;
    gas_fractions::GasFractions=GasFractions(),
    O2conversion::OxygenJoulesConversion=Typical(),
)
    (; metabolic, sum, minimum) = rates
    metabolic_flux, flux_sum, minimum_flux = metabolic, sum, minimum
    (; oxygen_extraction_efficiency, pant, respiratory_quotient, exhaled_temperature_offset, exhaled_relative_humidity) = resp_pars
    exit_air_temperature = lung_temperature  # exhaled air at lung temperature (original default)
    (; relative_humidity, atmospheric_pressure) = atmos
    (; oxygen_fraction=fO2, carbon_dioxide_fraction=fCO2, nitrogen_fraction=fN2) = gas_fractions
    # adjust O2 to ensure sum to 1
    if fO2 + fCO2 + fN2 != 1
        fO2 = 1 - (fN2 + fCO2)
    end

    resp_gen = max(metabolic_flux, minimum_flux)

    #Joule_m3_O2 = 20.1e6u"J/m^3" # joules of energy dissipated per m3 O2 consumed at STP (enthalpy of combustion)
    V_O2_STP = u"m^3/s"(Joules_to_O2(O2conversion, resp_gen, respiratory_quotient))

    # converting stp -> vol. of O2 at animal lung temperature, atm. press.
    O2_partial_pressure = atmospheric_pressure * fO2
    V_O2 = (V_O2_STP * O2_partial_pressure / 273.15u"K") * (lung_temperature / O2_partial_pressure)
    #n = PV/RT (ideal gas law: number of moles from press,vol,temp)
    J_O2 = uconvert(u"mol/s", atmospheric_pressure * V_O2 / (Unitful.R * lung_temperature)) # mol O2 consumed
    # moles/s of O2, N2, dry air at entrance [air flow = f(O2 consumption)]
    J_O2_in = J_O2 / oxygen_extraction_efficiency # actual oxygen flow in (moles/s), accounting for efficiency of extraction
    J_N2_in = J_O2_in * (fN2 / fO2) #  actual nitrogen flow in (moles/s), accounting for efficiency of extraction
    V_air = V_O2 / fO2 # air flow
    V_CO2 = fCO2 * V_air #O2 flow
    J_CO2_in = atmospheric_pressure * V_CO2 / (Unitful.R * lung_temperature)
    J_air_in = (J_O2_in + J_N2_in + J_CO2_in) * pant
    V_air = uconvert(u"m^3/s", (J_air_in * Unitful.R * 273.15u"K" / 101325u"Pa")) # air volume @ stp (m3/s)
    # computing the vapor pressure at saturation for the subsequent calculation of
    # actual moles of water based on actual relative humidity
    saturation_vapour_pressure = vapour_pressure(air_temperature)
    J_H2O_in = J_air_in * (saturation_vapour_pressure * relative_humidity) / (atmospheric_pressure - saturation_vapour_pressure * relative_humidity)
    # moles at exit
    J_O2_out = J_O2_in - J_O2 # remove consumed oxygen from the total
    J_N2_out = J_N2_in
    J_CO2_out = respiratory_quotient * J_O2 + J_CO2_in
    # total moles of air at exit will be approximately the same as at entrance, since
    # the moles of O2 removed = approx. the # moles of co2 added
    J_air_out = (J_O2_out + J_N2_out + J_CO2_out) * pant
    # assuming saturated air at exit
    exit_vapour_pressure = vapour_pressure(exit_air_temperature)
    J_H2O_out = J_air_out * (exit_vapour_pressure / (atmospheric_pressure - exit_vapour_pressure))
    # enthalpy = U2-U1, internal energy only, i.e. lat. heat of vap. only involved, since assume
    # P,T,V constant, so not significant flow energy, PV. (H = U + PV)

    # moles/s lost by breathing:
    J_evap = J_H2O_out - J_H2O_in
    # grams/s lost by breathing = moles lost * gram molecular weight of water:
    m_resp = J_evap * 18u"g/mol"

    # putting a cap on water loss for small animals in very cold conditions
    # by assuming they will seek more moderate conditions if they exceed
    # this cap. this will improve stability for solution method.
    # based on data from w.r. welch. 1980. evaporative water loss from
    # endotherms in thermally and hygically complex environments: an
    # empirical approach for interspecific comparisons.
    # j. comp. physiol. 139: 135-143. maximum value recorded for prairie
    # dogs was 0.6 g/(kg-h) = 1.667 x 10**-4 g/(kg-s)
    # highest recorded rate was a resting deer mouse at 8 g/kg-h =
    # 2.22e-03*
    # (edwards & haines 1978. j. comp. physiol. 128: 177-184 in welch 1980)
    # for a 0.01 kg animal, the max. rate would be 1.67 x 10^-6 g/s
    m_resp = min(m_resp, (2.22E-03 * ustrip(u"kg", mass) * 15)u"g/s")

    # get latent heat of vapourisation and compute heat exchange due to respiration
    L_v = enthalpy_of_vaporisation(lung_temperature)
    # heat loss by breathing (J/s)=(J/kg)*(kg/s)
    respiration_flux = uconvert(u"W", L_v * m_resp)

    # get latent heat of vapourisation and compute heat exchange due to respiration
    L_v = enthalpy_of_vaporisation(lung_temperature)
    (; molar_mass) = dry_air_properties(air_temperature, atmospheric_pressure; gas_fractions)
    M_a = molar_mass
    (; specific_heat) = wet_air_properties(air_temperature, relative_humidity, atmospheric_pressure; gas_fractions)
    c_p = specific_heat
    air_flux = c_p * J_air_in * M_a * (air_temperature - lung_temperature)
    respiration_flux = uconvert(u"W", L_v * m_resp) - air_flux
    net_flux_check = metabolic_flux - respiration_flux
    balance = net_flux_check - flux_sum
    molar_fluxes = MolarFluxes(
        J_air_in,
        J_air_out,
        J_H2O_in,
        J_H2O_out,
        J_O2_in,
        J_O2_out,
        J_CO2_in,
        J_CO2_out,
        J_N2_in,
        J_N2_out,
    )
    return (; balance, respiration_flux, m_resp, generated_flux=metabolic_flux, V_air, V_O2_STP, molar_fluxes)
end
