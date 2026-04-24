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
NamedTuple with balance, respiration_heat_flow, respiration_mass, generated_heat_flow, air_flow, oxygen_flow_standard, molar_fluxes
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
    metabolic_heat_flow, heat_flow_sum, minimum_heat_flow = metabolic, sum, minimum
    (; oxygen_extraction_efficiency, pant, respiratory_quotient, exhaled_relative_humidity) = resp_pars
    exit_air_temperature = lung_temperature
    exhaled_rh = exhaled_relative_humidity[]
    (; relative_humidity, atmospheric_pressure) = atmos
    fO2 = gas_fractions.oxygen
    fCO2 = gas_fractions.carbon_dioxide
    fN2 = gas_fractions.nitrogen
    # adjust O2 to ensure sum to 1
    if !isapprox(fO2 + fCO2 + fN2, 1; atol=1e-6)
        fO2 = 1 - (fN2 + fCO2)
    end

    resp_gen = max(metabolic_heat_flow, minimum_heat_flow)

    # joules of energy dissipated per m3 O2 consumed at standard temperature and pressure (enthalpy of combustion)
    oxygen_flow_standard = u"m^3/s"(Joules_to_O2(O2conversion, resp_gen, respiratory_quotient))

    # converting standard temperature and pressure -> vol. of O2 at animal lung temperature, atm. press.
    O2_partial_pressure = atmospheric_pressure * fO2
    O2_ref = 101325u"Pa" * fO2
    oxygen_flow = oxygen_flow_standard * (O2_ref / 273.15u"K") * (lung_temperature / O2_partial_pressure)    #n = PV/RT (ideal gas law: number of moles from press,vol,temp)
    oxygen_consumed = uconvert(u"mol/s", atmospheric_pressure * oxygen_flow / (Unitful.R * lung_temperature)) # mol O2 consumed
    # moles/s of O2, N2, dry air at entrance [air flow = f(O2 consumption)]
    oxygen_in = oxygen_consumed / oxygen_extraction_efficiency # actual oxygen flow in (moles/s), accounting for efficiency of extraction
    nitrogen_in = oxygen_in * (fN2 / fO2) #  actual nitrogen flow in (moles/s), accounting for efficiency of extraction
    air_flow_raw = oxygen_flow / fO2 # air flow
    co2_flow = fCO2 * air_flow_raw # CO2 flow
    carbon_dioxide_in = atmospheric_pressure * co2_flow / (Unitful.R * lung_temperature)
    air_in = (oxygen_in + nitrogen_in + carbon_dioxide_in) * pant
    air_flow = uconvert(u"m^3/s", (air_in * Unitful.R * 273.15u"K" / 101325u"Pa")) # air volume at standard temperature and pressure (m3/s)
    # computing the vapor pressure at saturation for the subsequent calculation of
    # actual moles of water based on actual relative humidity
    saturation_vapour_pressure = vapour_pressure(air_temperature)
    water_in = air_in * (saturation_vapour_pressure * relative_humidity) / 
        max(atmospheric_pressure - saturation_vapour_pressure * relative_humidity, 1e-6u"Pa")
    # moles at exit
    oxygen_out = oxygen_in - oxygen_consumed # remove consumed oxygen from the total
    nitrogen_out = nitrogen_in
    carbon_dioxide_out = respiratory_quotient * oxygen_consumed + carbon_dioxide_in
    # total moles of air at exit will be approximately the same as at entrance, since
    # the moles of O2 removed = approx. the # moles of co2 added
    air_out = (oxygen_out + nitrogen_out + carbon_dioxide_out) * pant
    exit_vapour_pressure = vapour_pressure(exit_air_temperature) * exhaled_rh
    water_out = air_out * (exit_vapour_pressure / (atmospheric_pressure - exit_vapour_pressure))
    # enthalpy = U2-U1, internal energy only, i.e. lat. heat of vap. only involved, since assume
    # P,T,V constant, so not significant flow energy, PV. (H = U + PV)

    # moles/s lost by breathing:
    water_evaporated = water_out - water_in
    # grams/s lost by breathing = moles lost * gram molecular weight of water:
    respiration_mass = water_evaporated * 18u"g/mol"

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
    respiration_mass = min(respiration_mass, (2.22E-03 * ustrip(u"kg", mass) * 15)u"g/s")

    # get latent heat of vapourisation and compute heat exchange due to respiration
    latent_heat_vaporisation = enthalpy_of_vaporisation(lung_temperature)
    (; molar_mass) = dry_air_properties(air_temperature, atmospheric_pressure; gas_fractions)
    molar_mass_air = molar_mass
    (; specific_heat) = wet_air_properties(air_temperature, relative_humidity, atmospheric_pressure; gas_fractions)
    specific_heat_capacity = specific_heat
    sensible_heat_flow = specific_heat_capacity * air_in * molar_mass_air * (air_temperature - lung_temperature)
    respiration_heat_flow = uconvert(u"W", latent_heat_vaporisation * respiration_mass) - sensible_heat_flow
    net_heat_flow_check = metabolic_heat_flow - respiration_heat_flow
    balance = net_heat_flow_check - heat_flow_sum
    molar_fluxes_in = MolarFluxes(air_in, water_in, oxygen_in, carbon_dioxide_in, nitrogen_in)
    molar_fluxes_out = MolarFluxes(air_out, water_out, oxygen_out, carbon_dioxide_out, nitrogen_out)
    return (; balance, respiration_heat_flow, respiration_mass, generated_heat_flow=metabolic_heat_flow, air_flow, oxygen_flow_standard, molar_fluxes_in, molar_fluxes_out)
end
