function respiration_endotherm(x; T_air_reference, fO2, fN2, fCO2, P_atmos, Q_min, rq,
     T_lung, mass, fO2_extract, rh, rh_exit, T_air_exit, pant, Q_sum)

    T_air = T_air_reference
    gen = x
 
    # adjust O2 to ensure sum to 1
    if fO2 + fCO2 + fN2 != 1
        fO2 = 1 - (fN2 + fCO2)
    end

    resp_gen = max(gen, Q_min)

    # respgen is assumed to be in J/s
    # o2stp will be in L/s
    # litres of O₂/S @ STP: (data from Kleiber, 1961)
    if rq ≥ 1.0
        # carbohydrate metabolism
        # carbohydrates ≈ 4193 cal/g
        # L/s = (J/s) * (cal/J) * (kcal/cal) * (L O2 / kcal)
        V_O2_STP = resp_gen * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/5.057u"kcal")
    end

    if rq ≤ 0.7
        # fat metabolism; fats ≈ 9400 cal/g
        V_O2_STP = resp_gen * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.7u"kcal")
    else
        # protein metabolism (RQ ≈ 0.8); ≈ 4300 cal/g
        V_O2_STP = resp_gen * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.5u"kcal")
    end

    P_O2 = P_atmos * fO2
    V_O2 = (V_O2_STP * P_O2 / 273.15u"K") * (T_lung / P_O2)
    #n = PV/RT (ideal gas law: number of moles from press,vol,temp)
    J_O2 = uconvert(u"mol/s", P_atmos * V_O2 / (Unitful.R * T_lung)) # mol O2 consumed
    # moles/s of O2, N2, dry air at entrance [air flow = f(O2 consumption)]
    J_O2_in = J_O2 / fO2_extract # actual oxygen flow in (moles/s), accounting for efficiency of extraction
    J_N2_in = J_O2_in * (fN2 / fO2) #  actual nitrogen flow in (moles/s), accounting for efficiency of extraction
    V_air = V_O2 / fO2 # air flow
    V_CO2 = fCO2 * V_air #O2 flow
    J_CO2_in = P_atmos * V_CO2 / (Unitful.R * T_lung)
    J_air_in = (J_O2_in + J_N2_in + J_CO2_in) * pant
    V_air = uconvert(u"m^3/s", (J_air_in * Unitful.R * 273.15u"K" / 101325u"Pa")) # air volume @ stp (m3/s)
    # computing the vapor pressure at saturation for the subsequent calculation of 
    # actual moles of water based on actual relative humidity
    P_vap_sat = vapour_pressure(T_air)
    J_H2O_in = J_air_in * (P_vap_sat * rh) / (P_atmos - P_vap_sat * rh)
    # moles at exit
    J_O2_out = J_O2_in - J_O2 # remove consumed oxygen from the total
    J_N2_out = J_N2_in
    J_CO2_out = rq * J_O2 + J_CO2_in
    # total moles of air at exit will be approximately the same as at entrance, since 
    # the moles of O2 removed = approx. the # moles of co2 added
    J_air_out = (J_O2_out + J_N2_out + J_CO2_out) * pant
    # calculate vapour pressue of air at exit
    wet_air_out = wet_air_properties(T_air_exit, rh_exit, P_atmos; fO2, fCO2, fN2)
    P_vap_exit = wet_air_out.P_vap
    J_H2O_out = J_air_out * (P_vap_exit / (P_atmos - P_vap_exit))
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
    L_v = enthalpy_of_vaporisation(T_lung)
    #L_v = (2.5012e6 - 2.3787e3 * (Unitful.ustrip(u"°C",  T_lung)))u"J/kg" # from wet_air_properties
    # heat loss by breathing (J/s)=(J/kg)*(kg/s)
    (; M_a) = dry_air_properties(T_air, P_atmos)
    (; c_p) = wet_air_properties(T_air, rh, P_atmos; fO2, fCO2, fN2)
    Q_air = c_p * J_air_in * M_a * (T_air - T_lung)
    Q_resp = uconvert(u"W", L_v * m_resp) - Q_air
    # note that there is no recovery of heat or moisture assumed in the nose
    Q_net_check = x - Q_resp
    balance = Q_net_check - Q_sum
    return (; balance, Q_resp, m_resp, Q_gen = x, V_air, V_O2_STP, J_air_in, J_air_out, J_H2O_in, J_H2O_out, J_O2_in, J_O2_out, J_CO2_in, J_CO2_out, J_N2_in, J_N2_out)
end