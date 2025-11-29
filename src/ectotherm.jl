using Roots

# Heat balance

"""
    heat_balance(T, organism::Union{Model,Organism}, pars::AbstractEnvironmentalPars, vars::AbstractEnvironmentalVars)

Calculate heat balance for an organism at temperature T.
"""
function heat_balance end

# A method dispatching on a `Model`
heat_balance(T_x, mod::Model, e_pars, vars) = heat_balance(T_x, stripparams(mod), stripparams(e_pars), vars)
# A generic method that expands dispatch to include the insulation
# this could be <:Ectotherm or <:Endotherm?
heat_balance(T_x, o::Organism, e_pars, vars) = heat_balance(T_x, insulation(o), o, integumentpars(o), physiopars(o), thermoregpars(o), e_pars, vars)
# A method for Naked organisms

function heat_balance(T_x, insulation::Naked, o, integumentpars, physiopars, thermoregpars, e_pars, vars)
    o_vars = vars.organism # make small function to get this
    e_vars = vars.environment # make small function to get this
    
    # compute areas for exchange
    A_total = get_total_area(o.body)
    A_convection = A_total * (1 - integumentpars.conduction_fraction)
    A_conduction = A_total * integumentpars.conduction_fraction
    #A_sil = silhouette_area(o.body, e_vars.zenith_angle)
    A_sil = silhouette_area(o.body.shape).normal # TODO make this work for zenith angle vs. normal/parallel

    # calculate heat fluxes
    metab_out = metabolic_rate(AndrewsPough2(), o.body.shape.mass, T_x)
    Q_metab = metab_out.Q_metab
    resp_out = respiration_ectotherm(T_x, Q_metab, physiopars.fO2_extract, thermoregpars.pant, physiopars.rq, 
        e_vars.T_air, e_vars.rh, e_vars.P_atmos, e_pars.fO2, e_pars.fCO2, e_pars.fN2)
    Q_resp = resp_out.Q_resp
    Q_gen_net = Q_metab - Q_resp
    Q_gen_spec = Q_gen_net / o.body.geometry.volume
    Tsurf_Tlung_out = Tsurf_and_Tlung(o.body, physiopars.k_flesh, Q_gen_spec, T_x)
    T_surface = Tsurf_Tlung_out.T_surface
    T_lung = Tsurf_Tlung_out.T_lung
    
    solar_out = solar(integumentpars.α_body_dorsal, integumentpars.α_body_ventral, A_sil, 
        A_total, A_conduction, integumentpars.F_substrate, integumentpars.F_sky, 
        e_pars.α_substrate, e_pars.shade, e_vars.solar_radiation, e_vars.direct_radiation, 
        e_vars.diffuse_radiation)
    Q_solar = solar_out.Q_solar
    ir_gain = radin(A_total, A_conduction, integumentpars.F_sky, integumentpars.F_substrate, 
        integumentpars.ϵ_body_dorsal, integumentpars.ϵ_body_ventral, e_pars.ϵ_substrate, 
        e_pars.ϵ_sky, e_vars.T_sky, e_vars.T_substrate)
    Q_ir_in = ir_gain.Q_ir_in
    ir_loss = radout(T_surface, A_total, A_conduction, integumentpars.F_sky, integumentpars.F_substrate, 
        integumentpars.ϵ_body_dorsal, integumentpars.ϵ_body_ventral)
    Q_ir_out = ir_loss.Q_ir_out
    Le = 0.025u"m"
    Q_cond = conduction(A_conduction, Le, T_surface, e_vars.T_substrate, e_vars.k_substrate)
    conv_out = convection(o.body, A_convection, e_vars.T_air, T_surface, e_vars.wind_speed, e_vars.P_atmos, e_pars.fluid, e_pars.fO2, e_pars.fCO2, e_pars.fN2)
    evap_out = evaporation(T_surface, o_vars.ψ_org, integumentpars.skin_wetness, A_convection, conv_out.hd, conv_out.hd_free, integumentpars.eye_fraction, 1.0, e_vars.T_air, e_vars.rh, e_vars.P_atmos, e_pars.fO2, e_pars.fCO2, e_pars.fN2)
    Q_conv = conv_out.Q_conv # convective heat loss
    Q_evap = evap_out.Q_evap # evaporative heat loss
    
    # calculate balance
    Q_in = Q_solar + Q_ir_in + Q_metab # energy in
    Q_out = Q_ir_out + Q_conv + Q_evap + Q_resp + Q_cond # energy out
    #@assert Q_in - Q_out = 0.0u"W" # this must balance
    Q_bal = Q_in - Q_out # this must balance

    enbal = (; Q_solar, Q_ir_in, Q_metab, Q_resp, Q_evap, Q_ir_out, Q_conv, Q_cond, Q_bal)
    masbal = (; V_O2 = metab_out.V_O2, m_resp = resp_out.m_resp, m_cut = evap_out.m_cut, m_eye = evap_out.m_eyes)
    (;Q_bal, T_core=T_x, T_surface, T_lung, enbal, masbal, resp_out, solar_out, ir_gain, ir_loss, conv_out, evap_out)

end
function heat_balance(T_x, insulation::Fur, pars, organism, vars) # A method for organisms with fur
    #....
end

function get_Tb(mod::Model, e_pars, vars)
    T_air = vars.environment.T_air
    T_c = find_zero(t -> heat_balance(t, mod, e_pars, vars), (T_air - 40K, T_air + 100K), Bisection())
    heat_balance(T_c, mod, e_pars, vars)
end

flip2vectors(x) = (; (k => getfield.(x, k) for k in keys(x[1]))...)

#function heat_balance(T_x)

#    # compute areas for exchange
#    A_convection = A_total * (1 - conduction_fraction)
#    A_sil = silhouette_area(geometric_pars, zenith_angle)

#    # calculate heat fluxes
#    metab_out = metabolic_rate(geometric_pars.shape.mass, T_x, M1, M2, M3)
#    Q_metab = metab_out.Q_metab
#    resp_out = respiration_ectotherm(T_x, Q_metab, fO2_extract, pant, rq, T_air, rh, elevation, P_atmos, fO2, fCO2, fN2)
#    Q_resp = resp_out.Q_resp
#    Q_gen_net = Q_metab - Q_resp
#    Q_gen_spec = Q_gen_net / geometric_pars.geometry.volume
#    Tsurf_Tlung_out = Tsurf_and_Tlung(geometric_pars, k_flesh, Q_gen_spec, T_x)
#    T_surface = Tsurf_Tlung_out.T_surface
#    T_lung = Tsurf_Tlung_out.T_lung
#    #Q_norm = Q_dir / cos(zenith_angle)
#    solar_out = solar(α_body_dorsal, α_body_ventral, A_sil, A_total, A_conduction, F_substrate, F_sky, α_substrate, Q_sol, Q_dir, Q_dif)
#    Q_solar = solar_out.Q_solar
#    ir_gain = radin(A_total, A_conduction, F_sky, F_substrate, ϵ_body_dorsal, ϵ_body_ventral, ϵ_substrate, ϵ_sky, T_sky, T_substrate)
#    Q_ir_in = ir_gain.Q_ir_in
#    ir_loss = radout(T_surface, A_total, A_conduction, F_sky, F_substrate, ϵ_body_dorsal, ϵ_body_ventral)
#    Q_ir_out = ir_loss.Q_ir_out
#    Q_cond = conduction(A_conduction, Le, T_surface, T_substrate, k_substrate)
#    conv_out = convection(geometric_pars, A_convection, T_air, T_surface, wind_speed, P_atmos, fluid)
#    evap_out = evaporation(T_surface, ψ_org, skin_wetness, A_convection, conv_out.hd, eye_fraction, T_air, rh, P_atmos)
#    Q_conv = conv_out.Q_conv
#    Q_evap = evap_out.Q_evap

#    # calculate balance
#    Q_in = Q_solar + Q_ir_in + Q_metab # energy in
#    Q_out = Q_ir_out + Q_conv + Q_evap + Q_resp + Q_cond # energy out
#    #Q_in - Q_out # this must balance
#    Q_bal = Q_in - Q_out # this must balance

#    enbal = [Q_solar, Q_ir_in, Q_metab, Q_resp, Q_evap, Q_ir_out, Q_conv, Q_cond, Q_bal]
#    #enbal = [Q_solar = Q_solar, Q_ir_in = Q_ir_in, Q_metab = Q_metab, Q_resp = Q_resp, Q_evap = Q_evap, Q_ir_out = Q_ir_out, Q_conv = Q_conv, Q_cond = Q_cond, Q_bal = Q_bal]
#    masbal = [metab_out.V_O2, resp_out.m_resp, evap_out.m_cut, evap_out.m_eyes]
#    (;Q_bal, T_core=T_x, T_surface, T_lung, enbal, masbal, resp_out, solar_out, ir_gain, ir_loss, conv_out, evap_out)
#end

"""
    respiration_ectotherm(; kw...)
    respiration_ectotherm(T_x, Q_metab, fO2_extract, pant, rq, T_air, rh, elevation, P_atmos, fO2, fCO2, fN2) 

Computes respiratory heat and water loss via mass flow through the lungs 
given gas concentrations, pressure, respiration rate and humidity for an ectotherm.
Note that there is no recovery of heat or moisture assumed in the nose.
If barometric preassure is known, elevation will be ignored. Otherwise, 
if atmospheric pressure is unknown, elevation will be used to estimate it.

# Keywords
- `T_x`: current core temperature guess, K
- `Q_metab`: metabolic rate, W
- `fO2_extract`: extraction efficiency, fractional
- `pant`: multiplier on breathing rate due to panting, -
- `rq`: respiratory quotient, (mol CO2 / mol O2)
- `T_air`: air temperature, K
- `rh`: relative humidity, fractional
- `elevation`: elevation, m
- `P_atmos`: barometric pressure, Pa
- `fO2`; fractional O2 concentration in atmosphere, -
- `fCO2`; fractional CO2 concentration in atmosphere, -
- `fN2`; fractional N2 concentration in atmosphere, -
"""
function respiration_ectotherm(;
    T_x,
    Q_metab,
    fO2_extract = 0.2,
    pant = 1.0,
    rq = 0.8,
    T_air,
    rh,
    P_atmos = 101325u"Pa",
    fO2 = 0.2095,
    fCO2 = 0.000412,
    fN2 = 0.7902,
)
    return respiration_ectotherm(T_x, Q_metab, fO2_extract, pant, rq, T_air, rh, P_atmos, fO2, fCO2, fN2)
end
function respiration_ectotherm(T_x, Q_metab, fO2_extract, pant, rq, T_air, rh, P_atmos, fO2, fCO2, fN2)
    # adjust O2 to ensure sum to 1
    if fO2 + fCO2 + fN2 != 1
        fO2 = 1 - (fN2 + fCO2)
    end
    P_O2 = P_atmos * fO2
    Joule_m3_O2 = 20.1e6u"J/m^3" # joules of energy dissipated per m3 O2 consumed at STP (enthalpy of combustion)
    V_O2_STP = uconvert(u"m^3/s", Q_metab / Joule_m3_O2)

    # converting stp -> vol. of O2 at animal lung temperature, atm. press.
    T_lung = T_x
    V_O2 = (V_O2_STP * P_O2 / 273.15u"K") * (T_lung / P_O2)
    #n = PV/RT (ideal gas law: number of moles from press,vol,temp)
    J_O2 = uconvert(u"mol/s", P_atmos * V_O2 / (Unitful.R * T_x)) # mol O2 consumed
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
    #wet_air_out = wet_air_properties(T_air, rh, P_atmos; fO2, fCO2, fN2)
    P_vap_sat = vapour_pressure(T_air)
    J_H2O_in = J_air_in * (P_vap_sat * rh) / (P_atmos - P_vap_sat * rh)
    # moles at exit
    J_O2_out = J_O2_in - J_O2 # remove consumed oxygen from the total
    J_N2_out = J_N2_in
    J_CO2_out = rq * J_O2 + J_CO2_in
    # total moles of air at exit will be approximately the same as at entrance, since 
    # the moles of O2 removed = approx. the # moles of co2 added
    J_air_out = (J_O2_out + J_N2_out + J_CO2_out) * pant
    # assuming saturated air at exit
    P_vap_sat = vapour_pressure(T_x)
    J_H2O_out = J_air_out * (P_vap_sat / (P_atmos - P_vap_sat))
    # enthalpy = U2-U1, internal energy only, i.e. lat. heat of vap. only involved, since assume 
    # P,T,V constant, so not significant flow energy, PV. (H = U + PV)

    # moles/s lost by breathing:
    J_evap = J_H2O_out - J_H2O_in
    # grams/s lost by breathing = moles lost * gram molecular weight of water:
    m_resp = J_evap * 18u"g/mol"
    # get latent heat of vapourisation and compute heat exchange due to respiration
    #L_v = (2.5012e6 - 2.3787e3 * (Unitful.ustrip(T_lung) - 273.15))J / kg # from wet_air_properties
    L_v = enthalpy_of_vaporisation(T_lung)
    # heat loss by breathing (J/s)=(J/kg)*(kg/s)
    Q_resp = uconvert(u"W", L_v * m_resp)

    return (; Q_resp, m_resp, J_air_in, J_air_out, J_H2O_in, J_H2O_out, J_O2_in, J_O2_out, J_CO2_in, J_CO2_out)
end

Tsurf_and_Tlung(body::AbstractBody, k_flesh, Q_gen_spec, T_core) = Tsurf_and_Tlung(shape(body), body, k_flesh, Q_gen_spec, T_core)
function Tsurf_and_Tlung(shape::Cylinder, body, k_flesh, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.length[2]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_flesh)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_flesh) + T_surface 

    return (; T_surface, T_lung)  
end
function Tsurf_and_Tlung(shape::DesertIguana, body, k_flesh, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.length[1]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_flesh)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_flesh) + T_surface 

    return (; T_surface, T_lung)  
end
function Tsurf_and_Tlung(shape::LeopardFrog, body, k_flesh, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.length[1]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_flesh)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_flesh) + T_surface 

    return (; T_surface, T_lung)  
end
function Tsurf_and_Tlung(shape::Ellipsoid, body, k_flesh, Q_gen_spec, T_core)
    a = body.geometry.length[1] ^ 2
    b = body.geometry.length[2] ^ 2
    c = body.geometry.length[3] ^ 2
    x = ((a * b * c) / (a * b + a * c + b * c))
    T_surface = T_core - (Q_gen_spec / (2 * k_flesh)) * x
    T_lung = (Q_gen_spec / (4 * k_flesh)) * x + T_surface

    return (; T_surface, T_lung)  
end
