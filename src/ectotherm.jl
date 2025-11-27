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
heat_balance(T_x, o::Organism, e_pars, vars) = heat_balance(T_x, insulation(o), o, morphopars(o), physiopars(o), e_pars, vars)
# A method for Naked organisms

function heat_balance(T_x, insulation::Naked, o, morphopars, physiopars, e_pars, vars)
    o_vars = vars.organism # make small function to get this
    e_vars = vars.environment # make small function to get this
    
    # compute areas for exchange
    A_total = get_total_area(o.body) #o.body.geometry.area.total
    A_convection = A_total * (1 - morphopars.conduction_fraction)
    A_conduction = A_total * morphopars.conduction_fraction
    #A_sil = silhouette_area(o.body, e_vars.zenith_angle)
    A_sil = silhouette_area(o.body.shape).normal # TODO make this work for zenith angle vs. normal/parallel

    # calculate heat fluxes
    metab_out = metabolism(o.body.shape.mass, T_x, physiopars.M1, physiopars.M2, physiopars.M3, physiopars.M4)
    Q_metab = metab_out.Q_metab
    resp_out = respiration(T_x, Q_metab, physiopars.fO2_extract, physiopars.pant, physiopars.rq, e_vars.T_air, e_vars.rh, e_vars.P_atmos, e_pars.fO2, e_pars.fCO2, e_pars.fN2)
    Q_resp = resp_out.Q_resp
    Q_gen_net = Q_metab - Q_resp
    Q_gen_spec = Q_gen_net / o.body.geometry.volume
    Tsurf_Tlung_out = Tsurf_and_Tlung(o.body, morphopars.k_body, Q_gen_spec, T_x)
    T_surface = Tsurf_Tlung_out.T_surface
    T_lung = Tsurf_Tlung_out.T_lung
    
    solar_out = solar(morphopars.α_body_dorsal, morphopars.α_body_ventral, A_sil, A_total, A_conduction, morphopars.F_substrate, morphopars.F_sky, e_pars.α_substrate, e_pars.shade, e_vars.solar_radiation, e_vars.direct_radiation, e_vars.diffuse_radiation)
    Q_solar = solar_out.Q_solar
    ir_gain = radin(A_total, A_conduction, morphopars.F_sky, morphopars.F_substrate, morphopars.ϵ_body_dorsal, morphopars.ϵ_body_ventral, e_pars.ϵ_substrate, e_pars.ϵ_sky, e_vars.T_sky, e_vars.T_substrate)
    Q_ir_in = ir_gain.Q_ir_in
    ir_loss = radout(T_surface, A_total, A_conduction, morphopars.F_sky, morphopars.F_substrate, morphopars.ϵ_body_dorsal, morphopars.ϵ_body_ventral)
    Q_ir_out = ir_loss.Q_ir_out
    Le = 0.025m
    Q_cond = conduction(A_conduction, Le, T_surface, e_vars.T_substrate, e_vars.k_substrate)
    conv_out = convection(o.body, A_convection, e_vars.T_air, T_surface, e_vars.wind_speed, e_vars.P_atmos, e_pars.fluid, e_pars.fO2, e_pars.fCO2, e_pars.fN2)
    evap_out = evaporation(T_surface, o_vars.ψ_org, morphopars.skin_wetness, A_convection, conv_out.hd, morphopars.eye_fraction, e_vars.T_air, e_vars.rh, e_vars.P_atmos, e_pars.fO2, e_pars.fCO2, e_pars.fN2)
    Q_conv = conv_out.Q_conv # convective heat loss
    Q_evap = evap_out.Q_evap # evaporative heat loss
    
    # calculate balance
    Q_in = Q_solar + Q_ir_in + Q_metab # energy in
    Q_out = Q_ir_out + Q_conv + Q_evap + Q_resp + Q_cond # energy out
    #Q_in - Q_out # this must balance
    Q_bal = Q_in - Q_out # this must balance

    enbal = (; Q_solar, Q_ir_in, Q_metab, Q_resp, Q_evap, Q_ir_out, Q_conv, Q_cond, Q_bal)
    #enbal = [Q_solar = Q_solar, Q_ir_in = Q_ir_in, Q_metab = Q_metab, Q_resp = Q_resp, Q_evap = Q_evap, Q_ir_out = Q_ir_out, Q_conv = Q_conv, Q_cond = Q_cond, Q_bal = Q_bal]
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
#    A_sil = silhouette_area(geometric_traits, zenith_angle)

#    # calculate heat fluxes
#    metab_out = metabolism(geometric_traits.shape.mass, T_x, M1, M2, M3)
#    Q_metab = metab_out.Q_metab
#    resp_out = respiration(T_x, Q_metab, fO2_extract, pant, rq, T_air, rh, elevation, P_atmos, fO2, fCO2, fN2)
#    Q_resp = resp_out.Q_resp
#    Q_gen_net = Q_metab - Q_resp
#    Q_gen_spec = Q_gen_net / geometric_traits.geometry.volume
#    Tsurf_Tlung_out = Tsurf_and_Tlung(geometric_traits, k_body, Q_gen_spec, T_x)
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
#    conv_out = convection(geometric_traits, A_convection, T_air, T_surface, wind_speed, P_atmos, fluid)
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


