
# Heat balance

function heat_balance() end
# A method dispatching on a `Model`
heat_balance(mod::Model, e_params, vars) = heat_balance(stripparams(mod), stripparams(e_params), vars)
heat_balance(T_x, mod::Model, e_params, vars) = heat_balance(T_x, stripparams(mod), stripparams(e_params), vars)
# A generic method that expands dispatch to include the insulation
# this could be <:Ectotherm or <:Endotherm?
heat_balance(T_x, o::Organism, e_pars, vars) = heat_balance(T_x, insulation(o), o, traits(o), e_pars, vars)
# A method for Naked organisms

function heat_balance(T_x, insulation::Naked, o, pars, e_pars, vars)
    o_vars = vars.organism # make small function to get this
    e_vars = vars.environment # make small function to get this
    
    # compute areas for exchange
    A_tot = o.body.geometry.area
    A_conv = A_tot * (1 - pars.p_cond)
    A_cond = A_tot * pars.p_cond
    A_sil = calc_silhouette_area(o.body, e_vars.zen)

    # calculate heat fluxes
    metab_out = metabolism(o.body.shape.mass, T_x, pars.M1, pars.M2, pars.M3)
    Q_metab = metab_out.Q_metab
    resp_out = respiration(T_x, Q_metab, pars.fO2_extract, pars.pant, pars.rq, e_vars.T_air, e_vars.rh, e_pars.elev, e_vars.P_atmos, e_pars.fO2, e_pars.fCO2, e_pars.fN2)
    Q_resp = resp_out.Q_resp
    Q_gen_net = Q_metab - Q_resp
    Q_gen_spec = Q_gen_net / o.body.geometry.volume
    Tsurf_Tlung_out = get_Tsurf_Tlung(o.body, pars.k_body, Q_gen_spec, T_x)
    T_surf = Tsurf_Tlung_out.T_surf
    
    #T_lung = Tsurf_Tlung_out.T_lung
    solar_out = solar(pars.α_org_dorsal, pars.α_org_ventral, A_sil, A_tot, A_cond, pars.F_sub, pars.F_sky, e_pars.α_sub, e_vars.Q_sol, e_vars.Q_dir, e_vars.Q_dif)
    Q_solar = solar_out.Q_solar
    ir_gain = radin(A_tot, A_cond, pars.F_sky, pars.F_sub, pars.ϵ_org_dorsal, pars.ϵ_org_ventral, e_pars.ϵ_sub, e_pars.ϵ_sky, e_vars.T_sky, e_vars.T_sub)
    Q_ir_in = ir_gain.Q_ir_in
    ir_loss = radout(T_surf, A_tot, A_cond, pars.F_sky, pars.F_sub, pars.ϵ_org_dorsal, pars.ϵ_org_ventral)
    Q_ir_out = ir_loss.Q_ir_out
    Le = 0.025m
    Q_cond = conduction(A_cond, Le, T_surf, e_vars.T_sub, e_vars.k_sub)
    conv_out = convection(o.body, A_conv, e_vars.T_air, T_surf, e_vars.vel, e_vars.P_atmos, e_pars.elev, e_pars.fluid)
    evap_out = evaporation(T_x, T_surf, resp_out.m_resp, o_vars.ψ_org, pars.p_wet, A_conv, conv_out.Hd, pars.p_eyes, e_vars.T_air, e_vars.rh, e_vars.P_atmos)
    Q_conv = conv_out.Q_conv # convective heat loss
    Q_evap = evap_out.Q_evap # evaporative heat loss
    
    # calculate balance
    Q_in = Q_solar + Q_ir_in + Q_metab # energy in
    Q_out = Q_ir_out + Q_conv + Q_evap + Q_resp + Q_cond # energy out
    #Q_in - Q_out # this must balance
    Q_bal = Q_in - Q_out # this must balance

    enbal = [Q_solar, Q_ir_in, Q_metab, Q_resp, Q_evap, Q_ir_out, Q_conv, Q_cond, Q_bal]
    #enbal = [Q_solar = Q_solar, Q_ir_in = Q_ir_in, Q_metab = Q_metab, Q_resp = Q_resp, Q_evap = Q_evap, Q_ir_out = Q_ir_out, Q_conv = Q_conv, Q_cond = Q_cond, Q_bal = Q_bal]
    masbal = [metab_out.V_O2, evap_out.m_resp, evap_out.m_cut, evap_out.m_eyes]
    (Q_bal=Q_bal, T_core=T_core, T_surf=T_surf, T_lung=T_lung, enbal=enbal, masbal=masbal, resp_out=resp_out, solar_out=solar_out, ir_gain=ir_gain, ir_loss=ir_loss, conv_out=conv_out, evap_out=evap_out)

end
function heat_balance(insulation::Fur, params, organism, vars) # A method for organisms with fur
    #....
end

function heat_balance(T_x)

    # compute areas for exchange
    A_conv = A_tot * (1 - p_cond)
    A_sil = calc_silhouette_area(geometric_traits, Z)

    # calculate heat fluxes
    metab_out = metabolism(geometric_traits.shape.mass, T_x, M1, M2, M3)
    Q_metab = metab_out.Q_metab
    resp_out = respiration(T_x, Q_metab, fO2_extract, pant, rq, T_air, rh, elev, P_atmos, fO2, fCO2, fN2)
    Q_resp = resp_out.Q_resp
    Q_gen_net = Q_metab - Q_resp
    Q_gen_spec = Q_gen_net / geometric_traits.geometry.volume
    Tsurf_Tlung_out = get_Tsurf_Tlung(geometric_traits, k_body, Q_gen_spec, T_x)
    T_surf = Tsurf_Tlung_out.T_surf
    #T_lung = Tsurf_Tlung_out.T_lung
    #Q_norm = Q_dir / cos(Z)
    solar_out = solar(α_org_dorsal, α_org_ventral, A_sil, A_tot, A_cond, F_sub, F_sky, α_sub, Q_sol, Q_dir, Q_dif)
    Q_solar = solar_out.Q_solar
    ir_gain = radin(A_tot, A_cond, F_sky, F_sub, ϵ_org_dorsal, ϵ_org_ventral, ϵ_sub, ϵ_sky, T_sky, T_sub)
    Q_ir_in = ir_gain.Q_ir_in
    ir_loss = radout(T_surf, A_tot, A_cond, F_sky, F_sub, ϵ_org_dorsal, ϵ_org_ventral)
    Q_ir_out = ir_loss.Q_ir_out
    Q_cond = conduction(A_cond, Le, T_surf, T_sub, k_sub)
    conv_out = convection(geometric_traits, A_conv, T_air, T_surf, vel, P_atmos, elev, fluid)
    evap_out = evaporation(T_x, T_surf, resp_out.m_resp, ψ_org, p_wet, A_conv, conv_out.Hd, p_eyes, T_air, rh, P_atmos)
    Q_conv = conv_out.Q_conv
    Q_evap = evap_out.Q_evap

    # calculate balance
    Q_in = Q_solar + Q_ir_in + Q_metab # energy in
    Q_out = Q_ir_out + Q_conv + Q_evap + Q_resp + Q_cond # energy out
    #Q_in - Q_out # this must balance
    Q_bal = Q_in - Q_out # this must balance

    enbal = [Q_solar, Q_ir_in, Q_metab, Q_resp, Q_evap, Q_ir_out, Q_conv, Q_cond, Q_bal]
    #enbal = [Q_solar = Q_solar, Q_ir_in = Q_ir_in, Q_metab = Q_metab, Q_resp = Q_resp, Q_evap = Q_evap, Q_ir_out = Q_ir_out, Q_conv = Q_conv, Q_cond = Q_cond, Q_bal = Q_bal]
    masbal = [metab_out.V_O2, evap_out.m_resp, evap_out.m_cut, evap_out.m_eyes]
    (Q_bal=Q_bal, T_core=T_core, T_surf=T_surf, T_lung=T_lung, enbal=enbal, masbal=masbal, resp_out=resp_out, solar_out=solar_out, ir_gain=ir_gain, ir_loss=ir_loss, conv_out=conv_out, evap_out=evap_out)
end
