
# Heat balance

function heat_balance() end
# A method dispatching on a `Model`
heat_balance(mod::Model, e_params, vars) = heat_balance(stripparams(mod), stripparams(e_params), vars)
# A generic method that expands dispatch to include the insulation
# this could be <:Ectotherm or <:Endotherm?
heat_balance(o::Organism, e_pars, vars) = heat_balance(insulation(o), o, params(o), e_pars, vars)
# A method for Naked organisms
function heat_balance(insulation::Naked, o, pars, e_pars, vars)
    o_vars = vars.organism # make small function to get this
    e_vars = vars.environment # make small function to get this
    
    T_x = o_vars.T_surf
    # compute areas for exchange
    A_tot = o.body.geometry.area
    A_c = A_tot * (1 - o_vars.p_cond)
    A_v = A_tot * o_vars.p_cond
    A_sil = calc_silhouette_area(o.body, e_vars.zen)
    A_up = A_tot / 2
    A_down = A_tot / 2

    Q_solar = solar(pars.α_org_dorsal, pars.α_org_ventral, A_sil, A_up, A_down, e_pars.α_sub, pars.F_sub, pars.F_sky, e_vars.Q_dir, e_vars.Q_dif)
    Q_IR_in = radin(A_tot, pars.F_sky, pars.F_sub, pars.ϵ_org_dorsal, pars.ϵ_org_ventral, e_pars.ϵ_sub, e_pars.ϵ_sky, e_vars.Tsky, e_vars.Tsub)
    Q_IR_out = radout(T_x, A_tot, pars.F_sky, pars.F_sub, pars.ϵ_org_dorsal, pars.ϵ_org_ventral)
    Q_metab = metabolism(o.body.shape.mass, o_vars.T_core, pars.M1, pars.M2, pars.M3)
    Le = 0.025m
    Q_cond = conduction(A_v, Le, o_vars.T_surf, e_vars.Tsub, e_vars.k_sub)
    
    conv_out = convection(o.body, A_c, e_vars.Ta, T_x, e_vars.vel, e_pars.P_atmos, e_pars.elev, e_pars.fluid)
    resp_out = respiration(T_x, Q_metab, pars.fO2_ext, o_vars.pant, pars.rq, e_vars.Ta, e_vars.rh, e_pars.elev, e_pars.P_atmos, e_pars.fO2, e_pars.fCO2, e_pars.fN2)
    m_resp = resp_out.m_resp
    evap_out = evaporation(T_x, o_vars.T_surf, m_resp, o_vars.ψ_org, o_vars.p_wet, A_tot, conv_out.Hd, pars.p_eyes, e_vars.Ta, e_vars.rh, e_pars.P_atmos)

    Q_conv = conv_out.Q_conv # convective heat loss
    Q_evap = evap_out.Q_evap # evaporative heat loss
    Q_resp = resp_out.Q_resp

    Q_in = Q_solar + Q_IR_in + Q_metab # energy in
    Q_out = Q_IR_out + Q_conv + Q_evap + Q_resp + Q_cond # energy out
    Q_in - Q_out # this must balance
end
function heat_balance(insulation::Fur, params, organism, vars) # A method for organisms with fur
    #....
end


function heat_balance(T_x)

    # compute areas for exchange
    A_c = A_tot * (1 - p_cond)
    A_sil = calc_silhouette_area(body_organism, Z)
    A_up = A_tot / 2
    A_down = A_tot / 2

    Q_solar = solar(α_org_dorsal, α_org_ventral, A_sil, A_up, A_down, α_sub, F_sub, F_sky, Q_dir, Q_dif)
    Q_IR_in = radin(A_tot, F_sky, F_sub, ϵ_org_dorsal, ϵ_org_ventral, ϵ_sub, ϵ_sky, T_sky, T_sub)
    Q_IR_out = radout(T_x, A_tot, F_sky, F_sub, ϵ_org_dorsal, ϵ_org_ventral)
    Q_metab = metabolism(body_organism.shape.mass, T_core, M1, M2, M3)
    Le = 0.025m
    Q_cond = conduction(A_v, Le, T_surf, T_sub, k_sub)
    
    conv_out = convection(body_organism, A_c, T_air, T_x, vel, P_atmos, elev, fluid)
    resp_out = respiration(T_x, Q_metab, fO2_ext, pant, rq, T_air, rh, elev, P_atmos, fO2, fCO2, fN2)
    m_resp = resp_out.m_resp
    evap_out = evaporation(T_x, T_surf, m_resp, ψ_org, p_wet, A_tot, conv_out.Hd, p_eyes, T_air, rh, P_atmos)

    Q_conv = conv_out.Q_conv # convective heat loss
    Q_evap = evap_out.Q_evap # evaporative heat loss

    Q_in = Q_solar + Q_IR_in # energy in
    Q_out = Q_IR_out + Q_conv + Q_evap # energy out
    Q_in - Q_out # this must balance

end