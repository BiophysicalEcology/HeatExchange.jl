
# Heat balance


function energy_balance() end

function energy_balance(T_x)

    # compute areas for exchange
    A_v = A_tot * p_cond
    A_t = A_tot * p_cont
    A_sil = sil_area_of_cylinder(body_organism.geometry.lengths[2]/2, body_organism.geometry.lengths[1], Z)
    A_up = A_tot / 2
    A_down = A_tot / 2
    
    Q_solar = solar(α_org_dorsal, α_org_ventral, A_sil, A_up, A_down, α_sub, F_sub, F_sky, Q_dir, Q_dif)
    Q_IR_in = radin(A_tot, F_sky, F_sub, ϵ_org, ϵ_sub, ϵ_sky, T_sky, T_sub)
    Q_IR_out = radout(T_x, A_tot, F_sky, F_sub, ϵ_org)
    
    conv_out = convection(body_organism, A_v, T_air, T_x, vel, P_atmos, elev, fluid)
    evap_out = evap(T_x, T_x, J_resp, ψ_org, p_wet, A_tot, conv_out.Hd, p_eyes, T_air, rh, P_atmos)
    
    Q_conv = conv_out.Q_conv # convective heat loss
    Q_evap = evap_out.Q_evap # evaporative heat loss
    
    Q_in = Q_solar + Q_IR_in # energy in
    Q_out = Q_IR_out + Q_conv + Q_evap # energy out
    Q_in - Q_out # this must balance
    
end

# function energy_balance(T_x, b::Body, i::Naked, p::Model, o::OrganismalVars, e::EnvironmentalVars)

#     # compute areas for exchange
#     A_v = A_tot * p_cond
#     A_t = A_tot * p_cont
#     A_sil = sil_area_of_cylinder(body_organism.geometry.lengths[2]/2, body_organism.geometry.lengths[1], Z)
#     A_up = A_tot / 2
#     A_down = A_tot / 2
    
#     Q_solar = solar(α_org_dorsal, α_org_ventral, A_sil, A_up, A_down, α_sub, F_sub, F_sky, Q_dir, Q_dif)
#     Q_IR_in = radin(A_tot, F_sky, F_sub, ϵ_org, ϵ_sub, ϵ_sky, T_sky, T_sub)
#     Q_IR_out = radout(T_x, A_tot, F_sky, F_sub, ϵ_org)
    
#     conv_out = convection(body_organism, A_v, T_air, T_x, vel, P_atmos, elev, fluid)
#     evap_out = evap(T_x, T_x, J_resp, ψ_org, p_wet, A_tot, conv_out.Hd, p_eyes, T_air, rh, P_atmos)
    
#     Q_conv = conv_out.Q_conv # convective heat loss
#     Q_evap = evap_out.Q_evap # evaporative heat loss
    
#     Q_in = Q_solar + Q_IR_in # energy in
#     Q_out = Q_IR_out + Q_conv + Q_evap # energy out
#     Q_in - Q_out # this must balance
    
# end        