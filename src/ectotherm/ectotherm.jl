# ectotherm heat balance

"""
    ectotherm(T_x, organism::Organism, e)

Calculate heat balance for an ectotherm at a given body temperature.

Computes all heat flux components (solar, infrared, convection, evaporation,
conduction, respiration, metabolism) and returns the energy balance.

# Arguments
- `T_x`: Body temperature to evaluate
- `organism::Organism`: Organism with body geometry and traits
- `e`: Environment containing `environment_pars` and `environment_vars`

# Returns
NamedTuple with:
- `Q_bal`: Net heat balance (should be zero at equilibrium)
- `T_core`: Core temperature (same as T_x for ectotherms)
- `T_surface`: Surface temperature
- `T_lung`: Lung temperature
- `enbal`: Energy balance components
- `masbal`: Mass balance components
- `resp_out`, `solar_out`, `ir_gain`, `ir_loss`, `conv_out`, `evap_out`: Detailed outputs
"""
function ectotherm end

# A method dispatching on a `Model`
#ectotherm(T_x, mod::Model, e) = ectotherm(T_x, stripparams(mod), 
#    stripparams(e.environment_pars), e.environment_vars)
# A generic method that expands dispatch to include the insulation
# this could be <:Ectotherm or <:Endotherm?
#ectotherm(T_x, o::Organism, e_pars, vars) = ectotherm(T_x, insulation(o), o, 
#    integumentpars(o), physiopars(o), thermoregpars(o), thermoregvars(o), e_pars, vars)
# A method for Naked organisms
ectotherm(T_x, o::Organism, e) = ectotherm(T_x, insulation(body(o)), o, e)
function ectotherm(T_x, insulation::Naked, o::Organism, e)
    e_pars = stripparams(e.environment_pars) # TODO make small function to get this, or extract all?
    e_vars = e.environment_vars # TODO make small function to get this, or extract all?
    cond_ex = conductionpars_external(o)
    cond_in = conductionpars_internal(o)
    #conv = convectionpars(o)
    rad = radiationpars(o)
    evap = evaporationpars(o)
    hyd = hydraulicpars(o)
    resp = respirationpars(o)
    metab = metabolismpars(o)

    # compute areas for exchange
    A_total = total_area(o.body)
    A_convection = A_total * (1 - cond_ex.conduction_fraction)
    A_conduction = A_total * cond_ex.conduction_fraction
    A_silhouette = rad.A_silhouette

    # calculate heat fluxes

    # metabolism
    Q_metab = metabolic_rate(metab.model, o.body.shape.mass, T_x)

    # respiration
    rates = MetabolicRates(; metabolic=Q_metab)
    atmos = AtmosphericConditions(e_vars)
    resp_out = respiration(
        rates,
        resp,
        atmos,
        o.body.shape.mass,
        T_x,  # T_lung
        e_vars.T_air;
        gasfrac=e_pars.gasfrac,
    )
    Q_resp = resp_out.Q_resp
    V_O2 = u"ml/hr"(Joules_to_O2(Q_metab))

    # net metabolic heat generation
    Q_gen_net = Q_metab - Q_resp
    Q_gen_spec = Q_gen_net / o.body.geometry.volume

    # resultant surfanec and lung temperature
    Tsurf_Tlung_out = Tsurf_and_Tlung(;
        body=o.body, k_flesh=cond_in.k_flesh, Q_gen_spec, T_core=T_x
    )
    T_surface = Tsurf_Tlung_out.T_surface
    T_lung = Tsurf_Tlung_out.T_lung

    # solar radiation
    absorptivities = Absorptivities(rad, e_pars)
    view_factors = ViewFactors(rad.F_sky, rad.F_ground, 0.0, 0.0)
    solar_conditions = SolarConditions(e_vars)
    solar_out = solar(
        o.body,
        absorptivities,
        view_factors,
        solar_conditions,
        A_silhouette,
        A_conduction,
    )
    Q_solar = solar_out.Q_solar

    # infrared in
    emissivities = Emissivities(rad, e_pars)
    env_temps = EnvironmentTemperatures(e_vars)
    ir_gain = radin(
        o.body,
        view_factors,
        emissivities,
        env_temps;
        conduction_fraction=cond_ex.conduction_fraction,
    )
    Q_ir_in = ir_gain.Q_ir_in

    # infrared out
    ir_loss = radout(
        o.body,
        view_factors,
        emissivities,
        cond_ex.conduction_fraction,
        T_surface,  # T_dorsal
        T_surface,  # T_ventral
    )
    Q_ir_out = ir_loss.Q_ir_out

    # conduction
    Q_cond = conduction(;
        A_conduction,
        L=e_pars.conduction_depth,
        T_surface,
        T_substrate=e_vars.T_substrate,
        k_substrate=e_vars.k_substrate,
    )

    # convection
    conv_out = convection(;
        body=o.body,
        area=A_convection,
        T_air=e_vars.T_air,
        T_surface,
        wind_speed=e_vars.wind_speed,
        P_atmos=e_vars.P_atmos,
        fluid=e_pars.fluid,
        gasfrac=e_pars.gasfrac,
        convection_enhancement=e_pars.convection_enhancement,
    )
    Q_conv = conv_out.Q_conv

    # evaporation
    transfer = TransferCoefficients(; heat=conv_out.hc, mass=conv_out.hd, mass_free=conv_out.hd_free)
    evap_out = evaporation(
        evap,
        transfer,
        atmos,
        A_convection,
        T_surface,
        e_vars.T_air;
        water_potential=hyd.water_potential,
        gasfrac=e_pars.gasfrac,
    )
    Q_evap = evap_out.Q_evap

    # heat balance
    Q_in = Q_solar + Q_ir_in + Q_metab # energy in
    Q_out = Q_ir_out + Q_conv + Q_evap + Q_resp + Q_cond # energy out
    #@assert Q_in - Q_out = 0.0u"W" # this must balance
    Q_bal = Q_in - Q_out # this must balance

    enbal = (; Q_solar, Q_ir_in, Q_metab, Q_resp, Q_evap, Q_ir_out, Q_conv, Q_cond, Q_bal)
    masbal = (; V_O2, m_resp=resp_out.m_resp, m_cut=evap_out.m_cut, m_eye=evap_out.m_eyes)
    (;
        Q_bal,
        T_core=T_x,
        T_surface,
        T_lung,
        enbal,
        masbal,
        resp_out,
        solar_out,
        ir_gain,
        ir_loss,
        conv_out,
        evap_out,
    )
end
function ectotherm(T_x, insulation::Fur, pars, organism, vars) # A method for organisms with fur
    #....
end

"""
    get_Tb(mod::Model, e_pars, vars)

Find the equilibrium body temperature for an ectotherm.

Uses root-finding (bisection) to find the body temperature where the heat
balance equals zero.

# Arguments
- `mod::Model`: Organism model
- `e_pars`: Environmental parameters
- `vars`: Environmental variables

# Returns
Full ectotherm output at the equilibrium body temperature.
"""
function get_Tb(mod::Model, e_pars, vars)
    T_air = vars.environment.T_air
    T_c = find_zero(
        t -> ectotherm(t, mod, e_pars, vars), (T_air - 40K, T_air + 100K), Bisection()
    )
    ectotherm(T_c, mod, e_pars, vars)
end

#flip2vectors(x) = (; (k => getfield.(x, k) for k in keys(x[1]))...)

#function ectotherm(T_x)

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
#    solar_out = solar(α_body_dorsal, α_body_ventral, A_sil, A_total, A_conduction, F_ground, F_sky, α_substrate, Q_sol, Q_dir, Q_dif)
#    Q_solar = solar_out.Q_solar
#    ir_gain = radin(A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral, ϵ_ground, ϵ_sky, T_sky, T_ground)
#    Q_ir_in = ir_gain.Q_ir_in
#    ir_loss = radout(T_surface, A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral)
#    Q_ir_out = ir_loss.Q_ir_out
#    Q_cond = conduction(A_conduction, Le, T_surface, T_substrate, k_substrate)
#    conv_out = convection(; body=geometric_pars, A_convection, T_air, T_surface, wind_speed, P_atmos, fluid)
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
