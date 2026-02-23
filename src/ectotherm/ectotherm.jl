# ectotherm heat balance

"""
    ectotherm(body_temperature, organism::Organism, e)

Calculate heat balance for an ectotherm at a given body temperature.

Computes all heat flux components (solar, infrared, convection, evaporation,
conduction, respiration, metabolism) and returns the energy balance.

# Arguments
- `body_temperature`: Body temperature to evaluate
- `organism::Organism`: Organism with body geometry and traits
- `e`: Environment containing `environment_pars` and `environment_vars`

# Returns
NamedTuple with:
- `heat_balance`: Net heat balance (should be zero at equilibrium)
- `core_temperature`: Core temperature (same as body_temperature for ectotherms)
- `surface_temperature`: Surface temperature
- `lung_temperature`: Lung temperature
- `enbal`: Energy balance components
- `masbal`: Mass balance components
- `resp_out`, `solar_out`, `ir_gain`, `ir_loss`, `conv`, `evap_out`: Detailed outputs
"""
function ectotherm end

# A method dispatching on a `Model`
#ectotherm(body_temperature, mod::Model, e) = ectotherm(body_temperature, stripparams(mod), 
#    stripparams(e.environment_pars), e.environment_vars)
# A generic method that expands dispatch to include the insulation
# this could be <:Ectotherm or <:Endotherm?
#ectotherm(body_temperature, o::Organism, environment_pars, vars) = ectotherm(body_temperature, insulation(o), o, 
#    integumentpars(o), physiopars(o), thermoregpars(o), thermoregvars(o), environment_pars, vars)
# A method for Naked organisms
ectotherm(body_temperature, o::Organism, e) = ectotherm(body_temperature, insulation(body(o)), o, e)
function ectotherm(body_temperature, insulation::Naked, o::Organism, e)
    environment_pars = stripparams(e.environment_pars) # TODO make small function to get this, or extract all?
    environment_vars = e.environment_vars # TODO make small function to get this, or extract all?
    external_conduction = conduction_pars_external(o)
    internal_conduction = conduction_pars_internal(o)
    #conv = convection_pars(o)
    rad_pars = radiation_pars(o)
    evap_pars = evaporation_pars(o)
    hyd_pars = hydraulic_pars(o)
    resp_pars = respiration_pars(o)
    metab_pars = metabolism_pars(o)

    # compute areas for exchange
    total_area = BiophysicalGeometry.total_area(o.body)
    convection_area = total_area * (1 - external_conduction.conduction_fraction)
    conduction_area = total_area * external_conduction.conduction_fraction
    silhouette_area = rad_pars.silhouette_area

    # calculate heat fluxes

    # metabolism
    metabolic_flux = metabolic_rate(metab_pars.model, o.body.shape.mass, body_temperature)

    # respiration
    rates = MetabolicRates(; metabolic=metabolic_flux)
    atmos = AtmosphericConditions(environment_vars)
    resp_out = respiration(
        rates,
        resp_pars,
        atmos,
        o.body.shape.mass,
        body_temperature,  # lung_temperature
        environment_vars.air_temperature;
        gas_fractions=environment_pars.gas_fractions,
    )
    respiration_flux = resp_out.respiration_flux
    oxygen_flow = u"ml/hr"(Joules_to_O2(metabolic_flux))

    # net metabolic heat generation
    net_generated_flux = metabolic_flux - respiration_flux
    generated_specific_flux = net_generated_flux / o.body.geometry.volume

    # resultant surface and lung temperature
    (; surface_temperature, lung_temperature) = surface_and_lung_temperature(;
        body=o.body, flesh_conductivity=internal_conduction.flesh_conductivity, generated_specific_flux, core_temperature=body_temperature
    )

    # solar radiation
    absorptivities = Absorptivities(rad_pars, environment_pars)
    view_factors = ViewFactors(rad_pars.sky_view_factor, rad_pars.ground_view_factor, 0.0, 0.0)
    solar_conditions = SolarConditions(environment_vars)
    solar_out = solar(
        o.body,
        absorptivities,
        view_factors,
        solar_conditions,
        silhouette_area,
        conduction_area,
    )
    solar_flux = solar_out.solar_flux

    # infrared in
    emissivities = Emissivities(rad_pars, environment_pars)
    env_temps = EnvironmentTemperatures(environment_vars)
    ir_gain = radin(
        o.body,
        view_factors,
        emissivities,
        env_temps;
        conduction_fraction=external_conduction.conduction_fraction,
    )
    infrared_in_flux = ir_gain.infrared_in_flux

    # infrared out
    ir_loss = radout(
        o.body,
        view_factors,
        emissivities,
        external_conduction.conduction_fraction,
        surface_temperature,  # dorsal_temperature
        surface_temperature,  # ventral_temperature
    )
    infrared_out_flux = ir_loss.infrared_out_flux

    # conduction
    conduction_flux = conduction(;
        conduction_area,
        L=environment_pars.conduction_depth,
        surface_temperature,
        substrate_temperature=environment_vars.substrate_temperature,
        substrate_conductivity=environment_vars.substrate_conductivity,
    )

    # convection
    conv = convection(;
        body=o.body,
        area=convection_area,
        air_temperature=environment_vars.air_temperature,
        surface_temperature,
        wind_speed=environment_vars.wind_speed,
        atmospheric_pressure=environment_vars.atmospheric_pressure,
        fluid=environment_pars.fluid,
        gas_fractions=environment_pars.gas_fractions,
        convection_enhancement=environment_pars.convection_enhancement,
    )
    convection_flux = conv.flux

    # evaporation
    evap_out = evaporation(
        evap_pars,
        conv.mass,
        atmos,
        convection_area,
        surface_temperature,
        environment_vars.air_temperature;
        water_potential=hyd_pars.water_potential,
        gas_fractions=environment_pars.gas_fractions,
    )
    evaporation_flux = evap_out.evaporation_flux

    # heat balance
    flux_in = solar_flux + infrared_in_flux + metabolic_flux # energy in
    flux_out = infrared_out_flux + convection_flux + evaporation_flux + respiration_flux + conduction_flux # energy out
    #@assert flux_in - flux_out = 0.0u"W" # this must balance
    heat_balance = flux_in - flux_out # this must balance

    enbal = (; solar_flux, infrared_in_flux, metabolic_flux, respiration_flux, evaporation_flux, infrared_out_flux, convection_flux, conduction_flux, heat_balance)
    masbal = (; oxygen_flow, respiration_mass=resp_out.respiration_mass, m_cut=evap_out.m_cut, m_eye=evap_out.m_eyes)
    (;
        heat_balance,
        core_temperature=body_temperature,
        surface_temperature,
        lung_temperature,
        enbal,
        masbal,
        resp_out,
        solar_out,
        ir_gain,
        ir_loss,
        conv,
        evap_out,
    )
end
function ectotherm(body_temperature, insulation::Fur, pars, organism, vars) # A method for organisms with fur
    #....
end

"""
    get_Tb(mod::Model, environment_pars, vars)

Find the equilibrium body temperature for an ectotherm.

Uses root-finding (bisection) to find the body temperature where the heat
balance equals zero.

# Arguments
- `mod::Model`: Organism model
- `environment_pars`: Environmental parameters
- `vars`: Environmental variables

# Returns
Full ectotherm output at the equilibrium body temperature.
"""
function get_Tb(mod::Model, environment_pars, vars)
    air_temperature = vars.environment.air_temperature
    body_temperature = find_zero(
        t -> ectotherm(t, mod, environment_pars, vars), (air_temperature - 40K, air_temperature + 100K), Bisection()
    )
    ectotherm(body_temperature, mod, environment_pars, vars)
end

#flip2vectors(x) = (; (k => getfield.(x, k) for k in keys(x[1]))...)

#function ectotherm(body_temperature)

#    # compute areas for exchange
#    A_convection = A_total * (1 - conduction_fraction)
#    A_sil = silhouette_area(geometric_pars, zenith_angle)

#    # calculate heat fluxes
#    metab_out = metabolic_rate(geometric_pars.shape.mass, body_temperature, mass_normalisation, mass_exponent, thermal_sensitivity)
#    Q_metab = metab_out.Q_metab
#    resp_out = respiration_ectotherm(body_temperature, Q_metab, fO2_extract, pant, rq, T_air, rh, elevation, P_atmos, fO2, fCO2, fN2)
#    Q_resp = resp_out.Q_resp
#    Q_gen_net = Q_metab - Q_resp
#    Q_gen_spec = Q_gen_net / geometric_pars.geometry.volume
#    (; surface_temperature, lung_temperature) = surface_and_lung_temperature(geometric_pars, k_flesh, Q_gen_spec, body_temperature)
#    #Q_norm = direct_flux / cos(zenith_angle)
#    solar_out = solar(α_body_dorsal, α_body_ventral, A_sil, A_total, A_conduction, F_ground, F_sky, α_substrate, solar_flux, direct_flux, diffuse_flux)
#    solar_fluxar = solar_out.solar_fluxar
#    ir_gain = radin(A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral, ϵ_ground, ϵ_sky, T_sky, T_ground)
#    infrared_in_flux = ir_gain.infrared_in_flux
#    ir_loss = radout(T_surface, A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral)
#    infrared_out_flux = ir_loss.infrared_out_flux
#    Q_cond = conduction(A_conduction, Le, T_surface, T_substrate, k_substrate)
#    conv_out = convection(; body=geometric_pars, A_convection, T_air, T_surface, wind_speed, P_atmos, fluid)
#    evap_out = evaporation(T_surface, ψ_org, skin_wetness, A_convection, conv_out.hd, eye_fraction, T_air, rh, P_atmos)
#    Q_conv = conv_out.Q_conv
#    Q_evap = evap_out.Q_evap

#    # calculate balance
#    Q_in = solar_fluxar + infrared_in_flux + Q_metab # energy in
#    Q_out = infrared_out_flux + Q_conv + Q_evap + Q_resp + Q_cond # energy out
#    #Q_in - Q_out # this must balance
#    Q_bal = Q_in - Q_out # this must balance

#    enbal = [solar_fluxar, infrared_in_flux, Q_metab, Q_resp, Q_evap, infrared_out_flux, Q_conv, Q_cond, Q_bal]
#    #enbal = [solar_fluxar = solar_fluxar, infrared_in_flux = infrared_in_flux, Q_metab = Q_metab, Q_resp = Q_resp, Q_evap = Q_evap, infrared_out_flux = infrared_out_flux, Q_conv = Q_conv, Q_cond = Q_cond, Q_bal = Q_bal]
#    masbal = [metab_out.V_O2, resp_out.respiration_mass, evap_out.m_cut, evap_out.m_eyes]
#    (;Q_bal, T_core=body_temperature, T_surface, T_lung, enbal, masbal, resp_out, solar_out, ir_gain, ir_loss, conv_out, evap_out)
#end
