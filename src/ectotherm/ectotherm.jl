# ectotherm heat balance

"""
    ectotherm(core_temperature, organism::Organism, e)

Calculate heat balance for an ectotherm at a given core temperature.

Computes all heat flow components (solar, longwave, convection, evaporation,
conduction, respiration, metabolism) and returns the energy balance.

# Arguments
- `core_temperature`: Core temperature to evaluate
- `organism::Organism`: Organism with body geometry and traits
- `e`: Environment containing `environment_pars` and `environment_vars`

# Returns
NamedTuple with:
- `heat_balance`: Net heat balance (should be zero at equilibrium)
- `core_temperature`: Core temperature
- `surface_temperature`: Surface temperature
- `lung_temperature`: Lung temperature
- `enbal`: Energy balance components
- `masbal`: Mass balance components
- `respiration_out`, `solar_out`, `longwave_gain_out`, `longwave_loss_out`, `convection_out`, `evaporation_out`: Detailed outputs
"""
function ectotherm end

# A method dispatching on a `Model`
#ectotherm(core_temperature, mod::Model, e) = ectotherm(core_temperature, stripparams(mod),
#    stripparams(e.environment_pars), e.environment_vars)
# A generic method that expands dispatch to include the insulation
# this could be <:Ectotherm or <:Endotherm?
#ectotherm(core_temperature, o::Organism, environment_pars, vars) = ectotherm(core_temperature, insulation(o), o,
#    integumentpars(o), physiopars(o), thermoregpars(o), thermoregvars(o), environment_pars, vars)
# A method for Naked organisms
ectotherm(core_temperature, o::Organism, e) = ectotherm(core_temperature, insulation(body(o)), o, e)
function ectotherm(core_temperature, insulation::Naked, o::Organism, e)
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

    # calculate heat flows

    # metabolism
    metabolic_heat_flow = metabolic_rate(metab_pars.model, o.body.shape.mass, core_temperature)

    # respiration
    rates = MetabolicRates(; metabolic=metabolic_heat_flow)
    atmos = AtmosphericConditions(environment_vars)
    respiration_out = respiration(
        rates,
        resp_pars,
        atmos,
        o.body.shape.mass,
        core_temperature,  # lung_temperature
        environment_vars.air_temperature;
        gas_fractions=environment_pars.gas_fractions,
    )
    respiration_heat_flow = respiration_out.respiration_heat_flow
    oxygen_flow = u"ml/hr"(Joules_to_O2(metabolic_heat_flow))

    # net metabolic heat generation
    net_generated_heat_flow = metabolic_heat_flow - respiration_heat_flow
    specific_generated_heat_flow = net_generated_heat_flow / o.body.geometry.volume

    # resultant surface and lung temperature
    (; surface_temperature, lung_temperature) = surface_and_lung_temperature(;
        body=o.body, flesh_conductivity=internal_conduction.flesh_conductivity, specific_generated_heat_flow, core_temperature
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
    solar_flow = solar_out.solar_flow

    # longwave in
    emissivities = Emissivities(rad_pars, environment_pars)
    environmental_temps = EnvironmentTemperatures(environment_vars)
    longwave_gain_out = radiation_in(
        o.body,
        view_factors,
        emissivities,
        environmental_temps;
        conduction_fraction=external_conduction.conduction_fraction,
    )
    longwave_flow_in = longwave_gain_out.longwave_flow_in

    # longwave out
    longwave_loss_out = radiation_out(
        o.body,
        view_factors,
        emissivities,
        external_conduction.conduction_fraction,
        surface_temperature,  # dorsal_temperature
        surface_temperature,  # ventral_temperature
    )
    longwave_flow_out = longwave_loss_out.longwave_flow_out

    # conduction
    conduction_flow = conduction(;
        conduction_area,
        L=environment_pars.conduction_depth,
        surface_temperature,
        substrate_temperature=environment_vars.substrate_temperature,
        substrate_conductivity=environment_vars.substrate_conductivity,
    )

    # convection
    convection_out = convection(;
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
    convection_heat_flow = convection_out.heat_flow

    # evaporation
    evaporation_out = evaporation(
        evap_pars,
        convection_out.mass,
        atmos,
        convection_area,
        surface_temperature,
        environment_vars.air_temperature;
        water_potential=hyd_pars.water_potential,
        gas_fractions=environment_pars.gas_fractions,
    )
    evaporation_heat_flow = evaporation_out.evaporation_heat_flow

    # heat balance
    heat_flow_in = solar_flow + longwave_flow_in + metabolic_heat_flow # energy in
    heat_flow_out = longwave_flow_out + convection_heat_flow + evaporation_heat_flow + respiration_heat_flow + conduction_flow # energy out
    #@assert heat_flow_in - heat_flow_out = 0.0u"W" # this must balance
    heat_balance = heat_flow_in - heat_flow_out # this must balance

    enbal = (; solar_flow, longwave_flow_in, metabolic_heat_flow, respiration_heat_flow, evaporation_heat_flow, longwave_flow_out, convection_heat_flow, conduction_flow, heat_balance)
    masbal = (; oxygen_flow, respiration_mass=respiration_out.respiration_mass, cutaneous_mass=evaporation_out.m_cut, eye_mass=evaporation_out.m_eyes)
    (;
        heat_balance,
        core_temperature,
        surface_temperature,
        lung_temperature,
        enbal,
        masbal,
        respiration_out,
        solar_out,
        longwave_gain_out,
        longwave_loss_out,
        convection_out,
        evaporation_out,
    )
end
function ectotherm(core_temperature, insulation::Fur, pars, organism, vars) # A method for organisms with fur
    #....
end

"""
    get_Tb(mod::Model, environment_pars, vars)

Find the equilibrium core temperature for an ectotherm.

Uses root-finding (bisection) to find the core temperature where the heat
balance equals zero.

# Arguments
- `mod::Model`: Organism model
- `environment_pars`: Environmental parameters
- `vars`: Environmental variables

# Returns
Full ectotherm output at the equilibrium core temperature.
"""
function get_Tb(mod::Model, environment_pars, vars)
    air_temperature = vars.environment.air_temperature
    core_temperature = find_zero(
        t -> ectotherm(t, mod, environment_pars, vars), (air_temperature - 40K, air_temperature + 100K), Bisection()
    )
    ectotherm(core_temperature, mod, environment_pars, vars)
end

#flip2vectors(x) = (; (k => getfield.(x, k) for k in keys(x[1]))...)

#function ectotherm(core_temperature)

#    # compute areas for exchange
#    A_convection = A_total * (1 - conduction_fraction)
#    A_sil = silhouette_area(geometric_pars, zenith_angle)

#    # calculate heat flows
#    metab_out = metabolic_rate(geometric_pars.shape.mass, core_temperature, mass_normalisation, mass_exponent, thermal_sensitivity)
#    metabolic_heat_flow = metab_out.metabolic_heat_flow
#    resp_out = respiration_ectotherm(core_temperature, metabolic_heat_flow, fO2_extract, pant, rq, T_air, rh, elevation, P_atmos, fO2, fCO2, fN2)
#    respiration_heat_flow = resp_out.respiration_heat_flow
#    net_generated_heat_flow = metabolic_heat_flow - respiration_heat_flow
#    specific_generated_heat_flow = net_generated_heat_flow / geometric_pars.geometry.volume
#    (; surface_temperature, lung_temperature) = surface_and_lung_temperature(geometric_pars, k_flesh, specific_generated_heat_flow, core_temperature)
#    #Q_norm = direct_flow / cos(zenith_angle)
#    solar_out = solar(α_body_dorsal, α_body_ventral, A_sil, A_total, A_conduction, F_ground, F_sky, α_substrate, solar_flow, direct_flow, diffuse_flow)
#    solar_flow = solar_out.solar_flow
#    longwave_gain = radiation_in(A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral, ϵ_ground, ϵ_sky, T_sky, T_ground)
#    longwave_flow_in = longwave_gain.longwave_flow_in
#    longwave_loss = radiation_out(T_surface, A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral)
#    longwave_flow_out = longwave_loss.longwave_flow_out
#    conduction_flow = conduction(A_conduction, Le, T_surface, T_substrate, k_substrate)
#    conv_out = convection(; body=geometric_pars, A_convection, T_air, T_surface, wind_speed, P_atmos, fluid)
#    evap_out = evaporation(T_surface, ψ_org, skin_wetness, A_convection, conv_out.hd, eye_fraction, T_air, rh, P_atmos)
#    convection_heat_flow = conv_out.convection_heat_flow
#    evaporation_heat_flow = evap_out.evaporation_heat_flow

#    # calculate balance
#    heat_flow_in = solar_flow + longwave_flow_in + metabolic_heat_flow # energy in
#    heat_flow_out = longwave_flow_out + convection_heat_flow + evaporation_heat_flow + respiration_heat_flow + conduction_flow # energy out
#    #heat_flow_in - heat_flow_out # this must balance
#    heat_balance = heat_flow_in - heat_flow_out # this must balance

#    enbal = [solar_flow, longwave_flow_in, metabolic_heat_flow, respiration_heat_flow, evaporation_heat_flow, longwave_flow_out, convection_heat_flow, conduction_flow, heat_balance]
#    masbal = [metab_out.oxygen_flow, resp_out.respiration_mass, evap_out.m_cut, evap_out.m_eyes]
#    (; heat_balance, core_temperature, surface_temperature, lung_temperature, enbal, masbal, resp_out, solar_out, longwave_gain, longwave_loss, conv_out, evap_out)
#end
