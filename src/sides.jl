"""
    sides(organism::Organism, environment; smoothing = HardBound())

The `dorsal` and `ventral` [`Part`](@ref)s by which a single body is solved, as in
NicheMapR: each is the whole body in that side's fur, seeing the sky (dorsal) or the ground
(ventral) over twice its share of the view and carrying that share of the body's heat.
"""
function sides(o::Organism, e; smoothing::SmoothingStrategy=HardBound())
    environment_pars = stripparams(e.environment_pars)
    environment_vars = e.environment_vars
    ins_pars = insulation_pars(o)
    (; conduction_fraction) = conduction_pars_external(o)
    internal_conduction = conduction_pars_internal(o)
    rad_pars = radiation_pars(o)
    evap_pars = evaporation_pars(o)
    view = _view_shares(rad_pars, environment_vars)

    solar_out = solar(o.body, Absorptivities(rad_pars, environment_pars), ViewFactors(view.sky, view.ground, 0.0, 0.0),
        SolarConditions(environment_vars), silhouette(o.body, rad_pars.solar_orientation),
        BiophysicalGeometry.total_area(o.body) * conduction_fraction)
    sunlit = solar_out.solar_flow > 0.0u"W"
    dorsal_solar = sunlit ? 2.0 * solar_out.solar_direct_flow + 2.0 * solar_out.solar_sky_flow : 0.0u"W"
    ventral_solar = sunlit ?
        (solar_out.solar_substrate_flow / (1.0 - view.sky - view.vegetation)) * (1.0 - 2.0 * conduction_fraction) : 0.0u"W"

    fibres = _fibres(ins_pars, rad_pars.ventral_fraction; smoothing)
    fat = FatLayer(internal_conduction.fat_fraction, internal_conduction.fat_density)
    temperature = EnvironmentTemperatures(
        environment_vars.air_temperature, environment_vars.sky_temperature,
        environment_vars.ground_temperature, environment_vars.reference_air_temperature,
        environment_vars.bush_temperature, environment_vars.substrate_temperature,
    )
    function side_part(side, view_factors, solar_flow, ϵ_body, weight)
        fibre = _side_value(fibres, side)
        body = Body(o.body.shape, CompositeInsulation(FibrousLayer(fibre.depth, fibre.diameter, fibre.density), fat))
        Part(;
            body,
            insulation = ins_pars,
            side,
            traits = (;
                flesh_conductivity = internal_conduction.flesh_conductivity,
                fat_conductivity = internal_conduction.fat_conductivity,
                ϵ_body,
                skin_wetness = evap_pars.skin_wetness,
                insulation_wetness = evap_pars.insulation_wetness,
                bare_skin_fraction = evap_pars.bare_skin_fraction,
                eye_fraction = evap_pars.eye_fraction,
            ),
            environment = (;
                temperature,
                view_factors,
                atmos = AtmosphericConditions(environment_vars),
                fluid = environment_pars.fluid,
                solar_flow,
                gas_fractions = environment_pars.gas_fractions,
                convection_enhancement = environment_pars.convection_enhancement,
            ),
            conduction_fraction,
            conductance_coefficient = _ground_conductance(side, body, conduction_fraction, environment_pars, environment_vars),
            ventral_fraction = rad_pars.ventral_fraction,
            characteristic_dimension = characteristic_dimension(convection_pars(o).characteristic_dimension_formula, body),
            weight,
        )
    end
    dorsal_share = view.sky + view.vegetation
    return (;
        dorsal = side_part(Dorsal(), ViewFactors(view.sky * 2.0, 0.0, 0.0, view.vegetation * 2.0),
            dorsal_solar, rad_pars.body_emissivity_dorsal, dorsal_share),
        ventral = side_part(Ventral(), ViewFactors(0.0, view.ground * 2.0, rad_pars.bush_view_factor * 2.0, 0.0),
            ventral_solar, rad_pars.body_emissivity_ventral, 1 - dorsal_share),
    )
end

# The shares of the view taken by open sky, by the vegetation shading it, and by the ground.
function _view_shares(rad_pars, environment_vars)
    vegetation = rad_pars.sky_view_factor * environment_vars.shade
    sky = rad_pars.sky_view_factor - vegetation
    (; sky, vegetation, ground = 1 - sky - vegetation)
end

# Only the belly touches the ground, over twice its share of the contact.
_ground_conductance(::Dorsal, body, conduction_fraction, environment_pars, environment_vars) = 0.0u"W/K"
_ground_conductance(::Ventral, body, conduction_fraction, environment_pars, environment_vars) =
    BiophysicalGeometry.total_area(body) * conduction_fraction * 2 * environment_vars.substrate_conductivity /
    environment_pars.conduction_depth

# The longwave a side emits, over twice its share of the view.
_side_lw_out(view_factor, emissivity, area, surface_temperature, σ) =
    2 * view_factor * emissivity * area * σ * surface_temperature^4

function _unpack_respiration(::Nothing)
    (; balance = nothing, respiration_heat_flow = 0.0u"W", respiration_mass_flow = nothing, air_flow = nothing,
        oxygen_flow_standard = nothing, molar_fluxes_in = nothing, molar_fluxes_out = nothing)
end
function _unpack_respiration(r)
    (; r.balance, r.respiration_heat_flow, r.respiration_mass_flow, r.air_flow, r.oxygen_flow_standard,
        r.molar_fluxes_in, r.molar_fluxes_out)
end

# The output of a single body from its solved `sides`, the flows of the two weighted by
# their shares.
function _organism_output(o, e, core_temperature, metabolic_heat_flow, respiration_out, sides, solved;
                          smoothing::SmoothingStrategy=HardBound())
    environment_vars = e.environment_vars
    ins_pars = insulation_pars(o)
    (; conduction_fraction) = conduction_pars_external(o)
    internal_conduction = conduction_pars_internal(o)
    rad_pars = radiation_pars(o)
    resp_pars = respiration_pars(o)
    metab_pars = metabolism_pars(o)
    evap_pars = evaporation_pars(o)
    view = _view_shares(rad_pars, environment_vars)
    dmult = sides.dorsal.weight
    vmult = sides.ventral.weight
    dorsal, ventral = solved.dorsal, solved.ventral
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    fibres = _fibres(ins_pars, rad_pars.ventral_fraction; smoothing)
    fat = FatLayer(internal_conduction.fat_fraction, internal_conduction.fat_density)
    geometry_avg = Body(o.body.shape, CompositeInsulation(
        FibrousLayer(fibres.average.depth, fibres.average.diameter, fibres.average.density), fat))
    total_area_final = BiophysicalGeometry.total_area(geometry_avg)
    area_skin_final = skin_area(geometry_avg)

    if _insulation_test(ins_pars) > zero(_insulation_test(ins_pars))
        surface_temperature_dorsal = dorsal.insulation_temperature
        surface_temperature_ventral = ventral.insulation_temperature
        area_rad_d = total_area_final
        area_rad_v = total_area_final * (1 - conduction_fraction)
    else
        surface_temperature_dorsal = dorsal.skin_temperature
        surface_temperature_ventral = ventral.skin_temperature
        area_rad_d = area_skin_final
        area_rad_v = area_skin_final
    end

    dorsal_radiation_out_flow = _side_lw_out(view.sky, rad_pars.body_emissivity_dorsal, area_rad_d, surface_temperature_dorsal, σ)
    ventral_radiation_out_flow = _side_lw_out(view.ground, rad_pars.body_emissivity_ventral, area_rad_v, surface_temperature_ventral, σ)
    dorsal_radiation_in_flow = dorsal_radiation_out_flow - dorsal.flows.longwave
    ventral_radiation_in_flow = ventral_radiation_out_flow - ventral.flows.longwave

    weighted(f) = f(dorsal) * dmult + f(ventral) * vmult
    solar_flow = weighted(s -> s.flows.solar)
    longwave_flow_in = dorsal_radiation_in_flow * dmult + ventral_radiation_in_flow * vmult
    longwave_flow_out = dorsal_radiation_out_flow * dmult + ventral_radiation_out_flow * vmult
    convection_heat_flow = weighted(s -> s.flows.convection)
    conduction_flow = weighted(s -> s.flows.conduction)
    evaporation_heat_flow =
        dorsal.flows.skin_evaporation * dmult + ventral.flows.skin_evaporation * vmult +
        dorsal.flows.insulation_evaporation * dmult + ventral.flows.insulation_evaporation * vmult

    (; balance, respiration_heat_flow, respiration_mass_flow, air_flow, oxygen_flow_standard,
       molar_fluxes_in, molar_fluxes_out) = _unpack_respiration(respiration_out)
    evaporation_heat_flow += respiration_heat_flow

    latent_heat_vaporisation = enthalpy_of_vaporisation(environment_vars.air_temperature)
    m_sweat = u"g/hr"((dorsal.flows.skin_evaporation + ventral.flows.skin_evaporation) * 0.5 / latent_heat_vaporisation)
    m_evap = isnothing(respiration_mass_flow) ? u"g/hr"(m_sweat) : u"g/hr"(respiration_mass_flow + m_sweat)

    skin_temperature = weighted(s -> s.skin_temperature)
    insulation_temperature = weighted(s -> s.insulation_temperature)
    insulation_final = insulation_properties(
        ins_pars, insulation_temperature * 0.7 + skin_temperature * 0.3, rad_pars.ventral_fraction; smoothing)

    heat_balance_val = solar_flow + longwave_flow_in + metabolic_heat_flow -
        longwave_flow_out - convection_heat_flow - evaporation_heat_flow - conduction_flow

    side_temperatures(s, depth) = (;
        skin_temperature = s.skin_temperature,
        insulation_temperature = s.insulation_temperature,
        insulation_conductivity = s.insulation_conductivity,
        insulation_depth = depth,
    )
    thermoregulation = (;
        core_temperature,
        skin_temperature,
        insulation_temperature,
        lung_temperature = (core_temperature + skin_temperature) * 0.5,
        insulation_depth = ins_pars.dorsal.depth * dmult + ins_pars.ventral.depth * vmult,
        axis_ratio_b = o.body.shape isa Sphere ? 1.0 : o.body.shape.axis_ratio_b,
        pant = resp_pars.pant,
        skin_wetness = evap_pars.skin_wetness,
        flesh_conductivity = internal_conduction.flesh_conductivity,
        insulation_conductivity_effective = insulation_final.conductivities.average,
        insulation_conductivity_compressed = insulation_final.conductivity_compressed,
        q10 = metab_pars.q10,
        dorsal = side_temperatures(dorsal, ins_pars.dorsal.depth),
        ventral = side_temperatures(ventral, ins_pars.ventral.depth),
    )

    morphology = (;
        total_area = total_area_final,
        area_skin = area_skin_final,
        area_evaporation = evaporation_area(o.body),
        area_convection = BiophysicalGeometry.total_area(o.body) * (1 - conduction_fraction),
        area_conduction = sides.ventral.area * conduction_fraction,
        area_silhouette = silhouette(geometry_avg, rad_pars.solar_orientation),
        sky_view_factor = view.sky,
        ground_view_factor = view.ground,
        volume = geometry_avg.geometry.volume,
        volume_flesh = flesh_volume(geometry_avg),
        characteristic_dimension = characteristic_dimension(VolumeCubeRoot(), geometry_avg),
        fat_mass = mass(geometry_avg.shape) * fat.fraction,
        geometry_avg.geometry.length...,
    )

    side_flows(s) = (;
        solar = s.flows.solar,
        convection = s.flows.convection,
        conduction = s.flows.conduction,
        net_metabolic = s.flows.net_metabolic,
        skin_evaporation = s.flows.skin_evaporation,
        insulation_evaporation = s.flows.insulation_evaporation,
        longwave = s.flows.longwave,
        sky_radiation = s.flows.sky_radiation,
        bush_radiation = s.flows.bush_radiation,
        vegetation_radiation = s.flows.vegetation_radiation,
        ground_radiation = s.flows.ground_radiation,
    )
    energy_flows = (;
        solar_flow,
        longwave_flow_in,
        metabolic_heat_flow,
        evaporation_heat_flow,
        respiration_heat_flow,
        longwave_flow_out,
        convection_heat_flow,
        conduction_flow,
        heat_balance = heat_balance_val,
        balance,
        ntry = max(dorsal.ntry, ventral.ntry),
        success = dorsal.success && ventral.success,
        dorsal = side_flows(dorsal),
        ventral = side_flows(ventral),
    )

    mass_flows = (;
        air_flow,
        oxygen_flow_standard,
        m_evap,
        respiration_mass_flow,
        m_sweat,
        molar_fluxes_in,
        molar_fluxes_out,
    )
    return ThermoregulationOutput(
        ThermoregulationState(thermoregulation),
        MorphologyState(morphology),
        EnergyFlowState(energy_flows),
        MassFlowState(mass_flows),
    )
end
