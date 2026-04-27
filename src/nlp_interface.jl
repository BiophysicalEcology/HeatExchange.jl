# NLP physics interface for IPOPT and other nonlinear solvers.
#
# Exposes three functions so that the calling package (e.g. BiophysicalBehaviour.jl)
# contains only solver policy (bounds, objective, variable unpacking) and no physics:
#
#   nlp_pack      — pack all fixed physics parameters once before the NLP solve
#   nlp_residuals — evaluate heat-balance constraint residuals on every NLP iteration
#   nlp_assemble_output — build the full output NamedTuple after the NLP solve
#
# Two strategies dispatch the above:
#   WeightedMeanNLP — dorsal/ventral mean-weighted single body (3 equality residuals)
#   MultiSidedNLP   — explicit per-side heat balance (6 equality residuals)

abstract type NLPStrategy end
struct WeightedMeanNLP <: NLPStrategy end
struct MultiSidedNLP   <: NLPStrategy end

struct WeightedMeanNLPPacked{P}
    p::P
end

struct MultiSidedNLPPacked{P}
    p::P
end

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

# Rebuild a trial body from piloerection (insulation_depth) and postural (shape_b) effectors.
function _nlp_rebuild_side_body(base_shape, body_is_sphere, fat, fibre_props, insulation_depth, shape_b)
    trial_shape = body_is_sphere ? base_shape : setproperties(base_shape; b = shape_b)
    trial_fur   = Fur(insulation_depth, fibre_props.diameter, fibre_props.density)
    return Body(trial_shape, CompositeInsulation(trial_fur, fat))
end

# ---------------------------------------------------------------------------
# nlp_pack — WeightedMeanNLP
# ---------------------------------------------------------------------------

"""
    nlp_pack(::WeightedMeanNLP, organism, environment, T_skin_init, T_ins_init)

Pack all fixed physics parameters for the weighted-mean single-body NLP formulation.
Mirrors the pre-IPOPT setup in `thermoregulate` (lines 64–175 of the legacy
BiophysicalBehaviour.jl code) using only HeatExchange functions.

Returns a `WeightedMeanNLPPacked` whose `p` field carries the packed parameters.
"""
function nlp_pack(::WeightedMeanNLP, organism::Organism, environment, T_skin_init, T_ins_init)
    env_pars = stripparams(environment.environment_pars)
    env_vars = environment.environment_vars

    ins_pars   = insulation_pars(organism)
    ext_cond   = conduction_pars_external(organism)
    int_cond   = conduction_pars_internal(organism)
    rad_pars   = radiation_pars(organism)
    evap_pars  = evaporation_pars(organism)
    resp_pars  = respiration_pars(organism)

    T_air        = env_vars.air_temperature
    T_vegetation = env_vars.reference_air_temperature

    # View-factor geometry (mirrors _pack_sides and solve_metabolic_rate)
    vegetation_view_factor = rad_pars.sky_view_factor * env_vars.shade
    sky_view_factor        = rad_pars.sky_view_factor - vegetation_view_factor
    ground_view_factor     = 1 - sky_view_factor - vegetation_view_factor
    bush_view_factor       = rad_pars.bush_view_factor
    dorsal_weight          = sky_view_factor + vegetation_view_factor
    ventral_weight         = 1 - dorsal_weight
    ventral_fraction       = rad_pars.ventral_fraction

    # Mean-weighted insulation geometry
    fat = Fat(int_cond.fat_fraction, int_cond.fat_density)

    mean_depth       = ins_pars.dorsal.depth       * dorsal_weight + ins_pars.ventral.depth       * ventral_weight
    mean_diameter    = ins_pars.dorsal.diameter    * dorsal_weight + ins_pars.ventral.diameter    * ventral_weight
    mean_density     = ins_pars.dorsal.density     * dorsal_weight + ins_pars.ventral.density     * ventral_weight
    mean_length      = ins_pars.dorsal.length      * dorsal_weight + ins_pars.ventral.length      * ventral_weight
    mean_reflectance = ins_pars.dorsal.reflectance * dorsal_weight + ins_pars.ventral.reflectance * ventral_weight
    mean_conductivity = ins_pars.dorsal.conductivity * dorsal_weight + ins_pars.ventral.conductivity * ventral_weight

    mean_fibre_props = FibreProperties(;
        diameter     = mean_diameter,
        length       = mean_length,
        density      = mean_density,
        depth        = mean_depth,
        reflectance  = mean_reflectance,
        conductivity = mean_conductivity,
    )
    mean_ins_pars = InsulationParameters(;
        dorsal                  = mean_fibre_props,
        ventral                 = mean_fibre_props,
        depth_compressed        = ins_pars.depth_compressed,
        longwave_depth_fraction = ins_pars.longwave_depth_fraction,
    )
    mean_fur  = Fur(mean_depth, mean_diameter, mean_density)
    mean_body = Body(organism.body.shape, CompositeInsulation(mean_fur, fat))

    # Solar heat flows
    reference_total_area      = BiophysicalGeometry.total_area(organism.body)
    reference_conduction_area = reference_total_area * ext_cond.conduction_fraction
    absorptivities            = Absorptivities(rad_pars, env_pars)
    solar_view_factors        = ViewFactors(sky_view_factor, ground_view_factor, 0.0, 0.0)
    solar_conds               = SolarConditions(env_vars)
    reference_silhouette_area = silhouette_area(organism.body, rad_pars.solar_orientation)
    solar_result              = solar(organism.body, absorptivities, solar_view_factors,
                                     solar_conds, reference_silhouette_area, reference_conduction_area)

    if solar_result.solar_flow > 0.0u"W"
        dorsal_solar_flow  = 2.0 * solar_result.direct_flow + solar_result.solar_sky_flow * 2.0
        ventral_solar_flow = (solar_result.solar_substrate_flow /
                              (1.0 - sky_view_factor - vegetation_view_factor)) *
                             (1.0 - 2.0 * ext_cond.conduction_fraction)
    else
        dorsal_solar_flow  = 0.0u"W"
        ventral_solar_flow = 0.0u"W"
    end
    mean_solar_flow = dorsal_solar_flow * dorsal_weight + ventral_solar_flow * ventral_weight

    # Mean view factors: dorsal sees sky+veg, ventral sees ground+bush
    mean_view_factors = ViewFactors(
        sky_view_factor        * 2.0 * dorsal_weight,
        ground_view_factor     * 2.0 * ventral_weight,
        bush_view_factor       * 2.0 * ventral_weight,
        vegetation_view_factor * 2.0,
    )

    # Mean substrate conductance (ventral side only)
    mean_body_total_area     = BiophysicalGeometry.total_area(mean_body)
    ventral_conduction_area  = mean_body_total_area * ext_cond.conduction_fraction * 2
    mean_conduction_coeff    = (ventral_conduction_area * env_vars.substrate_conductivity) /
                               env_pars.conduction_depth * ventral_weight

    # Pack environment and traits for heat_balance
    mean_ins_temperature_init = T_ins_init * 0.7 + T_skin_init * 0.3
    insulation_props_init     = insulation_properties(mean_ins_pars, mean_ins_temperature_init, ventral_fraction)

    env_temperatures = EnvironmentTemperatures(
        T_air, env_vars.sky_temperature, env_vars.ground_temperature,
        T_vegetation, env_vars.bush_temperature, env_vars.substrate_temperature,
    )
    atmos = AtmosphericConditions(env_vars)
    heat_balance_env = (;
        temperature            = env_temperatures,
        view_factors           = mean_view_factors,
        atmos,
        fluid                  = env_pars.fluid,
        solar_flow             = mean_solar_flow,
        gas_fractions          = env_pars.gas_fractions,
        convection_enhancement = env_pars.convection_enhancement,
    )
    mean_body_emissivity = rad_pars.body_emissivity_dorsal * dorsal_weight +
                           rad_pars.body_emissivity_ventral * ventral_weight
    heat_balance_traits = (;
        fat_conductivity   = int_cond.fat_conductivity,
        flesh_conductivity = int_cond.flesh_conductivity,
        ϵ_body             = mean_body_emissivity,
        skin_wetness       = evap_pars.skin_wetness,
        insulation_wetness = evap_pars.insulation_wetness,
        bare_skin_fraction = evap_pars.bare_skin_fraction,
        eye_fraction       = evap_pars.eye_fraction,
    )
    geometry_vars = GeometryVariables(;
        side                    = :dorsal,
        conductance_coefficient = mean_conduction_coeff,
        ventral_fraction,
        conduction_fraction     = ext_cond.conduction_fraction,
        longwave_depth_fraction = ins_pars.longwave_depth_fraction,
    )

    p = (;
        mean_fibre_props,
        mean_ins_pars,
        fat,
        body_shape     = organism.body.shape,
        body_is_sphere = organism.body.shape isa Sphere,
        ventral_fraction,
        ventral_weight,
        dorsal_weight,
        sky_view_factor,
        ground_view_factor,
        ext_cond,
        rad_pars,
        int_cond,
        resp_pars,
        heat_balance_env,
        heat_balance_traits,
        geometry_vars,
        substrate_conductivity = env_vars.substrate_conductivity,
        conduction_depth       = env_pars.conduction_depth,
        insulation_props_init,
    )
    return WeightedMeanNLPPacked(p)
end

# ---------------------------------------------------------------------------
# nlp_pack — MultiSidedNLP
# ---------------------------------------------------------------------------

"""
    nlp_pack(::MultiSidedNLP, organism, environment, T_skin_init, T_ins_init)

Pack per-side physics parameters for the two-sided NLP formulation.
Calls `_pack_sides` at the organism setpoint temperature to get initial temperature
guesses and per-side packed structures. No physics is duplicated.

Returns a `MultiSidedNLPPacked` whose `p` field carries the packed parameters.
"""
function nlp_pack(::MultiSidedNLP, organism::Organism, environment, T_skin_init, T_ins_init)
    env_pars = stripparams(environment.environment_pars)
    env_vars = environment.environment_vars

    ins_pars  = insulation_pars(organism)
    ext_cond  = conduction_pars_external(organism)
    int_cond  = conduction_pars_internal(organism)
    rad_pars  = radiation_pars(organism)
    resp_pars = respiration_pars(organism)
    metab_pars = metabolism_pars(organism)

    fat = Fat(int_cond.fat_fraction, int_cond.fat_density)

    # Call _pack_sides at the setpoint temperature to get initial guesses and
    # per-side packed structures (environment, traits, geometry_vars, bodies)
    packed = _pack_sides(organism, environment, metab_pars.core_temperature, T_skin_init, T_ins_init)

    p = (;
        ins_pars,
        body_shape     = organism.body.shape,
        body_is_sphere = organism.body.shape isa Sphere,
        fat,
        ext_cond,
        int_cond,
        rad_pars,
        resp_pars,
        ventral_fraction       = rad_pars.ventral_fraction,
        substrate_conductivity = env_vars.substrate_conductivity,
        conduction_depth       = env_pars.conduction_depth,
        sky_factor_ref         = packed.sky_factor_ref,
        ground_factor_ref      = packed.ground_factor_ref,
        vegetation_factor_ref  = packed.vegetation_factor_ref,
        area_evaporation       = packed.area_evaporation,
        area_convection        = packed.area_convection,
        fibres                 = packed.fibres,
        evap_pars              = packed.evap_pars,
        side_env_vars_d        = packed.side_env_vars[1],
        side_env_vars_v        = packed.side_env_vars[2],
        side_traits_d          = packed.side_traits[1],
        side_traits_v          = packed.side_traits[2],
        side_geometry_vars_d   = packed.side_geometry_vars[1],
        side_geometry_vars_v   = packed.side_geometry_vars[2],
        # Initial temperature guesses from _pack_sides
        initial_T_skin_d = packed.temps_out[1].skin_temperature,
        initial_T_ins_d  = packed.temps_out[1].insulation_temperature,
        initial_T_skin_v = packed.temps_out[2].skin_temperature,
        initial_T_ins_v  = packed.temps_out[2].insulation_temperature,
    )
    return MultiSidedNLPPacked(p)
end

# ---------------------------------------------------------------------------
# nlp_residuals — WeightedMeanNLP
# ---------------------------------------------------------------------------

"""
    nlp_residuals(p::WeightedMeanNLPPacked, T_core, T_skin, T_ins, Q_gen,
                  k_flesh, pant, skin_wetness, insulation_depth, shape_b)

Evaluate heat-balance constraint residuals for the weighted-mean formulation.
Rebuilds the trial body from `insulation_depth` and `shape_b`, then calls the
existing 4-arg `heat_balance`. No physics is duplicated.

Returns `(; residuals, heat_flows, trial_body, trial_ins_pars)`.
`residuals` is a 3-tuple: `(energy_balance [W], internal_conduction [W], skin_temperature [K])`.
"""
function nlp_residuals(p::WeightedMeanNLPPacked, T_core, T_skin, T_ins, Q_gen,
        k_flesh, pant, skin_wetness, insulation_depth, shape_b)
    pp = p.p

    # Rebuild insulation and body from piloerection and shape_b effectors
    trial_fibre_props = setproperties(pp.mean_fibre_props; depth = insulation_depth)
    trial_ins_pars    = setproperties(pp.mean_ins_pars; dorsal = trial_fibre_props, ventral = trial_fibre_props)
    trial_body        = _nlp_rebuild_side_body(pp.body_shape, pp.body_is_sphere, pp.fat, pp.mean_fibre_props, insulation_depth, shape_b)

    # Recompute temperature-dependent insulation properties
    mean_ins_temperature  = T_ins * 0.7 + T_skin * 0.3
    trial_insulation_props = insulation_properties(trial_ins_pars, mean_ins_temperature, pp.ventral_fraction)

    # Update conductance coefficient for new body size
    trial_conduction_area  = BiophysicalGeometry.total_area(trial_body) * pp.ext_cond.conduction_fraction * 2
    trial_conduction_coeff = trial_conduction_area * pp.substrate_conductivity / pp.conduction_depth * pp.ventral_weight
    trial_geometry_vars    = setproperties(pp.geometry_vars; conductance_coefficient = trial_conduction_coeff)

    balance = heat_balance(
        T_core, T_skin, T_ins, Q_gen;
        body             = trial_body,
        insulation_pars  = trial_ins_pars,
        insulation       = trial_insulation_props,
        geometry_vars    = trial_geometry_vars,
        environment_vars = pp.heat_balance_env,
        traits           = pp.heat_balance_traits,
        resp_pars        = pp.resp_pars,
        k_flesh,
        pant,
        skin_wetness,
    )

    return (;
        residuals  = (balance.residual_energy_balance, balance.residual_internal_conduction,
                      balance.residual_skin_temperature),
        heat_flows = balance,
        trial_body,
        trial_ins_pars,
    )
end

# ---------------------------------------------------------------------------
# nlp_residuals — MultiSidedNLP
# ---------------------------------------------------------------------------

"""
    nlp_residuals(p::MultiSidedNLPPacked, T_core, T_skin_d, T_ins_d, T_skin_v, T_ins_v,
                  Q_gen, k_flesh, pant, skin_wetness, insulation_depth, shape_b)

Evaluate heat-balance constraint residuals for the two-sided NLP formulation.
Calls the 4-arg `heat_balance` once for each body side.

Returns `(; residuals, heat_flows_dorsal, heat_flows_ventral, trial_body_d, trial_body_v, trial_ins_pars)`.
`residuals` is a 6-tuple: dorsal (energy_balance, internal_conduction, skin_temperature),
ventral (energy_balance, internal_conduction, skin_temperature).
"""
function nlp_residuals(p::MultiSidedNLPPacked, T_core, T_skin_d, T_ins_d, T_skin_v, T_ins_v,
        Q_gen, k_flesh, pant, skin_wetness, insulation_depth, shape_b)
    pp = p.p

    # Update insulation depth for both sides (shared depth effector)
    trial_dorsal_fibre  = setproperties(pp.ins_pars.dorsal;  depth = insulation_depth)
    trial_ventral_fibre = setproperties(pp.ins_pars.ventral; depth = insulation_depth)
    trial_ins_pars      = setproperties(pp.ins_pars; dorsal = trial_dorsal_fibre, ventral = trial_ventral_fibre)

    # Per-side bodies (use side-specific fibre diameter and density)
    trial_body_d = _nlp_rebuild_side_body(pp.body_shape, pp.body_is_sphere, pp.fat, pp.ins_pars.dorsal,  insulation_depth, shape_b)
    trial_body_v = _nlp_rebuild_side_body(pp.body_shape, pp.body_is_sphere, pp.fat, pp.ins_pars.ventral, insulation_depth, shape_b)

    # Per-side insulation properties at current temperatures
    trial_ins_props_d = insulation_properties(trial_ins_pars, T_ins_d * 0.7 + T_skin_d * 0.3, pp.ventral_fraction)
    trial_ins_props_v = insulation_properties(trial_ins_pars, T_ins_v * 0.7 + T_skin_v * 0.3, pp.ventral_fraction)

    # Update ventral conductance coefficient for new body size (dorsal is always 0)
    trial_cond_area_v   = BiophysicalGeometry.total_area(trial_body_v) * pp.ext_cond.conduction_fraction * 2
    trial_cond_coeff_v  = trial_cond_area_v * pp.substrate_conductivity / pp.conduction_depth
    trial_geom_vars_v   = setproperties(pp.side_geometry_vars_v; conductance_coefficient = trial_cond_coeff_v)

    balance_d = heat_balance(
        T_core, T_skin_d, T_ins_d, Q_gen;
        body             = trial_body_d,
        insulation_pars  = trial_ins_pars,
        insulation       = trial_ins_props_d,
        geometry_vars    = pp.side_geometry_vars_d,
        environment_vars = pp.side_env_vars_d,
        traits           = pp.side_traits_d,
        resp_pars        = pp.resp_pars,
        k_flesh, pant, skin_wetness,
    )
    balance_v = heat_balance(
        T_core, T_skin_v, T_ins_v, Q_gen;
        body             = trial_body_v,
        insulation_pars  = trial_ins_pars,
        insulation       = trial_ins_props_v,
        geometry_vars    = trial_geom_vars_v,
        environment_vars = pp.side_env_vars_v,
        traits           = pp.side_traits_v,
        resp_pars        = pp.resp_pars,
        k_flesh, pant, skin_wetness,
    )

    return (;
        residuals = (
            balance_d.residual_energy_balance, balance_d.residual_internal_conduction, balance_d.residual_skin_temperature,
            balance_v.residual_energy_balance, balance_v.residual_internal_conduction, balance_v.residual_skin_temperature,
        ),
        heat_flows_dorsal  = balance_d,
        heat_flows_ventral = balance_v,
        trial_body_d,
        trial_body_v,
        trial_ins_pars,
    )
end

# ---------------------------------------------------------------------------
# nlp_assemble_output — WeightedMeanNLP
# ---------------------------------------------------------------------------

"""
    nlp_assemble_output(p::WeightedMeanNLPPacked, organism, environment,
                        T_core, T_skin, T_ins, Q_gen,
                        k_flesh, pant, skin_wetness, insulation_depth, shape_b)

Assemble `(; thermoregulation, morphology, energy_flows, mass_flows)` from a converged
WeightedMeanNLP solution. Calls `nlp_residuals` once at the solution for final heat flows,
then delegates mass flows to `respiration` and geometry to BiophysicalGeometry functions.
Output field structure matches `solve_metabolic_rate`.
"""
function nlp_assemble_output(p::WeightedMeanNLPPacked, organism::Organism, environment,
        T_core, T_skin, T_ins, Q_gen, k_flesh, pant, skin_wetness, insulation_depth, shape_b)
    pp = p.p
    air_temperature = environment.environment_vars.air_temperature

    r = nlp_residuals(p, T_core, T_skin, T_ins, Q_gen, k_flesh, pant, skin_wetness, insulation_depth, shape_b)
    hf          = r.heat_flows
    sol_body    = r.trial_body
    sol_ins_pars = r.trial_ins_pars

    lung_temperature = (T_core + T_skin) / 2

    # Respiration mass flows at solution
    panting_resp_pars  = setproperties(pp.resp_pars; pant)
    respiration_result = respiration(
        MetabolicRates(; metabolic = Q_gen, sum = Q_gen, minimum = zero(Q_gen)),
        panting_resp_pars,
        pp.heat_balance_env.atmos,
        sol_body.shape.mass,
        lung_temperature,
        air_temperature;
        gas_fractions = pp.heat_balance_env.gas_fractions,
        O2conversion  = Kleiber1961(),
    )
    latent_heat_vap = enthalpy_of_vaporisation(air_temperature)
    m_sweat = u"g/hr"(hf.skin_evaporation_heat_flow / latent_heat_vap)
    m_evap  = u"g/hr"(respiration_result.respiration_mass + m_sweat)

    # Longwave flows (Stefan-Boltzmann at insulation surface)
    σ                    = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    sol_total_area       = BiophysicalGeometry.total_area(sol_body)
    sol_ventral_rad_area = sol_total_area * (1 - pp.ext_cond.conduction_fraction)
    dorsal_lw_out        = _side_lw_out(pp.sky_view_factor,    pp.rad_pars.body_emissivity_dorsal,  sol_total_area,       T_ins, σ)
    ventral_lw_out       = _side_lw_out(pp.ground_view_factor, pp.rad_pars.body_emissivity_ventral, sol_ventral_rad_area, T_ins, σ)
    longwave_flow_out    = dorsal_lw_out * pp.dorsal_weight + ventral_lw_out * pp.ventral_weight
    longwave_flow_in     = longwave_flow_out - hf.radiation_heat_flow

    # Insulation properties at solution temperatures
    sol_ins_props = insulation_properties(sol_ins_pars, T_ins * 0.7 + T_skin * 0.3, pp.ventral_fraction)
    sol_shape_b   = pp.body_is_sphere ? 1.0 : sol_body.shape.aspect_ratio_b

    thermoregulation = (;
        core_temperature               = T_core,
        skin_temperature               = T_skin,
        insulation_temperature         = T_ins,
        lung_temperature,
        skin_temperature_dorsal        = T_skin,
        skin_temperature_ventral       = T_skin,
        insulation_temperature_dorsal  = T_ins,
        insulation_temperature_ventral = T_ins,
        shape_b                        = sol_shape_b,
        pant,
        skin_wetness,
        flesh_conductivity             = k_flesh,
        insulation_conductivity_effective  = sol_ins_props.conductivities.average,
        insulation_conductivity_dorsal     = sol_ins_props.conductivities.dorsal,
        insulation_conductivity_ventral    = sol_ins_props.conductivities.ventral,
        insulation_conductivity_compressed = sol_ins_props.conductivity_compressed,
        insulation_depth_dorsal  = sol_ins_pars.dorsal.depth,
        insulation_depth_ventral = sol_ins_pars.ventral.depth,
        q10 = metabolism_pars(organism).q10,
        dorsal = (;
            skin_temperature        = T_skin,
            insulation_temperature  = T_ins,
            insulation_conductivity = sol_ins_props.conductivities.dorsal,
        ),
        ventral = (;
            skin_temperature        = T_skin,
            insulation_temperature  = T_ins,
            insulation_conductivity = sol_ins_props.conductivities.ventral,
        ),
    )

    area_skin        = skin_area(sol_body)
    area_evaporation = evaporation_area(sol_body)
    area_convection  = sol_total_area * (1 - pp.ext_cond.conduction_fraction)
    area_silhouette  = silhouette_area(sol_body, pp.rad_pars.solar_orientation)

    morphology = (;
        total_area         = sol_total_area,
        area_skin,
        area_evaporation,
        area_convection,
        area_conduction    = sol_total_area * pp.ext_cond.conduction_fraction / 2,
        area_silhouette,
        sky_view_factor    = pp.sky_view_factor,
        ground_view_factor = pp.ground_view_factor,
        volume             = sol_body.geometry.volume,
        volume_flesh       = flesh_volume(sol_body),
        characteristic_dimension = characteristic_dimension(VolumeCubeRoot(), sol_body),
        fat_mass           = sol_body.shape.mass * pp.fat.fraction,
        sol_body.geometry.length...,
    )

    evaporation_heat_flow = hf.skin_evaporation_heat_flow + hf.insulation_evaporation_heat_flow +
                            hf.respiration_heat_flow

    energy_flows = (;
        solar_flow            = hf.solar_heat_flow,
        longwave_flow_in,
        longwave_flow_out,
        generated_heat_flow   = Q_gen,
        respiration_heat_flow = hf.respiration_heat_flow,
        evaporation_heat_flow,
        convection_heat_flow  = hf.convection_heat_flow,
        conduction_flow       = hf.conduction_heat_flow,
        heat_balance          = hf.residual_energy_balance,
        balance               = hf.residual_energy_balance,
        ntry    = 1,
        success = true,
    )

    mass_flows = (;
        air_flow             = respiration_result.air_flow,
        oxygen_flow_standard = respiration_result.oxygen_flow_standard,
        m_evap,
        respiration_mass     = respiration_result.respiration_mass,
        m_sweat,
        molar_fluxes_in  = respiration_result.molar_fluxes_in,
        molar_fluxes_out = respiration_result.molar_fluxes_out,
    )

    return (; thermoregulation, morphology, energy_flows, mass_flows)
end

# ---------------------------------------------------------------------------
# nlp_assemble_output — MultiSidedNLP
# ---------------------------------------------------------------------------

"""
    nlp_assemble_output(p::MultiSidedNLPPacked, organism, environment,
                        T_core, T_skin_d, T_ins_d, T_skin_v, T_ins_v,
                        Q_gen, k_flesh, pant, skin_wetness, insulation_depth, shape_b,
                        respiration_out)

Assemble `(; thermoregulation, morphology, energy_flows, mass_flows)` from a converged
MultiSidedNLP solution. `respiration_out` is the result of a `respiration(...)` call
computed by the caller (BiophysicalBehaviour.jl), or `nothing` when respiration is disabled.
Output field structure matches `solve_metabolic_rate`, including per-side sub-fields.
"""
function nlp_assemble_output(p::MultiSidedNLPPacked, organism::Organism, environment,
        T_core, T_skin_d, T_ins_d, T_skin_v, T_ins_v,
        Q_gen, k_flesh, pant, skin_wetness, insulation_depth, shape_b,
        respiration_out)
    pp = p.p
    env_vars = environment.environment_vars

    r    = nlp_residuals(p, T_core, T_skin_d, T_ins_d, T_skin_v, T_ins_v,
                          Q_gen, k_flesh, pant, skin_wetness, insulation_depth, shape_b)
    hf_d = r.heat_flows_dorsal
    hf_v = r.heat_flows_ventral
    sol_body_d   = r.trial_body_d
    sol_body_v   = r.trial_body_v
    sol_ins_pars = r.trial_ins_pars

    dmult = pp.sky_factor_ref + pp.vegetation_factor_ref
    vmult = 1 - dmult

    # Insulation properties at solution temperatures
    sol_ins_props_d = insulation_properties(sol_ins_pars, T_ins_d * 0.7 + T_skin_d * 0.3, pp.ventral_fraction)
    sol_ins_props_v = insulation_properties(sol_ins_pars, T_ins_v * 0.7 + T_skin_v * 0.3, pp.ventral_fraction)

    # Whole-organism weighted means
    skin_temperature       = T_skin_d * dmult + T_skin_v * vmult
    insulation_temperature = T_ins_d  * dmult + T_ins_v  * vmult
    lung_temperature       = (T_core + skin_temperature) / 2

    # Longwave flows (Stefan-Boltzmann at per-side insulation surface temperatures)
    σ                    = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    sol_total_area       = BiophysicalGeometry.total_area(sol_body_d)
    sol_ventral_rad_area = sol_total_area * (1 - pp.ext_cond.conduction_fraction)
    dorsal_lw_out        = _side_lw_out(pp.sky_factor_ref,    pp.rad_pars.body_emissivity_dorsal,  sol_total_area,       T_ins_d, σ)
    ventral_lw_out       = _side_lw_out(pp.ground_factor_ref, pp.rad_pars.body_emissivity_ventral, sol_ventral_rad_area, T_ins_v, σ)
    dorsal_lw_in         = dorsal_lw_out  - hf_d.radiation_heat_flow
    ventral_lw_in        = ventral_lw_out - hf_v.radiation_heat_flow
    longwave_flow_out    = dorsal_lw_out * dmult + ventral_lw_out * vmult
    longwave_flow_in     = dorsal_lw_in  * dmult + ventral_lw_in  * vmult

    # Respiration
    (; respiration_heat_flow, respiration_mass, air_flow, oxygen_flow_standard,
       molar_fluxes_in, molar_fluxes_out) = _unpack_respiration(respiration_out)

    latent_heat_vap = enthalpy_of_vaporisation(env_vars.air_temperature)
    m_sweat = u"g/hr"((hf_d.skin_evaporation_heat_flow + hf_v.skin_evaporation_heat_flow) * 0.5 /
                       latent_heat_vap)
    m_evap  = !isnothing(respiration_mass) ? u"g/hr"(respiration_mass + m_sweat) : u"g/hr"(m_sweat)

    evaporation_heat_flow =
        hf_d.skin_evaporation_heat_flow       * dmult + hf_v.skin_evaporation_heat_flow       * vmult +
        hf_d.insulation_evaporation_heat_flow * dmult + hf_v.insulation_evaporation_heat_flow * vmult +
        respiration_heat_flow

    # Whole-organism insulation properties
    insulation_final = insulation_properties(sol_ins_pars,
        insulation_temperature * 0.7 + skin_temperature * 0.3, pp.ventral_fraction)
    sol_shape_b = pp.body_is_sphere ? 1.0 : sol_body_d.shape.aspect_ratio_b

    thermoregulation = (;
        core_temperature               = T_core,
        skin_temperature,
        insulation_temperature,
        lung_temperature,
        skin_temperature_dorsal        = T_skin_d,
        skin_temperature_ventral       = T_skin_v,
        insulation_temperature_dorsal  = T_ins_d,
        insulation_temperature_ventral = T_ins_v,
        shape_b                        = sol_shape_b,
        pant,
        skin_wetness,
        flesh_conductivity             = k_flesh,
        insulation_conductivity_effective  = insulation_final.conductivities.average,
        insulation_conductivity_dorsal     = sol_ins_props_d.conductivities.average,
        insulation_conductivity_ventral    = sol_ins_props_v.conductivities.average,
        insulation_conductivity_compressed = insulation_final.conductivity_compressed,
        insulation_depth_dorsal  = sol_ins_pars.dorsal.depth,
        insulation_depth_ventral = sol_ins_pars.ventral.depth,
        q10 = metabolism_pars(organism).q10,
        dorsal = (;
            skin_temperature        = T_skin_d,
            insulation_temperature  = T_ins_d,
            insulation_conductivity = sol_ins_props_d.conductivities.average,
        ),
        ventral = (;
            skin_temperature        = T_skin_v,
            insulation_temperature  = T_ins_v,
            insulation_conductivity = sol_ins_props_v.conductivities.average,
        ),
    )

    area_evaporation = evaporation_area(sol_body_d)
    area_convection  = sol_total_area * (1 - pp.ext_cond.conduction_fraction)
    area_skin        = skin_area(sol_body_d)
    area_silhouette  = silhouette_area(sol_body_d, pp.rad_pars.solar_orientation)

    morphology = (;
        total_area         = sol_total_area,
        area_skin,
        area_evaporation,
        area_convection,
        area_conduction    = BiophysicalGeometry.total_area(sol_body_v) * pp.ext_cond.conduction_fraction,
        area_silhouette,
        sky_view_factor    = pp.sky_factor_ref,
        ground_view_factor = pp.ground_factor_ref,
        volume             = sol_body_d.geometry.volume,
        volume_flesh       = flesh_volume(sol_body_d),
        characteristic_dimension = characteristic_dimension(VolumeCubeRoot(), sol_body_d),
        fat_mass           = sol_body_d.shape.mass * pp.fat.fraction,
        sol_body_d.geometry.length...,
    )

    solar_flow           = hf_d.solar_heat_flow * dmult + hf_v.solar_heat_flow * vmult
    convection_heat_flow = hf_d.convection_heat_flow * dmult + hf_v.convection_heat_flow * vmult
    conduction_flow      = hf_d.conduction_heat_flow * dmult + hf_v.conduction_heat_flow * vmult
    heat_balance_val     = solar_flow + longwave_flow_in + Q_gen -
                           longwave_flow_out - convection_heat_flow - evaporation_heat_flow - conduction_flow

    energy_flows = (;
        solar_flow,
        longwave_flow_in,
        longwave_flow_out,
        generated_heat_flow   = Q_gen,
        respiration_heat_flow,
        evaporation_heat_flow,
        convection_heat_flow,
        conduction_flow,
        heat_balance          = heat_balance_val,
        balance               = nothing,
        ntry    = 1,
        success = true,
        dorsal = (;
            solar                  = hf_d.solar_heat_flow,
            convection             = hf_d.convection_heat_flow,
            conduction             = hf_d.conduction_heat_flow,
            net_generated          = hf_d.net_metabolic_heat_internal,
            skin_evaporation       = hf_d.skin_evaporation_heat_flow,
            insulation_evaporation = hf_d.insulation_evaporation_heat_flow,
            longwave               = hf_d.radiation_heat_flow,
            sky_radiation          = hf_d.sky_radiation_flow,
            bush_radiation         = hf_d.bush_radiation_flow,
            vegetation_radiation   = hf_d.vegetation_radiation_flow,
            ground_radiation       = hf_d.ground_radiation_flow,
        ),
        ventral = (;
            solar                  = hf_v.solar_heat_flow,
            convection             = hf_v.convection_heat_flow,
            conduction             = hf_v.conduction_heat_flow,
            net_generated          = hf_v.net_metabolic_heat_internal,
            skin_evaporation       = hf_v.skin_evaporation_heat_flow,
            insulation_evaporation = hf_v.insulation_evaporation_heat_flow,
            longwave               = hf_v.radiation_heat_flow,
            sky_radiation          = hf_v.sky_radiation_flow,
            bush_radiation         = hf_v.bush_radiation_flow,
            vegetation_radiation   = hf_v.vegetation_radiation_flow,
            ground_radiation       = hf_v.ground_radiation_flow,
        ),
    )

    mass_flows = (;
        air_flow,
        oxygen_flow_standard,
        m_evap,
        respiration_mass,
        m_sweat,
        molar_fluxes_in,
        molar_fluxes_out,
    )

    return (; thermoregulation, morphology, energy_flows, mass_flows)
end
