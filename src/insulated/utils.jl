# Unpack a respiration result (or nothing) into a flat NamedTuple of scalar fields.
# Used by _assemble_multisided_output and nlp_assemble_output(::MultiSidedNLPPacked).
function _unpack_respiration(resp_out)
    isnothing(resp_out) && return (;
        balance               = nothing,
        respiration_heat_flow = 0.0u"W",
        respiration_mass_flow      = nothing,
        air_flow              = nothing,
        oxygen_flow_standard  = nothing,
        molar_fluxes_in       = nothing,
        molar_fluxes_out      = nothing,
    )
    return (;
        balance               = resp_out.balance,
        respiration_heat_flow = resp_out.respiration_heat_flow,
        respiration_mass_flow      = resp_out.respiration_mass_flow,
        air_flow              = resp_out.air_flow,
        oxygen_flow_standard  = resp_out.oxygen_flow_standard,
        molar_fluxes_in       = resp_out.molar_fluxes_in,
        molar_fluxes_out      = resp_out.molar_fluxes_out,
    )
end

# Stefan-Boltzmann outgoing longwave for one body side.
# The factor 2 converts a half-hemisphere view factor to a whole-body effective factor.
_side_lw_out(view_factor, emissivity, area, surface_temperature, σ) =
    2 * view_factor * emissivity * area * σ * surface_temperature^4

"""
    interpolate(a::FibreProperties, b::FibreProperties, t::Number)
    interpolate(a::Number, b::Number, t)

Linearly interpolate between two `FibreProperties` objects.

Returns `a * (1 - t) + b * t` for each field.
"""
function interpolate(a::FibreProperties, b::FibreProperties, t::Number)
    interpolated = map(getproperties(a), getproperties(b)) do a, b
        interpolate(a, b, t)
    end
    return setproperties(a; interpolated...)
end
interpolate(a::Number, b::Number, t::Number) = a * (1 - t) + b * t

"""
    _pack_sides(o::Organism, e, core_temperature, skin_temperature, insulation_temperature)

Set up per-side (dorsal/ventral) packed structures for insulated heat balance and run
`solve_temperatures` on each side. Encapsulates the pre-loop setup and per-side loop
shared by `solve_metabolic_rate` and `heat_balance(T, ::MultiSided, ...)`.

Returns a NamedTuple with:
- `temps_out` — 2-element vector of `solve_temperatures` outputs (dorsal, ventral)
- `side_bodies`, `side_geometry_vars`, `side_env_vars`, `side_traits` — per-side structures
- `sky_factor_ref`, `ground_factor_ref`, `bush_factor_ref`, `vegetation_factor_ref`
- `skin_temperature`, `insulation_temperature`, `temperature_error_tolerance` (converged)
- `area_evaporation`, `area_convection`, `total_area`
- `fibres` — insulation fibre properties (for geometry_avg construction)
- `evap_pars` — possibly adjusted evaporation parameters (bare_skin_fraction reset if no insulation)
"""
function _pack_sides(o::Organism, e, core_temperature, skin_temperature, insulation_temperature)
    environment_pars = stripparams(e.environment_pars)
    environment_vars = e.environment_vars
    opts = options(o)
    ins_pars = insulation_pars(o)
    external_conduction = conduction_pars_external(o)
    internal_conduction = conduction_pars_internal(o)
    rad_pars = radiation_pars(o)
    evap_pars = evaporation_pars(o)

    avg_insulation_temp = insulation_temperature * 0.7 + skin_temperature * 0.3
    insulation = insulation_properties(ins_pars, avg_insulation_temp, rad_pars.ventral_fraction)
    fibres = insulation.fibres
    if insulation.insulation_test <= 0.0u"m" && evap_pars.bare_skin_fraction < 1.0
        evap_pars = AnimalEvaporationParameters(;
            skin_wetness=evap_pars.skin_wetness,
            insulation_wetness=evap_pars.insulation_wetness,
            eye_fraction=evap_pars.eye_fraction,
            bare_skin_fraction=1.0,
            insulation_fraction=evap_pars.insulation_fraction,
        )
    end

    fat = FatLayer(internal_conduction.fat_fraction, internal_conduction.fat_density)

    vegetation_factor = rad_pars.sky_view_factor * environment_vars.shade
    sky_factor = rad_pars.sky_view_factor - vegetation_factor
    ground_factor = 1 - sky_factor - vegetation_factor

    area_silhouette = silhouette_area(o.body, rad_pars.solar_orientation)
    total_area = BiophysicalGeometry.total_area(o.body)
    area_conduction = total_area * external_conduction.conduction_fraction
    area_evaporation = evaporation_area(o.body)
    area_convection = total_area * (1 - external_conduction.conduction_fraction)
    absorptivities = Absorptivities(rad_pars, environment_pars)
    view_factors_solar = ViewFactors(sky_factor, ground_factor, 0.0, 0.0)
    solar_conditions = SolarConditions(environment_vars)
    (; solar_flow, solar_direct_flow, solar_sky_flow, solar_substrate_flow) = solar(
        o.body, absorptivities, view_factors_solar, solar_conditions, area_silhouette, area_conduction,
    )
    ventral_flow = solar_substrate_flow

    vegetation_temperature = environment_vars.reference_air_temperature
    bush_factor_ref = rad_pars.bush_view_factor
    sky_factor_ref = sky_factor
    ground_factor_ref = ground_factor
    vegetation_factor_ref = vegetation_factor

    temps_out        = Vector{NamedTuple}(undef, 2)
    side_bodies        = Vector{Any}(undef, 2)
    side_geometry_vars = Vector{GeometryVariables}(undef, 2)
    side_env_vars      = Vector{NamedTuple}(undef, 2)
    side_traits        = Vector{NamedTuple}(undef, 2)
    temperature_error_tolerance = opts.temperature_error_tolerance

    for (side_idx, side) in enumerate((:dorsal, :ventral))
        side_sky_factor, side_ground_factor, side_bush_factor, side_vegetation_factor = if solar_flow > 0.0u"W"
            if side == :dorsal
                solar_flow = 2.0 * solar_direct_flow + ((solar_sky_flow / sky_factor_ref) * sky_factor_ref * 2.0)
                (sky_factor_ref * 2.0, 0.0, 0.0, vegetation_factor_ref * 2.0)
            else
                solar_flow =
                    (ventral_flow / (1.0 - sky_factor_ref - vegetation_factor_ref)) *
                    (1.0 - (2.0 * external_conduction.conduction_fraction))
                (0.0, ground_factor_ref * 2.0, bush_factor_ref * 2.0, 0.0)
            end
        else
            solar_flow = 0.0u"W"
            if side == :dorsal
                (sky_factor_ref * 2.0, 0.0, 0.0, vegetation_factor_ref * 2.0)
            else
                (0.0, ground_factor_ref * 2.0, bush_factor_ref * 2.0, 0.0)
            end
        end

        ϵ_body = side == :dorsal ? rad_pars.body_emissivity_dorsal : rad_pars.body_emissivity_ventral

        side_fibres = getproperty(fibres, side)
        insulation_layer = FibrousLayer(side_fibres.depth, side_fibres.diameter, side_fibres.density)
        geometry_pars = Body(o.body.shape, CompositeInsulation(insulation_layer, fat))

        conductance_coefficient = if side == :ventral
            body_total_area = BiophysicalGeometry.total_area(geometry_pars)
            area_cond = body_total_area * external_conduction.conduction_fraction * 2
            (area_cond * environment_vars.substrate_conductivity) / environment_pars.conduction_depth
        else
            0.0u"W/K"
        end

        geometry_vars = GeometryVariables(;
            side,
            conductance_coefficient,
            ventral_fraction=rad_pars.ventral_fraction,
            conduction_fraction=external_conduction.conduction_fraction,
            longwave_depth_fraction=ins_pars.longwave_depth_fraction,
        )
        view_factors = ViewFactors(side_sky_factor, side_ground_factor, side_bush_factor, side_vegetation_factor)
        atmos = AtmosphericConditions(environment_vars)
        temperature = EnvironmentTemperatures(
            environment_vars.air_temperature, environment_vars.sky_temperature,
            environment_vars.ground_temperature, vegetation_temperature,
            environment_vars.bush_temperature, environment_vars.substrate_temperature,
        )
        env_packed = (;
            temperature,
            view_factors,
            atmos,
            fluid=environment_pars.fluid,
            solar_flow,
            gas_fractions=environment_pars.gas_fractions,
            convection_enhancement=environment_pars.convection_enhancement,
        )
        traits = (;
            core_temperature,
            flesh_conductivity=internal_conduction.flesh_conductivity,
            fat_conductivity=internal_conduction.fat_conductivity,
            ϵ_body,
            skin_wetness=evap_pars.skin_wetness,
            insulation_wetness=evap_pars.insulation_wetness,
            bare_skin_fraction=evap_pars.bare_skin_fraction,
            eye_fraction=evap_pars.eye_fraction,
        )

        insulation = insulation_properties(
            ins_pars, insulation_temperature * 0.7 + skin_temperature * 0.3, rad_pars.ventral_fraction
        )
        temps_out[side_idx] = solve_temperatures(;
            body=geometry_pars,
            insulation_pars=ins_pars,
            insulation,
            geometry_vars,
            environment_vars=env_packed,
            traits,
            temperature_tolerance=opts.temperature_error_tolerance,
            skin_temperature,
            insulation_temperature,
        )

        side_bodies[side_idx]        = geometry_pars
        side_geometry_vars[side_idx] = geometry_vars
        side_env_vars[side_idx]      = env_packed
        side_traits[side_idx]        = traits

        skin_temperature = temps_out[side_idx].skin_temperature
        insulation_temperature = temps_out[side_idx].insulation_temperature
        temperature_error_tolerance = temps_out[side_idx].tolerance
    end

    return (;
        temps_out,
        side_bodies, side_geometry_vars, side_env_vars, side_traits,
        sky_factor_ref, ground_factor_ref, bush_factor_ref, vegetation_factor_ref,
        skin_temperature, insulation_temperature, temperature_error_tolerance,
        area_evaporation, area_convection, total_area,
        fibres, evap_pars,
    )
end

"""
    _assemble_multisided_output(o, e, core_temperature, metabolic_heat_flow, respiration_out, packed)

Assemble the rich `(; thermoregulation, morphology, energy_flows, mass_flows)` output tuple
from `_pack_sides` results. Shared by `solve_metabolic_rate` and `solve_temperature(::MultiSided)`.

`metabolic_heat_flow` is the metabolic heat rate (solved by zbrent for endotherms, or from the
metabolic model for ectotherms). `respiration_out` is the result of a `respiration(...)` call,
or `nothing` when respiration is disabled.

`energy_flows` includes `metabolic_heat_flow` and per-side `dorsal`/`ventral` sub-fields.
`thermoregulation` includes per-side `dorsal`/`ventral` sub-fields alongside existing flat scalars.
"""
function _assemble_multisided_output(o, e, core_temperature, metabolic_heat_flow, respiration_out, packed)
    environment_vars = e.environment_vars
    ins_pars = insulation_pars(o)
    external_conduction = conduction_pars_external(o)
    internal_conduction = conduction_pars_internal(o)
    rad_pars = radiation_pars(o)
    resp_pars = respiration_pars(o)
    metab_pars = metabolism_pars(o)

    (;
        temps_out, side_bodies,
        sky_factor_ref, ground_factor_ref, vegetation_factor_ref,
        area_evaporation, area_convection, fibres, evap_pars,
    ) = packed

    dmult = sky_factor_ref + vegetation_factor_ref
    vmult = 1 - dmult

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    # Simple average for insulation test (determines whether fur is present)
    skin_temp_avg = (temps_out[1].skin_temperature + temps_out[2].skin_temperature) * 0.5
    ins_temp_avg  = (temps_out[1].insulation_temperature + temps_out[2].insulation_temperature) * 0.5
    insulation_for_test = insulation_properties(
        ins_pars, ins_temp_avg * 0.7 + skin_temp_avg * 0.3, rad_pars.ventral_fraction
    )
    has_insulation = insulation_for_test.insulation_test > 0.0u"m"

    fat = FatLayer(internal_conduction.fat_fraction, internal_conduction.fat_density)
    geometry_avg = Body(o.body.shape, CompositeInsulation(
        FibrousLayer(fibres.average.depth, fibres.average.diameter, fibres.average.density), fat
    ))
    total_area_final = BiophysicalGeometry.total_area(geometry_avg)
    area_skin_final  = skin_area(geometry_avg)

    if has_insulation
        surface_temperature_dorsal  = temps_out[1].insulation_temperature
        surface_temperature_ventral = temps_out[2].insulation_temperature
        area_rad_d = total_area_final
        area_rad_v = total_area_final * (1 - external_conduction.conduction_fraction)
    else
        surface_temperature_dorsal  = temps_out[1].skin_temperature
        surface_temperature_ventral = temps_out[2].skin_temperature
        area_rad_d = area_skin_final
        area_rad_v = area_skin_final
    end

    dorsal_radiation_out_flow  = _side_lw_out(sky_factor_ref,    rad_pars.body_emissivity_dorsal,  area_rad_d, surface_temperature_dorsal, σ)
    ventral_radiation_out_flow = _side_lw_out(ground_factor_ref, rad_pars.body_emissivity_ventral, area_rad_v, surface_temperature_ventral, σ)
    dorsal_radiation_in_flow   = dorsal_radiation_out_flow  - temps_out[1].flows.longwave
    ventral_radiation_in_flow  = ventral_radiation_out_flow - temps_out[2].flows.longwave

    solar_flow           = temps_out[1].flows.solar * dmult + temps_out[2].flows.solar * vmult
    longwave_flow_in     = dorsal_radiation_in_flow  * dmult + ventral_radiation_in_flow  * vmult
    longwave_flow_out    = dorsal_radiation_out_flow * dmult + ventral_radiation_out_flow * vmult
    convection_heat_flow = temps_out[1].flows.convection * dmult + temps_out[2].flows.convection * vmult
    conduction_flow      = temps_out[1].flows.conduction * dmult + temps_out[2].flows.conduction * vmult
    evaporation_heat_flow =
        temps_out[1].flows.skin_evaporation * dmult + temps_out[2].flows.skin_evaporation * vmult +
        temps_out[1].flows.insulation_evaporation * dmult + temps_out[2].flows.insulation_evaporation * vmult

    (; balance, respiration_heat_flow, respiration_mass_flow, air_flow, oxygen_flow_standard,
       molar_fluxes_in, molar_fluxes_out) = _unpack_respiration(respiration_out)

    evaporation_heat_flow += respiration_heat_flow

    latent_heat_vaporisation = enthalpy_of_vaporisation(environment_vars.air_temperature)
    m_sweat = u"g/hr"(
        (temps_out[1].flows.skin_evaporation + temps_out[2].flows.skin_evaporation) * 0.5 /
        latent_heat_vaporisation
    )
    m_evap = !isnothing(respiration_mass_flow) ? u"g/hr"(respiration_mass_flow + m_sweat) : u"g/hr"(m_sweat)

    fat_mass     = geometry_avg.shape.mass * fat.fraction
    volume       = geometry_avg.geometry.volume
    volume_flesh = flesh_volume(geometry_avg)

    skin_temperature       = temps_out[1].skin_temperature * dmult + temps_out[2].skin_temperature * vmult
    insulation_temperature = temps_out[1].insulation_temperature * dmult + temps_out[2].insulation_temperature * vmult
    insulation_depth       = ins_pars.dorsal.depth * dmult + ins_pars.ventral.depth * vmult
    insulation_final = insulation_properties(
        ins_pars, insulation_temperature * 0.7 + skin_temperature * 0.3, rad_pars.ventral_fraction
    )
    insulation_conductivity_effective  = insulation_final.conductivities.average
    insulation_conductivity_compressed = insulation_final.conductivity_compressed

    lung_temperature = (core_temperature + skin_temperature) * 0.5
    axis_ratio_b = o.body.shape isa Sphere ? 1.0 : o.body.shape.axis_ratio_b

    heat_balance_val = solar_flow + longwave_flow_in + metabolic_heat_flow -
        longwave_flow_out - convection_heat_flow - evaporation_heat_flow - conduction_flow

    thermoregulation = (;
        core_temperature,
        skin_temperature,
        insulation_temperature,
        lung_temperature,
        insulation_depth,
        axis_ratio_b,
        pant               = resp_pars.pant,
        skin_wetness       = evap_pars.skin_wetness,
        flesh_conductivity = internal_conduction.flesh_conductivity,
        insulation_conductivity_effective,
        insulation_conductivity_compressed,
        q10 = metab_pars.q10,
        dorsal = (;
            skin_temperature        = temps_out[1].skin_temperature,
            insulation_temperature  = temps_out[1].insulation_temperature,
            insulation_conductivity = temps_out[1].insulation_conductivity,
            insulation_depth        = ins_pars.dorsal.depth,
        ),
        ventral = (;
            skin_temperature        = temps_out[2].skin_temperature,
            insulation_temperature  = temps_out[2].insulation_temperature,
            insulation_conductivity = temps_out[2].insulation_conductivity,
            insulation_depth        = ins_pars.ventral.depth,
        ),
    )

    morphology = (;
        total_area         = total_area_final,
        area_skin          = area_skin_final,
        area_evaporation,
        area_convection,
        area_conduction    = BiophysicalGeometry.total_area(side_bodies[2]) * external_conduction.conduction_fraction,
        area_silhouette    = silhouette_area(geometry_avg, rad_pars.solar_orientation),
        sky_view_factor    = sky_factor_ref,
        ground_view_factor = ground_factor_ref,
        volume,
        volume_flesh,
        characteristic_dimension = characteristic_dimension(VolumeCubeRoot(), geometry_avg),
        fat_mass,
        geometry_avg.geometry.length...,
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
        ntry    = max(temps_out[1].ntry, temps_out[2].ntry),
        success = all((temps_out[1].success, temps_out[2].success)),
        dorsal = (;
            solar                  = temps_out[1].flows.solar,
            convection             = temps_out[1].flows.convection,
            conduction             = temps_out[1].flows.conduction,
            net_metabolic          = temps_out[1].flows.net_metabolic,
            skin_evaporation       = temps_out[1].flows.skin_evaporation,
            insulation_evaporation = temps_out[1].flows.insulation_evaporation,
            longwave               = temps_out[1].flows.longwave,
            sky_radiation          = temps_out[1].flows.sky_radiation,
            bush_radiation         = temps_out[1].flows.bush_radiation,
            vegetation_radiation   = temps_out[1].flows.vegetation_radiation,
            ground_radiation       = temps_out[1].flows.ground_radiation,
        ),
        ventral = (;
            solar                  = temps_out[2].flows.solar,
            convection             = temps_out[2].flows.convection,
            conduction             = temps_out[2].flows.conduction,
            net_metabolic          = temps_out[2].flows.net_metabolic,
            skin_evaporation       = temps_out[2].flows.skin_evaporation,
            insulation_evaporation = temps_out[2].flows.insulation_evaporation,
            longwave               = temps_out[2].flows.longwave,
            sky_radiation          = temps_out[2].flows.sky_radiation,
            bush_radiation         = temps_out[2].flows.bush_radiation,
            vegetation_radiation   = temps_out[2].flows.vegetation_radiation,
            ground_radiation       = temps_out[2].flows.ground_radiation,
        ),
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
    return (; thermoregulation, morphology, energy_flows, mass_flows)
end
