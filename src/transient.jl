# Transient (ODE, lumped-capacitance) body-temperature models.
#
# This file currently holds two parallel formulations, pending a follow-up pass that
# unifies them onto the same steady-state-derived physics:
# - Ectotherm-style (`ectotherm_onelump`/`ectotherm_twolump`): a near-literal port of
#   NicheMapR's onelump.R/onelump_var.R/twolump.R, with its own linearized
#   solar/convection/radiation physics and a bare-`body` calling convention.
# - Endotherm-style (`endotherm_onelump`): reuses the steady-state `heat_balance`/
#   `_pack_sides` residuals directly as the ODE numerator, `Organism`-based calling
#   convention.

_resolve_metabolic_heat(value, core_temperature) = value
_resolve_metabolic_heat(f::Function, core_temperature) = f(core_temperature)

# ---------------------------------------------------------------------------
# Ectotherm: onelump / twolump (NicheMapR-derived, linearized physics)
# ---------------------------------------------------------------------------

_lumped_characteristic_dimension_formula(shape) = VolumeCubeRoot()
_lumped_characteristic_dimension_formula(::Cylinder) = ScaledDimension(1.0, :length_skin)

_lumped_internal_gradient_correction(::Ellipsoid, shape_factor, flesh_conductivity, volumetric_metabolic_heat) =
    (convection = volumetric_metabolic_heat * shape_factor / (2 * flesh_conductivity),
     radiation = volumetric_metabolic_heat * shape_factor / (2 * flesh_conductivity))
_lumped_internal_gradient_correction(shape, shape_factor, flesh_conductivity, volumetric_metabolic_heat) =
    (convection = volumetric_metabolic_heat * shape_factor / (4 * flesh_conductivity),
     radiation = volumetric_metabolic_heat * shape_factor / (2 * flesh_conductivity))

"""Solar/convective/radiative flows shared by [`ectotherm_onelump`](@ref) and [`ectotherm_twolump`](@ref)."""
function _lumped_flows(
    body::AbstractBody, posture, environment_pars, environment_vars,
    surface_temperature, radiation_linearisation_temperature;
    body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
    smoothing::SmoothingStrategy=HardBound(),
)
    total_area = BiophysicalGeometry.total_area(body)

    absorptivities = Absorptivities(; body=DorsalVentral(body_absorptivity, body_absorptivity), ground=environment_pars.ground_albedo)
    view_factors = ViewFactors(sky_view_factor, ground_view_factor, 0.0, 0.0)
    solar_conditions = SolarConditions(environment_vars)
    silhouette_area = BiophysicalGeometry.silhouette_area(body, posture)
    solar_flow = solar(body, absorptivities, view_factors, solar_conditions, silhouette_area, 0.0u"m^2").solar_flow

    convection_out = convection(;
        body, area=total_area,
        air_temperature=environment_vars.air_temperature,
        surface_temperature,
        wind_speed=environment_vars.wind_speed,
        atmospheric_pressure=environment_vars.atmospheric_pressure,
        fluid=environment_pars.fluid,
        gas_fractions=environment_pars.gas_fractions,
        convection_enhancement=environment_pars.convection_enhancement,
        characteristic_dimension_formula=_lumped_characteristic_dimension_formula(shape(body)),
        smoothing,
    )
    convective_heat_transfer_coefficient = convection_out.heat_transfer_coefficient.combined

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    radiant_temperature = (environment_vars.sky_temperature + environment_vars.ground_temperature) / 2
    radiative_heat_transfer_coefficient = 4 * emissivity * σ * ((radiation_linearisation_temperature + radiant_temperature) / 2)^3

    return (; solar_flow, convective_heat_transfer_coefficient, radiative_heat_transfer_coefficient, total_area, radiant_temperature)
end

"""
    ectotherm_onelump(core_temperature, t, body, environment_pars, environment_vars;
            internal_conduction, posture=Intermediate(),
            body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
            metabolic_heat_volumetric, smoothing=HardBound())

One-lump transient body-temperature derivative `dTc/dt` under a time-varying environment
(`environment_vars` already evaluated at `t`). `metabolic_heat_volumetric` (W/m³) may be
a `Quantity` or a `Function` of `core_temperature`. Port of NicheMapR's `onelump_var.R`.

# Returns
`dTc/dt` (K/s).

# References
Kearney, Michael R., Warren P. Porter, and Raymond B. Huey. 2021. “Modelling the Joint Effects of
 Body Size and Microclimate on Heat Budgets and Foraging Opportunities of Ectotherms.”
 Methods in Ecology and Evolution 12 (13): 458–67. https://doi.org/10.1111/2041-210X.13528.
"""
function ectotherm_onelump(
    core_temperature, t, body::AbstractBody, environment_pars, environment_vars;
    internal_conduction::InternalConductionParameters,
    posture=Intermediate(), body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
    metabolic_heat_volumetric, smoothing::SmoothingStrategy=HardBound(),
)
    mass = shape(body).mass
    volume = body.geometry.volume
    thermal_capacitance = mass * internal_conduction.flesh_specific_heat
    volumetric_metabolic_heat = _resolve_metabolic_heat(metabolic_heat_volumetric, core_temperature)
    metabolic_heat_flow = volumetric_metabolic_heat * volume

    flows = _lumped_flows(
        body, posture, environment_pars, environment_vars,
        core_temperature + 0.1u"K", core_temperature;
        body_absorptivity, emissivity, sky_view_factor, ground_view_factor, smoothing,
    )
    (; solar_flow, convective_heat_transfer_coefficient, radiative_heat_transfer_coefficient, total_area, radiant_temperature) = flows

    shape_factor = internal_gradient_shape_factor(body)
    correction = _lumped_internal_gradient_correction(shape(body), shape_factor, internal_conduction.flesh_conductivity, volumetric_metabolic_heat)

    forcing_rate = (solar_flow + metabolic_heat_flow +
        convective_heat_transfer_coefficient * total_area * (correction.convection + environment_vars.air_temperature) +
        radiative_heat_transfer_coefficient * total_area * (correction.radiation + radiant_temperature)) / thermal_capacitance
    decay_rate = total_area * (convective_heat_transfer_coefficient + radiative_heat_transfer_coefficient) / thermal_capacitance

    return forcing_rate - decay_rate * core_temperature
end

"""
    ectotherm_onelump(t::AbstractVector, core_temperature_init, body, environment_pars, environment_vars;
            internal_conduction, posture=Intermediate(),
            body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
            metabolic_heat_volumetric, smoothing=HardBound())

Closed-form one-lump transient body temperature for a **constant** environment.
Port of NicheMapR's `onelump.R`.

# Returns
NamedTuple with `core_temperature` (vector, one per `t`), `final_core_temperature`,
`time_constant`, `initial_rate`.

# References
Kearney, Michael R., Warren P. Porter, and Raymond B. Huey. 2021. “Modelling the Joint Effects of
 Body Size and Microclimate on Heat Budgets and Foraging Opportunities of Ectotherms.”
 Methods in Ecology and Evolution 12 (13): 458–67. https://doi.org/10.1111/2041-210X.13528.
"""
function ectotherm_onelump(
    t::AbstractVector, core_temperature_init, body::AbstractBody, environment_pars, environment_vars;
    internal_conduction::InternalConductionParameters,
    posture=Intermediate(), body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
    metabolic_heat_volumetric, smoothing::SmoothingStrategy=HardBound(),
)
    mass = shape(body).mass
    volume = body.geometry.volume
    thermal_capacitance = mass * internal_conduction.flesh_specific_heat
    volumetric_metabolic_heat = _resolve_metabolic_heat(metabolic_heat_volumetric, core_temperature_init)
    metabolic_heat_flow = volumetric_metabolic_heat * volume

    flows = _lumped_flows(
        body, posture, environment_pars, environment_vars,
        core_temperature_init + 0.1u"K", core_temperature_init;
        body_absorptivity, emissivity, sky_view_factor, ground_view_factor, smoothing,
    )
    (; solar_flow, convective_heat_transfer_coefficient, radiative_heat_transfer_coefficient, total_area, radiant_temperature) = flows

    shape_factor = internal_gradient_shape_factor(body)
    correction = _lumped_internal_gradient_correction(shape(body), shape_factor, internal_conduction.flesh_conductivity, volumetric_metabolic_heat)

    forcing_rate = (solar_flow + metabolic_heat_flow +
        convective_heat_transfer_coefficient * total_area * (correction.convection + environment_vars.air_temperature) +
        radiative_heat_transfer_coefficient * total_area * (correction.radiation + radiant_temperature)) / thermal_capacitance
    decay_rate = total_area * (convective_heat_transfer_coefficient + radiative_heat_transfer_coefficient) / thermal_capacitance

    final_core_temperature = forcing_rate / decay_rate
    time_constant = 1 / decay_rate
    initial_rate = forcing_rate - decay_rate * core_temperature_init
    core_temperature = @. (core_temperature_init - final_core_temperature) * exp(-decay_rate * t) + final_core_temperature

    return (; core_temperature, final_core_temperature, time_constant, initial_rate)
end

function _twolump_core_geometry(::Ellipsoid, body::AbstractBody, shell_thickness)
    a, b, c = body.geometry.length[1], body.geometry.length[2], body.geometry.length[3]
    core_a, core_b, core_c = a - shell_thickness, b - shell_thickness, c - shell_thickness
    core_volume = (4 / 3) * π * core_a * core_b * core_c
    core_a_m, core_b_m, core_c_m = ustrip(u"m", core_a), ustrip(u"m", core_b), ustrip(u"m", core_c)
    eccentricity = sqrt(core_a_m^2 - core_c_m^2) / core_a_m
    core_area = BiophysicalGeometry.surface_area(shape(body), core_a_m, core_b_m, core_c_m, eccentricity)
    core_characteristic_radius = min(core_a, core_b, core_c)
    return (; core_volume, core_area, core_characteristic_radius)
end

function _twolump_core_geometry(::Cylinder, body::AbstractBody, shell_thickness)
    radius_skin = body.geometry.length.radius_skin
    length_skin = body.geometry.length.length_skin
    core_radius = radius_skin - shell_thickness
    core_length = length_skin - shell_thickness
    core_volume = π * core_radius^2 * core_length
    core_area = 2π * core_radius^2 + 2π * core_radius * core_length
    core_characteristic_radius = min(radius_skin, length_skin / 2) - shell_thickness
    return (; core_volume, core_area, core_characteristic_radius)
end

"""
    ectotherm_twolump(state::NamedTuple{(:core_temperature,:shell_temperature)}, t, body, environment_pars,
            environment_vars; internal_conduction, shell_thickness, posture=Intermediate(),
            body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
            metabolic_heat_volumetric, smoothing=HardBound())

Two-lump (core + shell) transient body-temperature derivatives, `Cylinder`/`Ellipsoid`
bodies only. The core compartment uses `internal_conduction`'s flesh conductivity/specific
heat; the shell uses its fat conductivity/specific heat. Port of NicheMapR's `twolump.R`.

Deviates from the R reference in one respect: `surface_temperature` is solved
algebraically from `shell_temperature` each call (eq. 64) rather than integrated as a
third ODE pseudo-state via a solver relaxation trick that exists there only to give
`deSolve` something to converge on.

# Returns
NamedTuple with `core_temperature_rate`, `shell_temperature_rate` (K/s), `surface_temperature`
(algebraic), `final_core_temperature` (closed-form steady state).

# References
Kearney, Michael R., Warren P. Porter, and Raymond B. Huey. 2021. “Modelling the Joint Effects of
 Body Size and Microclimate on Heat Budgets and Foraging Opportunities of Ectotherms.”
 Methods in Ecology and Evolution 12 (13): 458–67. https://doi.org/10.1111/2041-210X.13528.
"""
function ectotherm_twolump(
    state::NamedTuple{(:core_temperature, :shell_temperature)}, t, body::AbstractBody, environment_pars, environment_vars;
    internal_conduction::InternalConductionParameters, shell_thickness,
    posture=Intermediate(), body_absorptivity, emissivity, sky_view_factor, ground_view_factor,
    metabolic_heat_volumetric, smoothing::SmoothingStrategy=HardBound(),
)
    (; core_temperature, shell_temperature) = state
    density = shape(body).density
    volume = body.geometry.volume
    volumetric_metabolic_heat = _resolve_metabolic_heat(metabolic_heat_volumetric, core_temperature)
    metabolic_heat_flow = volumetric_metabolic_heat * volume

    core_geometry = _twolump_core_geometry(shape(body), body, shell_thickness)
    shell_volume = volume - core_geometry.core_volume
    shell_capacitance = shell_volume * density * internal_conduction.fat_specific_heat
    core_capacitance = core_geometry.core_volume * density * internal_conduction.flesh_specific_heat
    core_shell_resistance = core_geometry.core_characteristic_radius / (internal_conduction.flesh_conductivity * core_geometry.core_area)

    flows = _lumped_flows(
        body, posture, environment_pars, environment_vars, shell_temperature, shell_temperature;
        body_absorptivity, emissivity, sky_view_factor, ground_view_factor, smoothing,
    )
    (; solar_flow, convective_heat_transfer_coefficient, radiative_heat_transfer_coefficient, total_area, radiant_temperature) = flows
    convective_resistance = 1 / (convective_heat_transfer_coefficient * total_area)
    radiative_resistance = 1 / (radiative_heat_transfer_coefficient * total_area)
    shell_resistance = (shell_thickness / 2) / (internal_conduction.fat_conductivity * total_area)

    air_temperature = environment_vars.air_temperature
    shell_balance_factor = shell_resistance / (2 * radiative_resistance) + shell_resistance / (2 * convective_resistance) + 1
    surface_temperature = (solar_flow + radiant_temperature / radiative_resistance + air_temperature / convective_resistance + 2 * shell_temperature / shell_resistance) /
        (1 / radiative_resistance + 1 / convective_resistance + 2 / shell_resistance)

    core_relaxation_rate = 1 / (core_capacitance * core_shell_resistance)
    core_metabolic_heating_rate = metabolic_heat_flow / core_capacitance
    shell_relaxation_rate = (1 / shell_capacitance) * (1 / core_shell_resistance + 2 / shell_resistance - (2 / shell_resistance) / shell_balance_factor)
    shell_core_coupling_rate = 1 / (core_shell_resistance * shell_capacitance)
    shell_environment_forcing_rate = (solar_flow + radiant_temperature / radiative_resistance + air_temperature / convective_resistance) / (shell_balance_factor * shell_capacitance)
    final_core_temperature = -(core_relaxation_rate * shell_environment_forcing_rate + shell_relaxation_rate * core_metabolic_heating_rate) /
        (core_relaxation_rate * (shell_core_coupling_rate - shell_relaxation_rate))

    core_temperature_rate = (metabolic_heat_flow - (core_temperature - shell_temperature) / core_shell_resistance) / core_capacitance
    shell_temperature_rate = ((core_temperature - shell_temperature) / core_shell_resistance - (shell_temperature - surface_temperature) / (shell_resistance / 2)) / shell_capacitance

    return (; core_temperature_rate, shell_temperature_rate, surface_temperature, final_core_temperature)
end

# ---------------------------------------------------------------------------
# Endotherm: onelump (reuses steady-state heat_balance / _pack_sides residuals)
# ---------------------------------------------------------------------------

_thermal_capacitance(organism::Organism) =
    flesh_volume(body(organism)) * shape(body(organism)).density *
    conduction_pars_internal(organism).flesh_specific_heat

"""
    endotherm_onelump(core_temperature, t, organism::Organism, e; kw...)

One-lump transient core-temperature derivative for an endotherm, reusing the existing
steady-state heat-balance physics as the ODE numerator. Dispatches on `insulation(body(organism))`:

- `Naked`: reuses `heat_balance(core_temperature, organism, e)` verbatim; metabolic heat comes
  from `organism`'s own `metabolism_pars.model` (a genuine forward model, not solved for).
- `FibrousLayer`/`CompositeInsulation`: needs a `metabolic_heat_flow` keyword (`Quantity` or
  `Function(core_temperature)`) since the steady-state path only has this as a zbrent-solved
  unknown, not a forward model. Skin/insulation temperature are solved algebraically each call
  via the existing `_pack_sides` machinery (not a second ODE state).

All other thermoregulatory effectors are read from `organism`'s traits and held fixed.

# Returns
NamedTuple with `core_temperature_rate` (K/s), `skin_temperature`, `insulation_temperature`,
`energy_flows`, `mass_flows`.
"""
endotherm_onelump(core_temperature, t, organism::Organism, e; kw...) =
    endotherm_onelump(core_temperature, t, insulation(body(organism)), organism, e; kw...)

function endotherm_onelump(
    core_temperature, t, ::Naked, organism::Organism, e; smoothing::SmoothingStrategy=HardBound(),
)
    out = heat_balance(core_temperature, organism, e; smoothing)
    core_temperature_rate = out.energy_balance.heat_balance / _thermal_capacitance(organism)
    return (;
        core_temperature_rate,
        skin_temperature=out.skin_temperature,
        insulation_temperature=out.insulation_temperature,
        energy_flows=out.energy_balance,
        mass_flows=out.mass_balance,
    )
end

function endotherm_onelump(
    core_temperature, t, ::Union{FibrousLayer,CompositeInsulation}, organism::Organism, e;
    metabolic_heat_flow,
    skin_temperature_guess=core_temperature - 3u"K",
    insulation_temperature_guess=e.environment_vars.air_temperature,
    minimum_metabolic_heat=metabolism_pars(organism).metabolic_heat_flow,
    smoothing::SmoothingStrategy=HardBound(),
)
    mhf = _resolve_metabolic_heat(metabolic_heat_flow, core_temperature)

    packed = _pack_sides(organism, e, core_temperature, skin_temperature_guess, insulation_temperature_guess; smoothing)
    (; temps_out, sky_factor_ref, vegetation_factor_ref) = packed
    dmult = sky_factor_ref + vegetation_factor_ref
    vmult = 1 - dmult

    net_metabolic_heat_internal = temps_out[1].flows.net_metabolic * dmult + temps_out[2].flows.net_metabolic * vmult
    skin_temperature = temps_out[1].skin_temperature * dmult + temps_out[2].skin_temperature * vmult
    lung_temperature = (core_temperature + skin_temperature) / 2

    environment_pars = stripparams(e.environment_pars)
    resp_atmos = AtmosphericConditions(e.environment_vars)
    respiration_out = respiration(
        MetabolicRates(; metabolic=mhf, sum=net_metabolic_heat_internal, minimum=minimum_metabolic_heat),
        respiration_pars(organism), resp_atmos, body(organism).shape.mass, lung_temperature,
        e.environment_vars.air_temperature;
        gas_fractions=environment_pars.gas_fractions, O2conversion=Kleiber1961(), smoothing,
    )
    respiration_heat_flow = respiration_out.respiration_heat_flow

    core_temperature_rate = (mhf - respiration_heat_flow - net_metabolic_heat_internal) / _thermal_capacitance(organism)

    out = _assemble_multisided_output(organism, e, core_temperature, mhf, respiration_out, packed; smoothing)
    return (;
        core_temperature_rate,
        skin_temperature=out.thermoregulation.skin_temperature,
        insulation_temperature=out.thermoregulation.insulation_temperature,
        net_metabolic_heat_internal,
        energy_flows=out.energy_flows,
        mass_flows=out.mass_flows,
    )
end
