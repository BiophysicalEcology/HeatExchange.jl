# Transient (ODE, lumped-capacitance) body-temperature models.
#
# One physics core, shared with steady state via `heat_balance`/`_radiative_convective_flows`.
# `onelump` dispatches on `insulation(body(organism))` like `heat_balance` itself; `twolump`
# (Naked/`Cylinder`/`Ellipsoid`-only) takes a caller-selectable `SurfaceSolveStrategy`.

_resolve_metabolic_heat(value, core_temperature) = value
_resolve_metabolic_heat(f::Function, core_temperature) = f(core_temperature)

"""
    SurfaceSolveStrategy

Policy for how [`twolump`](@ref) resolves its shell-surface temperature each call.
Mirrors [`SmoothingStrategy`](@ref) in style.

- [`LinearizedSurface`](@ref) (default): fast, no root-finding, has `final_core_temperature`.
- [`RootFindSurface`](@ref): exact root-find each call, no `final_core_temperature`.
"""
abstract type SurfaceSolveStrategy end

"""
    LinearizedSurface() <: SurfaceSolveStrategy

Default `twolump` surface-solve strategy: linearizes solar/longwave/convective exchange
about the current shell temperature. No root-finding; yields a closed-form
`final_core_temperature`.
"""
struct LinearizedSurface <: SurfaceSolveStrategy end

"""
    RootFindSurface() <: SurfaceSolveStrategy

`twolump` surface-solve strategy that root-finds (`zbrent`) the exact nonlinear
surface-temperature residual each call. No closed-form `final_core_temperature`.
"""
struct RootFindSurface <: SurfaceSolveStrategy end

"""
    _with_solar_orientation(organism::Organism, posture)

Copy of `organism` with `radiation_pars.solar_orientation` overridden to `posture`.
Kept as an explicit kwarg (rather than a trait override in `heat_balance.jl`) because
behavioral drivers vary orientation per call in a way a static trait can't capture.

Overrides via the accessor interface, not `setproperties` on `traits(organism)` directly —
composed traits types don't expose `radiation_pars` as a top-level field.
"""
struct _SolarOrientationOverride{T<:AbstractFunctionalTraits,S} <: AbstractFunctionalTraits
    traits::T
    solar_orientation::S
end
shape_pars(t::_SolarOrientationOverride) = shape_pars(t.traits)
insulation_pars(t::_SolarOrientationOverride) = insulation_pars(t.traits)
conduction_pars_external(t::_SolarOrientationOverride) = conduction_pars_external(t.traits)
conduction_pars_internal(t::_SolarOrientationOverride) = conduction_pars_internal(t.traits)
convection_pars(t::_SolarOrientationOverride) = convection_pars(t.traits)
radiation_pars(t::_SolarOrientationOverride) =
    setproperties(radiation_pars(t.traits), (; solar_orientation=t.solar_orientation))
evaporation_pars(t::_SolarOrientationOverride) = evaporation_pars(t.traits)
hydraulic_pars(t::_SolarOrientationOverride) = hydraulic_pars(t.traits)
respiration_pars(t::_SolarOrientationOverride) = respiration_pars(t.traits)
metabolism_pars(t::_SolarOrientationOverride) = metabolism_pars(t.traits)
options(t::_SolarOrientationOverride) = options(t.traits)

_with_solar_orientation(organism::Organism, posture) =
    Organism(body(organism), _SolarOrientationOverride(traits(organism), posture))

# ---------------------------------------------------------------------------
# onelump
# ---------------------------------------------------------------------------

_thermal_capacitance(organism::Organism) =
    flesh_volume(body(organism)) * shape(body(organism)).density *
    conduction_pars_internal(organism).flesh_specific_heat

"""
    onelump(core_temperature, t, organism::Organism, e; kw...)

One-lump transient core-temperature derivative, reusing the steady-state heat-balance
physics as the ODE numerator. Dispatches on `insulation(body(organism))`:

- `Naked`: `heat_balance(core_temperature, organism, e)` (with `posture` applied); metabolic
  heat is a forward model from `organism`'s own `metabolism_pars.model`; evaporative and
  respiratory heat loss are included.
- `FibrousLayer`/`CompositeInsulation`: needs a `metabolic_heat_flow` keyword (`Quantity`
  or `Function(core_temperature)`), since the steady-state path only has this as a
  zbrent-solved unknown. Skin/insulation temperature are solved algebraically each call
  via `_pack_sides` (not a second ODE state). `posture` is accepted but unused here, for
  signature symmetry with the `Naked` branch.

# Returns
NamedTuple with `core_temperature_rate` (K/s), `skin_temperature`, `insulation_temperature`,
`energy_flows`, `mass_flows`.
"""
onelump(core_temperature, t, organism::Organism, e; kw...) =
    onelump(core_temperature, t, insulation(body(organism)), organism, e; kw...)

function onelump(
    core_temperature, t, ::Naked, organism::Organism, e;
    posture=radiation_pars(organism).solar_orientation, smoothing::SmoothingStrategy=HardBound(),
)
    out = heat_balance(core_temperature, _with_solar_orientation(organism, posture), e; smoothing)
    core_temperature_rate = out.energy_balance.heat_balance / _thermal_capacitance(organism)
    return (;
        core_temperature_rate,
        skin_temperature=out.skin_temperature,
        insulation_temperature=out.insulation_temperature,
        energy_flows=out.energy_balance,
        mass_flows=out.mass_balance,
    )
end

function onelump(
    core_temperature, t, ::Union{FibrousLayer,CompositeInsulation}, organism::Organism, e;
    metabolic_heat_flow, # TODO: rename this kwarg (Quantity or Function(core_temperature)) - tracked as a follow-up issue
    posture=radiation_pars(organism).solar_orientation,
    skin_temperature_guess=core_temperature - 3u"K",
    insulation_temperature_guess=e.environment_vars.air_temperature,
    minimum_metabolic_heat=metabolism_pars(organism).metabolic_heat_flow,
    smoothing::SmoothingStrategy=HardBound(),
)
    resolved_metabolic_heat_flow = _resolve_metabolic_heat(metabolic_heat_flow, core_temperature)

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
        MetabolicRates(; metabolic=resolved_metabolic_heat_flow, sum=net_metabolic_heat_internal, minimum=minimum_metabolic_heat),
        respiration_pars(organism), resp_atmos, body(organism).shape.mass, lung_temperature,
        e.environment_vars.air_temperature;
        gas_fractions=environment_pars.gas_fractions, O2conversion=Kleiber1961(), smoothing,
    )
    respiration_heat_flow = respiration_out.respiration_heat_flow

    core_temperature_rate = (resolved_metabolic_heat_flow - respiration_heat_flow - net_metabolic_heat_internal) / _thermal_capacitance(organism)

    out = _assemble_multisided_output(organism, e, core_temperature, resolved_metabolic_heat_flow, respiration_out, packed; smoothing)
    return (;
        core_temperature_rate,
        skin_temperature=out.thermoregulation.skin_temperature,
        insulation_temperature=out.thermoregulation.insulation_temperature,
        net_metabolic_heat_internal,
        energy_flows=out.energy_flows,
        mass_flows=out.mass_flows,
    )
end

# ---------------------------------------------------------------------------
# twolump
# ---------------------------------------------------------------------------

"""
    _linearized_surface_flows(shell_temperature, organism::Organism, environment_pars, environment_vars;
                               posture, smoothing=HardBound())

Solar/longwave/convective flows for [`twolump`](@ref)'s [`LinearizedSurface`](@ref) mode,
built on `_radiative_convective_flows`. The longwave term is linearized about
`shell_temperature` via the exact local slope of `σT⁴` (`d(σT⁴)/dT = 4σT³ = 4·(σT⁴)/T`),
giving an equivalent `radiant_temperature` matching both value and slope there.
"""
function _linearized_surface_flows(shell_temperature, organism::Organism, environment_pars, environment_vars;
                                    posture, smoothing::SmoothingStrategy=HardBound())
    organism_with_posture = _with_solar_orientation(organism, posture)
    flows = _radiative_convective_flows(shell_temperature, organism_with_posture, environment_pars, environment_vars; smoothing)
    (; solar_flow, longwave_flow_in, longwave_flow_out, convection_heat_flow, convection_out, convection_area) = flows
    total_area = BiophysicalGeometry.total_area(organism.body)
    convective_heat_transfer_coefficient = convection_out.heat_transfer_coefficient.combined

    radiative_heat_transfer_coefficient = 4 * longwave_flow_out / shell_temperature / total_area

    radiant_temperature = shell_temperature +
        (longwave_flow_in - longwave_flow_out) / (radiative_heat_transfer_coefficient * total_area)

    return (; solar_flow, convective_heat_transfer_coefficient, radiative_heat_transfer_coefficient,
              total_area, convection_area, radiant_temperature)
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
    twolump(state::NamedTuple{(:core_temperature,:shell_temperature)}, t, organism::Organism, e;
            shell_thickness, posture=radiation_pars(organism).solar_orientation,
            surface_solve::SurfaceSolveStrategy=LinearizedSurface(), smoothing::SmoothingStrategy=HardBound())

Two-lump (core + shell) transient body-temperature derivatives, `Cylinder`/`Ellipsoid`
bodies only. Core uses `conduction_pars_internal(organism)`'s flesh conductivity/specific
heat; shell uses its fat conductivity/specific heat. Surface exchange is computed by
`surface_solve` ([`LinearizedSurface`](@ref) by default, or [`RootFindSurface`](@ref)).

Silently ignores any fur/insulation present — core/shell is a flesh/fat split, orthogonal
to fur. Does **not** include respiration/evaporation.

# Returns
NamedTuple with `core_temperature_rate`, `shell_temperature_rate` (K/s), `surface_temperature`
(algebraic), and — `LinearizedSurface` only — `final_core_temperature` (closed-form steady state).
"""
twolump(state::NamedTuple{(:core_temperature,:shell_temperature)}, t, organism::Organism, e;
    shell_thickness, posture=radiation_pars(organism).solar_orientation,
    surface_solve::SurfaceSolveStrategy=LinearizedSurface(), smoothing::SmoothingStrategy=HardBound()) =
    twolump(surface_solve, state, t, organism, e; shell_thickness, posture, smoothing)

function twolump(
    ::LinearizedSurface, state, t, organism::Organism, e;
    shell_thickness, posture, smoothing::SmoothingStrategy=HardBound(),
)
    (; core_temperature, shell_temperature) = state
    (; environment_pars, environment_vars) = e
    b = body(organism)
    internal_conduction = conduction_pars_internal(organism)
    density = shape(b).density
    volume = b.geometry.volume
    metabolic_heat_flow = metabolic_rate(metabolism_pars(organism).model, shape(b).mass, core_temperature)

    core_geometry = _twolump_core_geometry(shape(b), b, shell_thickness)
    shell_volume = volume - core_geometry.core_volume
    shell_capacitance = shell_volume * density * internal_conduction.fat_specific_heat
    core_capacitance = core_geometry.core_volume * density * internal_conduction.flesh_specific_heat
    core_shell_resistance = core_geometry.core_characteristic_radius / (internal_conduction.flesh_conductivity * core_geometry.core_area)

    flows = _linearized_surface_flows(shell_temperature, organism, environment_pars, environment_vars; posture, smoothing)
    (; solar_flow, convective_heat_transfer_coefficient, radiative_heat_transfer_coefficient, total_area, convection_area, radiant_temperature) = flows
    convective_resistance = 1 / (convective_heat_transfer_coefficient * convection_area)
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

function twolump(
    ::RootFindSurface, state, t, organism::Organism, e;
    shell_thickness, posture, smoothing::SmoothingStrategy=HardBound(),
)
    (; core_temperature, shell_temperature) = state
    (; environment_pars, environment_vars) = e
    b = body(organism)
    internal_conduction = conduction_pars_internal(organism)
    metabolic_heat_flow = metabolic_rate(metabolism_pars(organism).model, shape(b).mass, core_temperature)

    core_geometry = _twolump_core_geometry(shape(b), b, shell_thickness)
    density = shape(b).density
    shell_volume = b.geometry.volume - core_geometry.core_volume
    shell_capacitance = shell_volume * density * internal_conduction.fat_specific_heat
    core_capacitance = core_geometry.core_volume * density * internal_conduction.flesh_specific_heat
    core_shell_resistance = core_geometry.core_characteristic_radius / (internal_conduction.flesh_conductivity * core_geometry.core_area)
    total_area = BiophysicalGeometry.total_area(b)
    shell_resistance = (shell_thickness / 2) / (internal_conduction.fat_conductivity * total_area)

    organism_with_posture = _with_solar_orientation(organism, posture)
    residual(Tsk) = begin
        flows = _radiative_convective_flows(Tsk * u"K", organism_with_posture, environment_pars, environment_vars; smoothing)
        net_env = flows.solar_flow + flows.longwave_flow_in - flows.longwave_flow_out - flows.convection_heat_flow
        conduction_from_shell = (Tsk * u"K" - shell_temperature) / (shell_resistance / 2)
        ustrip(u"W", net_env - conduction_from_shell)
    end
    surface_temperature = try
        zbrent(residual, 240.0, 340.0, 1e-3) * u"K"
    catch
        shell_temperature
    end

    core_temperature_rate = (metabolic_heat_flow - (core_temperature - shell_temperature) / core_shell_resistance) / core_capacitance
    shell_temperature_rate = ((core_temperature - shell_temperature) / core_shell_resistance -
                               (shell_temperature - surface_temperature) / (shell_resistance / 2)) / shell_capacitance
    return (; core_temperature_rate, shell_temperature_rate, surface_temperature)
end
