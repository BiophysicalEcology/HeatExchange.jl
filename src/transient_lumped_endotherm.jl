# Lumped-capacitance transient core temperature for endotherms. Reuses the existing
# steady-state residual physics (heat_balance, _pack_sides) as the ODE numerator instead of
# driving it to zero - see docs/transient_body_temperature.md for the mapping. Thermoregulatory
# effectors (insulation, panting, skin wetness, flesh conductivity) are held fixed; only core
# temperature evolves.

_endotherm_metabolic_heat(f::Function, core_temperature) = f(core_temperature)
_endotherm_metabolic_heat(q, core_temperature) = q

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
    mhf = _endotherm_metabolic_heat(metabolic_heat_flow, core_temperature)

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
