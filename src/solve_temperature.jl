# Passive thermoregulation solver (ectotherms, torpid endotherms, leaves)

"""
    solve_temperature(organism, environment; temperature_bracket=(270.0u"K", 370.0u"K"))

Find the steady-state temperature of an organism (animal or leaf) by root-finding on its
heat balance. A naked animal or a leaf returns the output of `heat_balance`; an insulated
animal is solved by its [`sides`](@ref) and returns a [`ThermoregulationOutput`](@ref).

Returns air temperature as the fallback `core_temperature` if root-finding fails.
"""
function solve_temperature(organism, environment; temperature_bracket=(270.0u"K", 370.0u"K"),
                            smoothing::SmoothingStrategy=HardBound())
    solve_temperature(insulation(body(organism)), organism, environment; temperature_bracket, smoothing)
end

function solve_temperature(::Naked, organism, environment; temperature_bracket=(270.0u"K", 370.0u"K"),
                            smoothing::SmoothingStrategy=HardBound())
    lo = ustrip(u"K", temperature_bracket[1])
    hi = ustrip(u"K", temperature_bracket[2])
    fallback_temperature = environment.environment_vars.air_temperature
    try
        solved_temperature = zbrent(
            T -> ustrip(u"W", heat_balance(T * u"K", organism, environment; smoothing).heat_balance),
            lo, hi, 1e-3,
        )
        heat_balance(solved_temperature * u"K", organism, environment; smoothing)
    catch
        heat_balance(fallback_temperature, organism, environment; smoothing)
    end
end

function solve_temperature(::Union{FibrousLayer,CompositeInsulation}, organism, environment;
                            temperature_bracket=(270.0u"K", 370.0u"K"), smoothing::SmoothingStrategy=HardBound())
    environment_pars = stripparams(environment.environment_pars)
    environment_vars = environment.environment_vars
    resp_pars = respiration_pars(organism)
    metab_pars = metabolism_pars(organism)
    tolerance = options(organism).temperature_error_tolerance
    lung_mass = mass(organism.body.shape)
    atmos = AtmosphericConditions(environment_vars)
    parts = sides(organism, environment; smoothing)

    initial_skin_temperature = metab_pars.core_temperature - 5u"K"
    initial_insulation_temperature = metab_pars.core_temperature - 3u"K"
    surfaces(core_temperature) = _solve_surfaces(parts, map(part -> part.neighbours, parts), core_temperature,
        initial_skin_temperature, initial_insulation_temperature, tolerance, smoothing)
    respiration_at(core_temperature, metabolic_heat_flow, lung_temperature) = respiration(
        MetabolicRates(; metabolic = metabolic_heat_flow, sum = metabolic_heat_flow, minimum = metabolic_heat_flow),
        resp_pars, atmos, lung_mass, lung_temperature, environment_vars.air_temperature;
        gas_fractions = environment_pars.gas_fractions, O2conversion = Kleiber1961(), smoothing,
    )

    function residual(temperature)
        core_temperature = temperature * u"K"
        solved = surfaces(core_temperature)
        net_metabolic = solved.dorsal.net_metabolic * parts.dorsal.weight + solved.ventral.net_metabolic * parts.ventral.weight
        lung_temperature = (core_temperature + (solved.dorsal.skin_temperature + solved.ventral.skin_temperature) * 0.5) * 0.5
        metabolic_heat_flow = metabolic_rate(metab_pars.model, lung_mass, core_temperature)
        respiration_heat_flow = respiration_at(core_temperature, metabolic_heat_flow, lung_temperature).respiration_heat_flow
        ustrip(u"W", metabolic_heat_flow - respiration_heat_flow - net_metabolic)
    end

    lo = ustrip(u"K", temperature_bracket[1])
    hi = ustrip(u"K", temperature_bracket[2])
    equilibrium_temperature = try
        zbrent(residual, lo, hi, 1e-3) * u"K"
    catch
        environment_vars.air_temperature
    end

    solved = surfaces(equilibrium_temperature)
    skin_mean = solved.dorsal.skin_temperature * parts.dorsal.weight + solved.ventral.skin_temperature * parts.ventral.weight
    lung_temperature = (equilibrium_temperature + skin_mean) * 0.5
    metabolic_heat_flow = metabolic_rate(metab_pars.model, lung_mass, equilibrium_temperature)
    respiration_out = respiration_at(equilibrium_temperature, metabolic_heat_flow, lung_temperature)
    _organism_output(organism, environment, equilibrium_temperature, respiration_out.metabolic_heat_flow,
        respiration_out, parts, solved; smoothing)
end
