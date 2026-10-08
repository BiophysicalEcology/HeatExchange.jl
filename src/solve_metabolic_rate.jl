# Active thermoregulation solver (endotherms, thermogenic plants)
# Finds the metabolic_heat_flow that closes the heat budget at a fixed core temperature.

# Named wrappers around the four output NamedTuples.
# Wrapping enables Base.show dispatch so that e.g. `result.thermoregulation` displays
# human-readably in the REPL instead of as a raw NamedTuple dump.
# getproperty/propertynames forward transparently, so field access is unchanged.
abstract type AbstractOutputSection end
Base.getproperty(s::AbstractOutputSection, f::Symbol) =
    f === :data ? getfield(s, :data) : getproperty(getfield(s, :data), f)
Base.propertynames(s::AbstractOutputSection, private::Bool=false) =
    propertynames(getfield(s, :data))

struct ThermoregulationState{NT<:NamedTuple} <: AbstractOutputSection; data::NT; end
struct MorphologyState{NT<:NamedTuple}       <: AbstractOutputSection; data::NT; end
struct EnergyFlowState{NT<:NamedTuple}       <: AbstractOutputSection; data::NT; end
struct MassFlowState{NT<:NamedTuple}         <: AbstractOutputSection; data::NT; end

"""
    ThermoregulationOutput(thermoregulation, morphology, energy_flows, mass_flows)

Canonical return type for `solve_metabolic_rate` and the IPOPT/NLP path.
Each field is a named wrapper around a NamedTuple (see `ThermoregulationState`,
`MorphologyState`, `EnergyFlowState`, `MassFlowState`).
Field access is transparent: `result.thermoregulation.core_temperature` works as expected.

- `thermoregulation`: temperature states and effector values at the solution.
- `morphology`: body surface areas and geometric measures.
- `energy_flows`: heat fluxes (solar, longwave, generated, evaporation, etc.).
- `mass_flows`: respiration and sweat mass flows.
"""
struct ThermoregulationOutput{T, M, E, F}
    thermoregulation::T
    morphology::M
    energy_flows::E
    mass_flows::F
end

"""
    SolveMetabolicRateOptions

Options controlling the endotherm metabolic rate solver.

# Fields
- `respire`: Whether to include respiration in the heat balance (default: true)
- `temperature_error_tolerance`: Convergence tolerance for simultaneous temperature solution (default: 1e-3 K)
- `resp_tolerance`: Relative convergence tolerance for respiration root finding (default: 1e-5)
"""
Base.@kwdef struct SolveMetabolicRateOptions{RE,ST,BT} <: AbstractModelParameters
    respire::RE = Param(true)
    temperature_error_tolerance::ST = Param(1e-3u"K")
    resp_tolerance::BT = Param(1e-5)
end

"""
    solve_metabolic_rate(o::Organism, e, skin_temperature, insulation_temperature; smoothing = HardBound())

The metabolic heat flow that holds the core of `o` at the `core_temperature` of its
metabolism, solving its [`sides`](@ref) from the given skin and insulation temperatures.
`e` holds the `environment_pars` and `environment_vars`.

Returns a [`ThermoregulationOutput`](@ref), with the results of each side under `dorsal`
and `ventral`.
"""
function solve_metabolic_rate(o::Organism, e, skin_temperature, insulation_temperature;
                              smoothing::SmoothingStrategy=HardBound())
    parts = sides(o, e; smoothing)
    solution = solve_metabolic_rate(parts, traits(o), e, skin_temperature, insulation_temperature;
        lung_mass = mass(o.body.shape), smoothing)
    return _organism_output(o, e, metabolism_pars(o).core_temperature, solution.metabolic_heat_flow,
        solution.respiration_out, parts, solution.parts; smoothing)
end

"""
    solve_metabolic_rate(parts::NamedTuple, traits, e, skin_temperature, insulation_temperature;
        lung_mass, core_temperature, neighbours = nothing, extra_net_metabolic = 0 W,
        smoothing = HardBound())

The metabolic heat flow that holds every one of the [`Part`](@ref)s at `core_temperature`, by
default that of the metabolism in `traits`. The parts' surfaces are solved from the given
temperatures, each starting from the temperatures of the part before it. Their heat, weighted
by each part's `weight`, is met by metabolism less the respiration of lungs of `lung_mass`.

`neighbours` is keyed like `parts`, each entry the fractions of that part's view taken by every
part, as the `neighbours` of `silhouette_factors`. With them the parts exchange longwave
radiation, and are solved again until their temperatures settle. `extra_net_metabolic` is heat
the parts lose by other routes, such as conduction to other compartments.

Returns `(; metabolic_heat_flow, parts, net_metabolic_total, skin_temperature,
insulation_temperature, lung_temperature, respiration_out)`, where `parts` holds each part's
[`solve_surface`](@ref) result and the temperatures are means over the parts.
"""
function solve_metabolic_rate(parts::NamedTuple, traits::AbstractFunctionalTraits, e,
    skin_temperature, insulation_temperature;
    lung_mass,
    core_temperature = metabolism_pars(traits).core_temperature,
    neighbours = nothing,
    extra_net_metabolic = zero(metabolism_pars(traits).metabolic_heat_flow),
    smoothing::SmoothingStrategy=HardBound(),
)
    environment_pars = stripparams(e.environment_pars)
    environment_vars = e.environment_vars
    opts = options(traits)
    resp_pars = respiration_pars(traits)
    setpoint = metabolism_pars(traits).metabolic_heat_flow
    solved = _solve_surfaces(parts, neighbours, core_temperature, skin_temperature, insulation_temperature,
        opts.temperature_error_tolerance, smoothing)

    net_metabolic_total = sum(map((part, s) -> s.net_metabolic * part.weight, parts, solved)) + extra_net_metabolic
    skin_mean = sum(s -> s.skin_temperature, solved) / length(solved)
    insulation_mean = sum(s -> s.insulation_temperature, solved) / length(solved)
    max_skin_temperature = maximum(s -> s.skin_temperature, solved)
    lung_temperature = (core_temperature + skin_mean) * 0.5

    if opts.respire
        atmos = AtmosphericConditions(environment_vars)
        flux_lower = setpoint * (-2.0)
        flux_upper = max_skin_temperature >= core_temperature ? setpoint * 1.01 : setpoint * 10.0
        respiration_at(metabolic_heat_flow) = respiration(
            MetabolicRates(; metabolic = metabolic_heat_flow, sum = net_metabolic_total, minimum = setpoint),
            resp_pars, atmos, lung_mass, lung_temperature, environment_vars.air_temperature;
            gas_fractions = environment_pars.gas_fractions, O2conversion = Kleiber1961(), smoothing,
        )
        metabolic_heat_flow = zbrent(
            x -> ustrip(u"W", respiration_at(x * u"W").balance),
            ustrip(u"W", flux_lower),
            ustrip(u"W", flux_upper),
            opts.resp_tolerance * ustrip(u"W", setpoint),
        ) * u"W"
        respiration_out = respiration_at(metabolic_heat_flow)
        metabolic_heat_flow = respiration_out.metabolic_heat_flow
    else
        metabolic_heat_flow = net_metabolic_total
        respiration_out = nothing
    end

    return (;
        metabolic_heat_flow,
        parts = solved,
        net_metabolic_total,
        skin_temperature = skin_mean,
        insulation_temperature = insulation_mean,
        lung_temperature,
        respiration_out,
    )
end

_solve_surfaces(parts::NamedTuple{K}, ::Nothing, core_temperature, skin_temperature, insulation_temperature, tolerance, smoothing) where {K} =
    NamedTuple{K}(_solve_in_turn(values(parts), core_temperature, skin_temperature, insulation_temperature, tolerance, smoothing))

_solve_in_turn(::Tuple{}, args...) = ()
function _solve_in_turn(parts::Tuple, core_temperature, skin_temperature, insulation_temperature, tolerance, smoothing)
    solved = solve_surface(first(parts), core_temperature; skin_temperature, insulation_temperature, tolerance, smoothing)
    return (solved, _solve_in_turn(Base.tail(parts), core_temperature,
        solved.skin_temperature, solved.insulation_temperature, tolerance, smoothing)...)
end

# Each part sees its neighbours' surfaces, so the parts are solved again with their
# neighbours at the last surface temperatures until those stop changing.
function _solve_surfaces(parts::NamedTuple, neighbours::NamedTuple, core_temperature, skin_temperature, insulation_temperature, tolerance, smoothing)
    beside(temperatures) = map(parts, neighbours) do part, fractions
        seen = values(map((fraction, temperature) -> (; fraction, temperature), fractions, temperatures))
        solve_surface(part, core_temperature; skin_temperature, insulation_temperature, tolerance, neighbours = seen, smoothing)
    end
    solved = beside(map(_ -> insulation_temperature, parts))
    for _ in 1:49
        temperatures = map(s -> s.insulation_temperature, solved)
        solved = beside(temperatures)
        maximum(map((s, t) -> abs(s.insulation_temperature - t), solved, temperatures)) < tolerance && break
    end
    return solved
end
