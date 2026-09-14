using HeatExchange
using BiophysicalGeometry
using Unitful
using Test

# Self-consistency checks for onelump: no R reference exists for this (transient
# endotherm body temperature is new, not a NicheMapR port), so the primary check is
# convergence to solve_metabolic_rate's independently-computed equilibrium.

environment_vars = example_environment_vars(;
    air_temperature=u"K"(10.0u"°C"), wind_speed=0.5u"m/s",
)
environment_pars = example_environment_pars()
e = (; environment_pars, environment_vars)

@testset "onelump: Naked" begin
    shape_pars = example_shape_pars(; mass=0.5u"kg")
    body = Body(shape_pars, Naked())
    traits = example_heat_exchange_traits(;
        shape_pars, metabolism_pars=example_metabolism_pars(; core_temperature=u"K"(38.0u"°C"), model=Kleiber()),
    )
    organism = Organism(body, traits)

    core_temperature = u"K"(38.0u"°C")
    out = onelump(core_temperature, 0.0u"s", organism, e)

    # exact cross-check: the Naked branch reuses heat_balance's own residual verbatim
    hb = heat_balance(core_temperature, organism, e)
    @test ustrip(u"W", out.energy_flows.heat_balance) == ustrip(u"W", hb.energy_balance.heat_balance)

    # sign sanity: colder core than the surrounding environment should warm (dT/dt > 0),
    # hotter core should cool (dT/dt < 0)
    rate_hot = onelump(u"K"(45.0u"°C"), 0.0u"s", organism, e).core_temperature_rate
    rate_cold = onelump(u"K"(5.0u"°C"), 0.0u"s", organism, e).core_temperature_rate
    @test rate_hot < 0u"K/s"
    @test rate_cold > 0u"K/s"

    # thermal-mass sanity: larger specific heat slows the rate without changing its sign
    traits_heavy = example_heat_exchange_traits(;
        shape_pars, metabolism_pars=example_metabolism_pars(; core_temperature=u"K"(38.0u"°C"), model=Kleiber()),
        conduction_pars_internal=example_conduction_pars_internal(; flesh_specific_heat=3.0e4u"J/kg/K"),
    )
    organism_heavy = Organism(body, traits_heavy)
    rate_heavy = onelump(core_temperature, 0.0u"s", organism_heavy, e).core_temperature_rate
    rate_normal = out.core_temperature_rate
    @test sign(rate_heavy) == sign(rate_normal)
    @test abs(ustrip(u"K/s", rate_heavy)) < abs(ustrip(u"K/s", rate_normal))
end

@testset "onelump: Insulated" begin
    shape_pars = example_shape_pars(; mass=0.0337u"kg")
    insulation_pars = example_insulation_pars()
    conduction_pars_internal = example_conduction_pars_internal(; fat_fraction=0.05)
    radiation_pars = example_radiation_pars()
    fat = FatLayer(conduction_pars_internal.fat_fraction, conduction_pars_internal.fat_density)
    pven = radiation_pars.ventral_fraction
    mean_depth = insulation_pars.dorsal.depth * (1 - pven) + insulation_pars.ventral.depth * pven
    mean_diam = insulation_pars.dorsal.diameter * (1 - pven) + insulation_pars.ventral.diameter * pven
    mean_density = insulation_pars.dorsal.density * (1 - pven) + insulation_pars.ventral.density * pven
    fur = FibrousLayer(mean_depth, mean_diam, mean_density)
    body = Body(shape_pars, CompositeInsulation(fur, fat))

    metabolism_pars = example_metabolism_pars(;
        core_temperature=u"K"(38.0u"°C"), model=McKechnieWolf(),
        metabolic_heat_flow=metabolic_rate(McKechnieWolf(), shape_pars.mass),
    )
    traits = example_heat_exchange_traits(;
        shape_pars, insulation_pars, conduction_pars_internal, radiation_pars, metabolism_pars,
    )
    organism = Organism(body, traits)
    core_temperature = u"K"(38.0u"°C")

    # steady-state convergence: feeding the exact equilibrium metabolic_heat_flow that
    # solve_metabolic_rate independently finds should give core_temperature_rate ≈ 0
    solved = solve_metabolic_rate(organism, e, core_temperature - 3.0u"K", environment_vars.air_temperature)
    equilibrium_metabolic_heat_flow = solved.energy_flows.metabolic_heat_flow
    out_eq = onelump(core_temperature, 0.0u"s", organism, e; metabolic_heat_flow=equilibrium_metabolic_heat_flow)
    @test ustrip(u"K/hr", out_eq.core_temperature_rate) ≈ 0.0 atol = 1e-3

    # flow-conservation check at that equilibrium point
    ef = out_eq.energy_flows
    @test ustrip(u"W", ef.metabolic_heat_flow - ef.respiration_heat_flow - out_eq.net_metabolic_heat_internal) ≈ 0.0 atol = 1e-6

    # below-equilibrium metabolic_heat_flow should cool the core
    low_metabolic_heat_flow = equilibrium_metabolic_heat_flow * 0.5
    out_low = onelump(core_temperature, 0.0u"s", organism, e; metabolic_heat_flow=low_metabolic_heat_flow)
    @test out_low.core_temperature_rate < 0u"K/s"

    # a Function(core_temperature) metabolic_heat_flow is supported (Q10-style)
    out_fn = onelump(core_temperature, 0.0u"s", organism, e; metabolic_heat_flow=T -> equilibrium_metabolic_heat_flow)
    @test ustrip(u"K/hr", out_fn.core_temperature_rate) ≈ ustrip(u"K/hr", out_eq.core_temperature_rate) rtol = 1e-8
end
