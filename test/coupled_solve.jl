using HeatExchange
using BiophysicalGeometry
using Unitful, UnitfulMoles
using FluidProperties
using Test

using HeatExchange: solve_coupled_metabolic_rate,
    EnvironmentTemperatures, ViewFactors, AtmosphericConditions,
    example_environment_vars, example_environment_pars, example_respiration_pars

# Two half-cylinders joined flat to flat, as a torso's dorsal and ventral halves.
ρ = 1000.0u"kg/m^3"
fat = FatLayer(0.10, 901.0u"kg/m^3")
fur = FibrousLayer(0.01u"m", 3e-5u"m", 3e7u"1/m^2")

part_fibres = FibreProperties(;
    diameter = 3e-5u"m", length = 0.01u"m", density = 3e7u"1/m^2",
    depth = 0.01u"m", reflectance = 0.2, conductivity = 0.209u"W/m/K",
)
part_insulation_pars = InsulationParameters(;
    dorsal = part_fibres, ventral = part_fibres,
    depth_compressed = 0.01u"m", longwave_depth_fraction = 1.0,
)

env_vars = example_environment_vars(; air_temperature = u"K"(20.0u"°C"), wind_speed = 1.0u"m/s")
env_pars = example_environment_pars()
core_temperature = u"K"(37.0u"°C")

traits = (;
    core_temperature,
    flesh_conductivity = 0.5u"W/m/K",
    fat_conductivity   = 0.2u"W/m/K",
    ϵ_body             = 0.99,
    skin_wetness       = 0.01,
    insulation_wetness = 0.0,
    bare_skin_fraction = 0.0,
    eye_fraction       = 0.0,
)

packed_environment = (;
    temperature = EnvironmentTemperatures(
        env_vars.air_temperature, env_vars.sky_temperature,
        env_vars.ground_temperature, env_vars.vegetation_temperature,
        env_vars.bush_temperature, env_vars.substrate_temperature,
    ),
    view_factors = ViewFactors(0.5, 0.5, 0.0, 0.0),
    atmos = AtmosphericConditions(env_vars),
    fluid = env_pars.fluid,
    solar_flow = 0.0u"W",
    gas_fractions = env_pars.gas_fractions,
    convection_enhancement = env_pars.convection_enhancement,
)

half(mass) = Body(HalfCylinder(; mass, density = ρ, axis_ratio_b = 3.0), CompositeInsulation(fur, fat))
flat_join = Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover()))
torso = CompositeBody(; parts = (; dorsal = half(5.0u"kg"), ventral = half(5.0u"kg")), joins = (flat_join,))
setup = (;
    insulation_pars = part_insulation_pars,
    traits,
    environment_vars = packed_environment,
    conduction_fraction = 0.0,
    conductance_coefficient = 0.0u"W/K",
    ventral_fraction = 0.5,
    longwave_depth_fraction = 1.0,
)
setups = (; dorsal = setup, ventral = setup)

resp_pars = example_respiration_pars()
atmos = AtmosphericConditions(env_vars)

common = (;
    core_temperature,
    skin_temperature = core_temperature - 5u"K",
    insulation_temperature = env_vars.air_temperature + 2u"K",
    temperature_tolerance = 1e-3u"K",
    respiration_pars = resp_pars,
    lung_mass = 10.0u"kg",
    air_temperature = env_vars.air_temperature,
    atmos,
    gas_fractions = env_pars.gas_fractions,
    metabolic_heat_flow_setpoint = 10.0u"W",
    resp_tolerance = 1e-5,
)

@testset "solve_coupled_metabolic_rate — mechanics and self-consistency" begin
    result = solve_coupled_metabolic_rate(torso, setups; respire = true, common...)
    @test length(result.parts) == 2
    @test all(part -> part.success, result.parts)
    # net_metabolic_total is the plain sum of per-part core→skin flow
    @test result.net_metabolic_total ≈ sum(part -> part.net_metabolic, result.parts)
    # both parts identical → each carries half the total
    @test result.parts[1].net_metabolic ≈ result.parts[2].net_metabolic
    # metabolic heat flow is positive and respiration was applied
    @test result.metabolic_heat_flow > 0.0u"W"
    @test result.respiration_out !== nothing
    # mean skin sits between air and core
    @test env_vars.air_temperature < result.skin_temperature < core_temperature
    @test result.lung_temperature ≈ (core_temperature + result.skin_temperature) * 0.5
end

@testset "solve_coupled_metabolic_rate — respire=false is the identity" begin
    result = solve_coupled_metabolic_rate(torso, setups; respire = false, common...)
    # With no respiration, the metabolic heat flow is exactly the internal heat sum
    @test result.metabolic_heat_flow == result.net_metabolic_total
    @test result.respiration_out === nothing
end

@testset "solve_coupled_metabolic_rate — one part" begin
    alone = CompositeBody(; parts = (; dorsal = half(10.0u"kg")), joins = ())
    result = solve_coupled_metabolic_rate(alone, (; dorsal = setup); respire = false, common...)
    @test result.net_metabolic_total > 0.0u"W"
end

@testset "solve_coupled_metabolic_rate — neighbours from silhouette_factors" begin
    apart = solve_coupled_metabolic_rate(torso, setups; respire = false, common...)
    neighbours = map(view -> view.neighbours, silhouette_factors(torso, Sky(0.5); ndirections = 64, resolution = 32))
    beside = solve_coupled_metabolic_rate(torso, setups; neighbours, respire = false, common...)
    @test all(part -> part.success, beside.parts)
    # The halves of a cylinder hardly see each other.
    @test beside.net_metabolic_total ≈ apart.net_metabolic_total rtol = 0.01
end
