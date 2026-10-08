using HeatExchange
using BiophysicalGeometry
using Unitful, UnitfulMoles
using FluidProperties
using Test

using HeatExchange: EnvironmentTemperatures, ViewFactors, AtmosphericConditions, OrganismTemperatures,
    example_environment_vars, example_environment_pars, example_respiration_pars

shape = Ellipsoid(mass = 1.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 2.0, axis_ratio_c = 2.0)
fat = FatLayer(0.10, 901.0u"kg/m^3")
fur = FibrousLayer(0.01u"m", 3e-5u"m", 3e7u"1/m^2")
fibres = FibreProperties(;
    diameter = 3e-5u"m", length = 0.01u"m", density = 3e7u"1/m^2",
    depth = 0.01u"m", reflectance = 0.2, conductivity = 0.209u"W/m/K",
)
insulation = InsulationParameters(; dorsal = fibres, ventral = fibres, depth_compressed = 0.01u"m", longwave_depth_fraction = 1.0)

env_vars = example_environment_vars(; air_temperature = u"K"(20.0u"°C"), wind_speed = 1.0u"m/s", global_radiation = 0.0u"W/m^2")
env_pars = example_environment_pars()
environment = (;
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
traits = (;
    flesh_conductivity = 0.5u"W/m/K",
    fat_conductivity = 0.2u"W/m/K",
    ϵ_body = 0.99,
    skin_wetness = 0.01,
    insulation_wetness = 0.0,
    bare_skin_fraction = 0.0,
    eye_fraction = 0.0,
)
part = Part(; body = Body(shape, CompositeInsulation(fur, fat)), insulation, traits, environment, ventral_fraction = 0.5)

core_temperature = u"K"(37.0u"°C")
start = (; skin_temperature = core_temperature - 5u"K", insulation_temperature = env_vars.air_temperature + 2u"K",
    tolerance = 1e-3u"K")
result = solve_surface(part, core_temperature; start...)

@testset "solve_surface" begin
    @test result.success
    @test env_vars.air_temperature < result.insulation_temperature < result.skin_temperature < core_temperature
    @test result.net_metabolic > 0.0u"W"
    @test result.flesh_conductance * (core_temperature - result.skin_temperature) ≈ result.net_metabolic
    @test unit(result.flesh_conductance) == u"W/K"
    hotter = solve_surface(part, core_temperature + 3u"K"; start...)
    @test hotter.net_metabolic > result.net_metabolic
    @test hotter.skin_temperature > result.skin_temperature
end

@testset "heat_balance vanishes at the solved surface" begin
    temperatures = OrganismTemperatures(core_temperature, result.skin_temperature, result.insulation_temperature)
    balance(metabolic_heat_flow) = heat_balance(temperatures, metabolic_heat_flow, part; respiration_pars = example_respiration_pars())
    b = balance(1.0u"W")
    @test abs(ustrip(u"W", b.surface_balance)) < 1e-2
    @test abs(ustrip(u"K", b.residual_skin_temperature)) < 1e-2
    @test b.net_metabolic_heat_internal ≈ result.net_metabolic rtol = 1e-3
    # Metabolism cancels out of the surface balance.
    @test balance(2.0u"W").surface_balance ≈ b.surface_balance atol = 1e-8u"W"
end

@testset "A hotter neighbour warms the surface, a colder one cools it" begin
    beside(neighbours) = solve_surface(part, core_temperature; start..., neighbours)
    @test beside(((; fraction = 0.0, temperature = 400.0u"K"),)).insulation_temperature ≈ result.insulation_temperature
    hot = beside(((; fraction = 0.4, temperature = core_temperature),))
    cold = beside(((; fraction = 0.4, temperature = env_vars.sky_temperature),))
    @test hot.success && cold.success
    @test hot.insulation_temperature > result.insulation_temperature > cold.insulation_temperature
    @test hot.net_metabolic < result.net_metabolic < cold.net_metabolic
end

@testset "Neighbours of bare skin" begin
    bare_fibres = FibreProperties(; diameter = 3e-5u"m", length = 0.0u"m", density = 0.0u"1/m^2",
        depth = 0.0u"m", reflectance = 0.2, conductivity = 0.209u"W/m/K")
    bare = Part(;
        body = Body(Sphere(mass = 1.0u"kg", density = 1000.0u"kg/m^3"), Naked()),
        insulation = InsulationParameters(; dorsal = bare_fibres, ventral = bare_fibres,
            depth_compressed = 0.0u"m", longwave_depth_fraction = 1.0),
        traits, environment, ventral_fraction = 0.5,
    )
    beside(neighbours) = solve_surface(bare, core_temperature; start..., neighbours)
    base = beside(())
    @test base.success
    @test beside(((; fraction = 0.0, temperature = 400.0u"K"),)).skin_temperature ≈ base.skin_temperature
    hot = beside(((; fraction = 0.4, temperature = core_temperature),))
    cold = beside(((; fraction = 0.4, temperature = env_vars.sky_temperature),))
    @test hot.success && cold.success
    @test hot.skin_temperature > base.skin_temperature > cold.skin_temperature
end
