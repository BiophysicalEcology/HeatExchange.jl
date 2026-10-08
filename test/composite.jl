using HeatExchange
using BiophysicalGeometry
using Unitful, UnitfulMoles
using FluidProperties
using Test

using HeatExchange: EnvironmentTemperatures, ViewFactors, AtmosphericConditions,
    example_environment_vars, example_environment_pars, example_heat_exchange_traits,
    example_metabolism_pars, example_metabolic_rate_options

# Two half-cylinders joined flat to flat, as a torso's dorsal and ventral halves.
fat = FatLayer(0.10, 901.0u"kg/m^3")
fur = FibrousLayer(0.01u"m", 3e-5u"m", 3e7u"1/m^2")
fibres = FibreProperties(;
    diameter = 3e-5u"m", length = 0.01u"m", density = 3e7u"1/m^2",
    depth = 0.01u"m", reflectance = 0.2, conductivity = 0.209u"W/m/K",
)
insulation = InsulationParameters(; dorsal = fibres, ventral = fibres, depth_compressed = 0.01u"m", longwave_depth_fraction = 1.0)

env_vars = example_environment_vars(; air_temperature = u"K"(20.0u"°C"), wind_speed = 1.0u"m/s")
env_pars = example_environment_pars()
e = (; environment_pars = env_pars, environment_vars = env_vars)
core_temperature = u"K"(37.0u"°C")

setup = (;
    insulation,
    traits = (;
        flesh_conductivity = 0.5u"W/m/K",
        fat_conductivity = 0.2u"W/m/K",
        ϵ_body = 0.99,
        skin_wetness = 0.01,
        insulation_wetness = 0.0,
        bare_skin_fraction = 0.0,
        eye_fraction = 0.0,
    ),
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
    ),
    ventral_fraction = 0.5,
)

half(mass) = Body(HalfCylinder(; mass, density = 1000.0u"kg/m^3", axis_ratio_b = 3.0), CompositeInsulation(fur, fat))
torso = CompositeBody(;
    parts = (; dorsal = half(5.0u"kg"), ventral = half(5.0u"kg")),
    joins = (Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),),
)
halves = parts(torso, (; dorsal = setup, ventral = setup))

traits(respire) = example_heat_exchange_traits(;
    metabolism_pars = example_metabolism_pars(; core_temperature, metabolic_heat_flow = 10.0u"W"),
    options = example_metabolic_rate_options(; respire),
)
solve(parts, respire; kw...) = solve_metabolic_rate(parts, traits(respire), e,
    core_temperature - 5u"K", env_vars.air_temperature + 2u"K"; lung_mass = 10.0u"kg", kw...)

@testset "The joined faces exchange nothing" begin
    @test halves.dorsal.area ≈ BiophysicalGeometry.total_area(torso.parts.dorsal) -
        BiophysicalGeometry.covered_areas(torso.parts, torso.joins).dorsal
end

@testset "solve_metabolic_rate over parts" begin
    result = solve(halves, true)
    @test all(part -> part.success, result.parts)
    @test result.net_metabolic_total ≈ sum(part -> part.net_metabolic, result.parts)
    @test result.parts.dorsal.net_metabolic ≈ result.parts.ventral.net_metabolic
    @test result.metabolic_heat_flow > 0.0u"W"
    @test result.respiration_out !== nothing
    @test env_vars.air_temperature < result.skin_temperature < core_temperature
    @test result.lung_temperature ≈ (core_temperature + result.skin_temperature) * 0.5

    unrespired = solve(halves, false)
    @test unrespired.metabolic_heat_flow == unrespired.net_metabolic_total
    @test unrespired.respiration_out === nothing

    alone = CompositeBody(; parts = (; dorsal = half(10.0u"kg")), joins = ())
    @test solve(parts(alone, (; dorsal = setup)), false).net_metabolic_total > 0.0u"W"
end

@testset "Neighbours from silhouette_factors" begin
    apart = solve(halves, false)
    neighbours = map(view -> view.neighbours, silhouette_factors(torso, Sky(0.5); ndirections = 64, resolution = 32))
    beside = solve(halves, false; neighbours)
    @test all(part -> part.success, beside.parts)
    # The halves of a cylinder hardly see each other.
    @test beside.net_metabolic_total ≈ apart.net_metabolic_total rtol = 0.01
end
