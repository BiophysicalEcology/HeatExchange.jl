using HeatExchange
using BiophysicalGeometry
using Unitful, UnitfulMoles
using FluidProperties
using Test

using HeatExchange: solve_part_surface, part_surface_residuals,
    EnvironmentTemperatures, ViewFactors, AtmosphericConditions,
    example_environment_vars, example_environment_pars, example_respiration_pars

# ---------------------------------------------------------------------------
# A single insulated part: an Ellipsoid with a uniform fur layer.
# ---------------------------------------------------------------------------
shape = Ellipsoid(1.0u"kg", 1000.0u"kg/m^3", 2.0, 2.0)
fat   = FatLayer(0.10, 901.0u"kg/m^3")
fur   = FibrousLayer(0.01u"m", 3e-5u"m", 3e7u"1/m^2")
part_body = Body(shape, CompositeInsulation(fur, fat))

# Single-valued part insulation: identical dorsal/ventral fibres so the surface
# solve's side label is immaterial (the part has one insulation).
part_fibres = FibreProperties(;
    diameter = 3e-5u"m", length = 0.01u"m", density = 3e7u"1/m^2",
    depth = 0.01u"m", reflectance = 0.2, conductivity = 0.209u"W/m/K",
)
part_insulation_pars = InsulationParameters(;
    dorsal = part_fibres, ventral = part_fibres,
    depth_compressed = 0.01u"m", longwave_depth_fraction = 1.0,
)

env_vars = example_environment_vars(;
    air_temperature = u"K"(20.0u"°C"),
    wind_speed = 1.0u"m/s",
    global_radiation = 0.0u"W/m^2",
)
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

# Packed per-part environment: view factors + boundaries come from the part's own
# pose/exposure (here a simple upright-ish exposure), independent of insulation.
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

result = solve_part_surface(;
    body = part_body,
    insulation_pars = part_insulation_pars,
    traits,
    environment_vars = packed_environment,
    conduction_fraction = 0.0,
    conductance_coefficient = 0.0u"W/K",
    ventral_fraction = 0.5,
    longwave_depth_fraction = 1.0,
    skin_temperature = core_temperature - 5u"K",
    insulation_temperature = env_vars.air_temperature + 2u"K",
    temperature_tolerance = 1e-3u"K",
)

@testset "solve_part_surface — physical sanity" begin
    @test result.success
    # Skin sits between the air and the core for a heat-generating endotherm
    @test env_vars.air_temperature < result.skin_temperature < core_temperature
    # Insulation surface is cooler than skin (heat flows outward through the fur)
    @test result.insulation_temperature < result.skin_temperature
    # Heat is conducted core → skin
    @test result.net_metabolic > 0.0u"W"
end

@testset "solve_part_surface — flesh conductance extraction is consistent" begin
    # G_flesh · (core − skin) reconstructs net_metabolic exactly
    reconstructed = result.flesh_conductance * (core_temperature - result.skin_temperature)
    @test reconstructed ≈ result.net_metabolic
    @test unit(result.flesh_conductance) == u"W/K"
end

@testset "part_surface_residuals — vanishes at the converged surface solve" begin
    # The non-iterative residual twin must report ≈0 surface and skin-temperature
    # residuals when handed the temperatures `solve_part_surface` converged to.
    setup = (;
        body = part_body,
        insulation_pars = part_insulation_pars,
        traits,
        environment_vars = packed_environment,
        conduction_fraction = 0.0,
        conductance_coefficient = 0.0u"W/K",
        ventral_fraction = 0.5,
        longwave_depth_fraction = 1.0,
        covered_area = 0.0u"m^2",
        characteristic_dim = HeatExchange.characteristic_dimension(
            HeatExchange.VolumeCubeRoot(), part_body),
    )
    res = part_surface_residuals(
        setup, core_temperature, result.skin_temperature, result.insulation_temperature,
        1.0u"W";
        k_flesh = traits.flesh_conductivity,
        pant = 1.0,
        skin_wetness = traits.skin_wetness,
        resp_pars = example_respiration_pars(),
    )
    # Surface balance and skin-temperature residual are driven to zero at the root.
    @test abs(ustrip(u"W", res.surface_balance)) < 1e-2
    @test abs(ustrip(u"K", res.residual_skin_temperature)) < 1e-2
    # Flesh-conducted heat matches solve_part_surface's net_metabolic (both are the
    # same shape-dispatched flesh conduction at the converged temperatures).
    @test res.net_metabolic_heat_internal ≈ result.net_metabolic rtol=1e-3
    # Metabolic/respiration cancel in the surface balance: doubling metabolic input
    # leaves the surface residual unchanged. Compared with an absolute tolerance because both
    # are evaluated at the converged root where surface_balance ≈ 0 (~1e-13 W); a broken
    # cancellation would move it by O(metabolic) = 1 W, far above this atol.
    res2 = part_surface_residuals(
        setup, core_temperature, result.skin_temperature, result.insulation_temperature,
        2.0u"W";
        k_flesh = traits.flesh_conductivity, pant = 1.0,
        skin_wetness = traits.skin_wetness, resp_pars = example_respiration_pars(),
    )
    @test res2.surface_balance ≈ res.surface_balance atol=1e-8u"W"
end

@testset "solve_part_surface — inter-part neighbour exchange" begin
    # The neighbour term is a longwave exchange with a sibling part's surface over a
    # view fraction, carried in `environment_vars.neighbours` as `(; fraction,
    # temperature)` entries. It must vanish for an empty / zero-fraction list and be
    # correctly signed: a hotter neighbour warms this part's surface, a colder one cools it.
    with_neighbours(nbrs) = solve_part_surface(;
        body = part_body,
        insulation_pars = part_insulation_pars,
        traits,
        environment_vars = merge(packed_environment, (; neighbours = nbrs)),
        conduction_fraction = 0.0,
        conductance_coefficient = 0.0u"W/K",
        ventral_fraction = 0.5,
        longwave_depth_fraction = 1.0,
        skin_temperature = core_temperature - 5u"K",
        insulation_temperature = env_vars.air_temperature + 2u"K",
        temperature_tolerance = 1e-3u"K",
    )

    # No neighbours, or a zero-fraction neighbour, reproduces the base solve exactly.
    @test with_neighbours(()).insulation_temperature ≈ result.insulation_temperature
    @test with_neighbours(((; fraction = 0.0, temperature = 400.0u"K"),)).insulation_temperature ≈
          result.insulation_temperature

    # A neighbour hotter than the part's surface pushes heat in → warmer insulation
    # surface and less core→skin flow; a cold neighbour does the opposite.
    hot  = with_neighbours(((; fraction = 0.4, temperature = core_temperature),))
    cold = with_neighbours(((; fraction = 0.4, temperature = env_vars.sky_temperature),))
    @test hot.success && cold.success
    @test hot.insulation_temperature  > result.insulation_temperature
    @test cold.insulation_temperature < result.insulation_temperature
    @test hot.net_metabolic < result.net_metabolic     # less internal heat needed to shed
    @test cold.net_metabolic > result.net_metabolic
end

@testset "solve_part_surface — hotter core drives more heat out" begin
    hotter = solve_part_surface(;
        body = part_body,
        insulation_pars = part_insulation_pars,
        traits = merge(traits, (; core_temperature = core_temperature + 3u"K")),
        environment_vars = packed_environment,
        conduction_fraction = 0.0,
        conductance_coefficient = 0.0u"W/K",
        ventral_fraction = 0.5,
        longwave_depth_fraction = 1.0,
        skin_temperature = core_temperature - 5u"K",
        insulation_temperature = env_vars.air_temperature + 2u"K",
        temperature_tolerance = 1e-3u"K",
    )
    @test hotter.net_metabolic > result.net_metabolic
    @test hotter.skin_temperature > result.skin_temperature
end
