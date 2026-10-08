# In a symmetric setting without sun or ground contact, a cylinder solved by its dorsal
# and ventral sides gives the same metabolic heat flow as the cylinder as a single part,
# and as two half-cylinders joined at their flat faces.

using HeatExchange
using BiophysicalGeometry
using Unitful, UnitfulMoles
using FluidProperties
using Test

using HeatExchange: EnvironmentTemperatures, ViewFactors, AtmosphericConditions,
    example_environment_vars, example_environment_pars, example_insulation_pars,
    example_conduction_pars_external, example_conduction_pars_internal, example_evaporation_pars,
    example_radiation_pars, example_heat_exchange_traits

mass = 1.0u"kg"
density = 1000.0u"kg/m^3"
insulation = example_insulation_pars()
fur = FibrousLayer(insulation.dorsal.depth, insulation.dorsal.diameter, insulation.dorsal.density)
fat = FatLayer(0.0, 901.0u"kg/m^3")
internal = example_conduction_pars_internal()
evaporation = example_evaporation_pars()
radiation = example_radiation_pars()
env_vars = example_environment_vars()
env_pars = example_environment_pars()
e = (; environment_pars = env_pars, environment_vars = env_vars)

full_body = Body(Cylinder(; mass, density, axis_ratio_b = 3.0), CompositeInsulation(fur, fat))
traits = example_heat_exchange_traits(;
    shape_pars = full_body.shape,
    insulation_pars = insulation,
    conduction_pars_external = example_conduction_pars_external(; conduction_fraction = 0.0),
    conduction_pars_internal = internal,
    radiation_pars = radiation,
    evaporation_pars = evaporation,
)
skin_temperature = metabolism_pars(traits).core_temperature - 5u"K"
insulation_temperature = env_vars.air_temperature + 2u"K"
sides_solution = solve_metabolic_rate(Organism(full_body, traits), e, skin_temperature, insulation_temperature)

setup(view_factors, ϵ_body; kw...) = (;
    insulation,
    traits = (;
        flesh_conductivity = internal.flesh_conductivity,
        fat_conductivity = internal.fat_conductivity,
        ϵ_body,
        skin_wetness = evaporation.skin_wetness,
        insulation_wetness = evaporation.insulation_wetness,
        bare_skin_fraction = evaporation.bare_skin_fraction,
        eye_fraction = evaporation.eye_fraction,
    ),
    environment = (;
        temperature = EnvironmentTemperatures(
            env_vars.air_temperature, env_vars.sky_temperature,
            env_vars.ground_temperature, env_vars.vegetation_temperature,
            env_vars.bush_temperature, env_vars.substrate_temperature,
        ),
        view_factors,
        atmos = AtmosphericConditions(env_vars),
        fluid = env_pars.fluid,
        solar_flow = 0.0u"W",
        gas_fractions = env_pars.gas_fractions,
        convection_enhancement = env_pars.convection_enhancement,
    ),
    ventral_fraction = 0.5,
    kw...,
)
solve(parts) = solve_metabolic_rate(parts, traits, e, skin_temperature, insulation_temperature; lung_mass = mass)

sky = radiation.sky_view_factor
whole = solve(parts(CompositeBody(; parts = (; whole = full_body), joins = ()),
    (; whole = setup(ViewFactors(sky, 1 - sky, 0.0, 0.0), radiation.body_emissivity_dorsal))))

half = Body(HalfCylinder(; mass = mass / 2, density, axis_ratio_b = 3.0), CompositeInsulation(fur, fat))
torso = CompositeBody(;
    parts = (; dorsal = half, ventral = half),
    joins = (Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),),
)
halves_parts = parts(torso, (;
    dorsal = setup(ViewFactors(sky * 2, 0.0, 0.0, 0.0), radiation.body_emissivity_dorsal),
    ventral = setup(ViewFactors(0.0, (1 - sky) * 2, 0.0, 0.0), radiation.body_emissivity_ventral),
))
halves = solve(halves_parts)

@testset "Two halves seeing what the whole sees carry its heat" begin
    open = setup(ViewFactors(sky, 1 - sky, 0.0, 0.0), radiation.body_emissivity_dorsal)
    flat_faces = (Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover())),)
    for (whole_shape, half_shape) in (
        (Cylinder(; mass, density, axis_ratio_b = 3.0), HalfCylinder(; mass = mass / 2, density, axis_ratio_b = 3.0)),
        (Sphere(; mass, density), HalfSphere(; mass = mass / 2, density)),
    )
        body(shape) = Body(shape, CompositeInsulation(fur, fat))
        one = solve(parts(CompositeBody(; parts = (; whole = body(whole_shape)), joins = ()), (; whole = open)))
        two = solve(parts(CompositeBody(; parts = (; dorsal = body(half_shape), ventral = body(half_shape)), joins = flat_faces),
            (; dorsal = open, ventral = open)))
        @test two.net_metabolic_total ≈ one.net_metabolic_total rtol = 1e-6
        @test two.skin_temperature ≈ one.skin_temperature rtol = 1e-6
    end
end

@testset "Halves convect as the whole" begin
    @test halves_parts.dorsal.characteristic_dimension ≈ characteristic_dimension(VolumeCubeRoot(), full_body)
end

@testset "The whole body as one part" begin
    @test whole.metabolic_heat_flow ≈ sides_solution.energy_flows.metabolic_heat_flow rtol = 1e-3
    @test whole.skin_temperature ≈ sides_solution.thermoregulation.skin_temperature rtol = 1e-4
    @test whole.parts.whole.net_metabolic ≈ sides_solution.energy_flows.dorsal.net_metabolic rtol = 1e-4
end

@testset "Two half-cylinders" begin
    @test halves.parts.dorsal.net_metabolic ≈ halves.parts.ventral.net_metabolic
    @test halves.skin_temperature ≈ sides_solution.thermoregulation.skin_temperature rtol = 1e-3
    @test halves.insulation_temperature ≈ sides_solution.thermoregulation.insulation_temperature rtol = 1e-3
    @test halves.metabolic_heat_flow ≈ sides_solution.energy_flows.metabolic_heat_flow rtol = 1e-3
end
