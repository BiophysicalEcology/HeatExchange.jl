using HeatExchange
using BiophysicalGeometry
using Test
using Unitful

density = 1000.0u"kg/m^3"
dorsal = Body(HalfCylinder(; mass = 10u"kg", density, axis_ratio_b = 3.0), Naked())
ventral = Body(HalfCylinder(; mass = 10u"kg", density, axis_ratio_b = 3.0), Naked())
head = Body(Sphere(; mass = 2u"kg", density), Naked())
torso_join = Join(dorsal = Attachment(Flat(), FullCover()), ventral = Attachment(Flat(), FullCover()))
head_join = Join(dorsal = Attachment(EndA(0.0u"m", 0.0), Disc(0.03u"m")), head = Attachment(Radial(0.0, 0.0), Disc(0.03u"m")))
animal = CompositeBody(; parts = (; dorsal, ventral, head), joins = (torso_join, head_join))
k = 0.5u"W/m/K"
flesh_conductivities = (; dorsal = k, ventral = k, head = k)

@testset "conductance" begin
    area = 2.0u"cm^2"
    d = 1.0u"cm"
    @test conductance(SharedCore(), area, (d, d), (k, k)) == 0.0u"W/K"
    @test conductance(Insulated(), area, (d, d), (k, k)) == 0.0u"W/K"
    @test conductance(ConductiveCoupling(), area, (d, d), (k, k)) ≈ k * area / 2d
    @test conductance(ConductiveCoupling(), area, (1.0u"cm", 3.0u"cm"), (k, 0.6u"W/m/K")) ≈
        1 / (1.0u"cm" / (k * area) + 3.0u"cm" / (0.6u"W/m/K" * area))
    h = 10.0u"W/m^2/K"
    @test conductance(ConductiveCoupling(h), area, (d, d), (k, k)) ≈ h * area
end

@testset "compartments come from the joins" begin
    @test ThermalNetwork(animal, (SharedCore(), Insulated()), flesh_conductivities).compartment == (1, 1, 3)
    @test ThermalNetwork(animal, (Insulated(), SharedCore()), flesh_conductivities).compartment == (1, 2, 1)
    @test ThermalNetwork(animal, (SharedCore(), SharedCore()), flesh_conductivities).compartment == (1, 1, 1)
    @test ThermalNetwork(animal, (Insulated(), Insulated()), flesh_conductivities).compartment == (1, 2, 3)
    network = ThermalNetwork(animal, (SharedCore(), ConductiveCoupling()), flesh_conductivities)
    distance = internal_distance(head_join, animal)
    @test network.conductances[2] == (1, 3, conductance(ConductiveCoupling(), join_area(head_join, animal),
        (distance.parent, distance.child), (k, k)))
end

generation = (; dorsal = 3.0u"W", ventral = 2.0u"W", head = 1.0u"W")
flesh_conductance = (; dorsal = 1.5u"W/K", ventral = 0.5u"W/K", head = 1.0u"W/K")
skin = (; dorsal = 301.0u"K", ventral = 299.0u"K", head = 300.0u"K")

@testset "separate cores" begin
    network = ThermalNetwork(animal, (Insulated(), Insulated()), flesh_conductivities)
    cores = core_temperatures(network, generation, flesh_conductance, skin)
    @test cores.dorsal ≈ 303.0u"K"
    @test cores.ventral ≈ 303.0u"K"
    @test cores.head ≈ 301.0u"K"
end

@testset "a shared core balances its parts' heat together" begin
    network = ThermalNetwork(animal, (SharedCore(), Insulated()), flesh_conductivities)
    cores = core_temperatures(network, generation, flesh_conductance, skin)
    @test cores.dorsal == cores.ventral ≈ (5.0u"W" + 1.5u"W/K" * 301.0u"K" + 0.5u"W/K" * 299.0u"K") / 2.0u"W/K"
    @test cores.head ≈ 301.0u"K"
end

@testset "conduction between compartments" begin
    network = ThermalNetwork(animal, (SharedCore(), ConductiveCoupling(1000.0u"W/m^2/K")), flesh_conductivities)
    G = network.conductances[2][3]
    cores = core_temperatures(network, generation, flesh_conductance, skin)
    across = G * (cores.dorsal - cores.head)
    @test across > 0.0u"W"
    @test 1.5u"W/K" * (cores.dorsal - skin.dorsal) + 0.5u"W/K" * (cores.ventral - skin.ventral) + across ≈ 5.0u"W"
    @test 1.0u"W/K" * (cores.head - skin.head) ≈ 1.0u"W" + across
end

@testset "a regulated core is held, and the rest float" begin
    network = ThermalNetwork(animal, (SharedCore(), ConductiveCoupling(1000.0u"W/m^2/K")), flesh_conductivities)
    G = network.conductances[2][3]
    cores = core_temperatures(network, generation, flesh_conductance, skin; regulated = (; ventral = 310.0u"K"))
    @test cores.dorsal == cores.ventral == 310.0u"K"
    @test 1.0u"W/K" * (cores.head - skin.head) ≈ 1.0u"W" + G * (cores.dorsal - cores.head)
    @test @inferred(core_temperatures(network, generation, flesh_conductance, skin)) isa NamedTuple
end
