using HeatExchange
using Test
using Unitful

using HeatExchange: HeatCoupling, SharedCore, ConductiveCoupling,
    CompartmentGraph, compartment_graph, num_compartments,
    compartment_part_names, compartment_of, parts_in_compartment,
    contribution_to_conductance, contribution_to_heat_load, build_conductance_matrix

@testset "HeatCoupling constructors" begin
    @test SharedCore() isa HeatCoupling
    # Derived-mode conductive coupling carries a Nothing type parameter
    @test ConductiveCoupling() isa ConductiveCoupling{Nothing}
    @test ConductiveCoupling().interface_conductivity === nothing
    # Explicit-override mode lifts the value type into the parameter
    k = 0.5u"W/m/K"
    @test ConductiveCoupling(k) isa ConductiveCoupling{typeof(k)}
    @test ConductiveCoupling(k).interface_conductivity == k
end

@testset "Single-part degeneracy" begin
    graph = compartment_graph((:body,), ())
    @test num_compartments(graph) == 1
    @test compartment_part_names(graph) == (:body,)
    @test compartment_of(graph, :body) == 1
    @test parts_in_compartment(graph, 1) == (:body,)
    @test graph isa CompartmentGraph{(:body,),1,1}
end

@testset "No couplings — every part its own compartment" begin
    graph = compartment_graph((:head, :torso, :leg), ())
    @test num_compartments(graph) == 3
    @test compartment_of(graph, :head) == 1
    @test compartment_of(graph, :torso) == 2
    @test compartment_of(graph, :leg) == 3
    @test parts_in_compartment(graph, 2) == (:torso,)
end

@testset "Dorsal/ventral SharedCore contracts to one compartment" begin
    graph = compartment_graph((:dorsal, :ventral), ((:dorsal, :ventral),))
    @test num_compartments(graph) == 1
    @test compartment_of(graph, :dorsal) == compartment_of(graph, :ventral)
    @test Set(parts_in_compartment(graph, 1)) == Set((:dorsal, :ventral))
end

@testset "All-SharedCore composite → single compartment" begin
    names = (:a, :b, :c, :d)
    pairs = ((:a, :b), (:b, :c), (:c, :d))
    graph = compartment_graph(names, pairs)
    @test num_compartments(graph) == 1
    @test all(compartment_of(graph, n) == 1 for n in names)
end

@testset "Mixed topology — dog torso SharedCore, limbs separate" begin
    # torso halves share a core; head + four legs are independent compartments
    names = (:dorsal, :ventral, :head, :leg_fl, :leg_fr, :leg_bl, :leg_br)
    pairs = ((:dorsal, :ventral),)   # only the torso halves are SharedCore
    graph = compartment_graph(names, pairs)
    @test num_compartments(graph) == 6
    @test compartment_of(graph, :dorsal) == compartment_of(graph, :ventral)
    # each of the other five parts is its own compartment
    others = (:head, :leg_fl, :leg_fr, :leg_bl, :leg_br)
    @test length(unique(compartment_of(graph, n) for n in others)) == 5
    @test Set(parts_in_compartment(graph, compartment_of(graph, :dorsal))) == Set((:dorsal, :ventral))
end

@testset "Transitive SharedCore groups (two disjoint pools)" begin
    names = (:a, :b, :c, :x, :y)
    pairs = ((:a, :b), (:b, :c), (:x, :y))   # {a,b,c} and {x,y}
    graph = compartment_graph(names, pairs)
    @test num_compartments(graph) == 2
    @test compartment_of(graph, :a) == compartment_of(graph, :c)
    @test compartment_of(graph, :x) == compartment_of(graph, :y)
    @test compartment_of(graph, :a) != compartment_of(graph, :x)
end

@testset "Type stability of graph accessors" begin
    graph = compartment_graph((:dorsal, :ventral, :head), ((:dorsal, :ventral),))
    @test @inferred(num_compartments(graph)) == 2
    @test @inferred(compartment_of(graph, :head)) == 2
    @test @inferred(compartment_part_names(graph)) == (:dorsal, :ventral, :head)
end

@testset "Unknown part errors" begin
    graph = compartment_graph((:body,), ())
    @test_throws ArgumentError compartment_of(graph, :nonexistent)
end

@testset "contribution_to_conductance" begin
    area = 2.0u"cm^2"
    dist = 1.0u"cm"
    k = 0.5u"W/m/K"

    # SharedCore edges are contracted — no matrix entry
    @test contribution_to_conductance(SharedCore(), area, dist, dist, k, k) === nothing

    # Derived series resistance: symmetric parts → G = 1 / (2·d/(k·A))
    G = contribution_to_conductance(ConductiveCoupling(), area, dist, dist, k, k)
    expected = 1 / (2 * dist / (k * area))
    @test G ≈ expected
    @test unit(G) == u"W/K"
    # Series resistance: two equal resistors → half the conductance of one alone
    G_one = 1 / (dist / (k * area))
    @test G ≈ G_one / 2

    # Asymmetric distances/conductivities still add in series
    G2 = contribution_to_conductance(ConductiveCoupling(), area, 1.0u"cm", 3.0u"cm", 0.5u"W/m/K", 0.6u"W/m/K")
    r2 = 1.0u"cm"/(0.5u"W/m/K"*area) + 3.0u"cm"/(0.6u"W/m/K"*area)
    @test G2 ≈ 1 / r2

    # Explicit override: interface conductance coefficient × area
    h = 10.0u"W/m^2/K"
    Goverride = contribution_to_conductance(ConductiveCoupling(h), area, dist, dist, k, k)
    @test Goverride ≈ h * area
    @test unit(Goverride) == u"W/K"
end

@testset "contribution_to_heat_load default is zero" begin
    @test contribution_to_heat_load(SharedCore()) == 0.0u"W"
    @test contribution_to_heat_load(ConductiveCoupling()) == 0.0u"W"
end

@testset "build_conductance_matrix — two compartments, one conductive join" begin
    # head + torso, no SharedCore → 2 compartments; one conductive join between them
    graph = compartment_graph((:head, :torso), ())
    G = 0.4u"W/K"
    K = build_conductance_matrix(graph, ((1, 2, G),))
    @test size(K) == (2, 2)
    @test K[1, 1] ≈ G
    @test K[2, 2] ≈ G
    @test K[1, 2] ≈ -G
    @test K[2, 1] ≈ -G
    # Graph-Laplacian: rows sum to zero
    @test K[1, 1] + K[1, 2] ≈ 0.0u"W/K"
    @test K[2, 1] + K[2, 2] ≈ 0.0u"W/K"
end

@testset "build_conductance_matrix — three compartments, chain" begin
    graph = compartment_graph((:a, :b, :c), ())
    Gab = 0.3u"W/K"
    Gbc = 0.5u"W/K"
    K = build_conductance_matrix(graph, ((1, 2, Gab), (2, 3, Gbc)))
    @test K[1, 1] ≈ Gab
    @test K[2, 2] ≈ Gab + Gbc
    @test K[3, 3] ≈ Gbc
    @test K[1, 2] ≈ -Gab
    @test K[2, 3] ≈ -Gbc
    @test K[1, 3] ≈ 0.0u"W/K"
    # every row of a Laplacian sums to zero
    for row in 1:3
        @test sum(K[row, :]) ≈ 0.0u"W/K"
    end
end

@testset "build_conductance_matrix — empty entries (no conductive joins)" begin
    # dorsal/ventral SharedCore → single compartment, zero conductive network
    graph = compartment_graph((:dorsal, :ventral), ((:dorsal, :ventral),))
    K = build_conductance_matrix(graph, ())
    @test size(K) == (1, 1)
    @test K[1, 1] == 0.0u"W/K"
end
