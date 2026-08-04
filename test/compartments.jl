using HeatExchange
using Test
using Unitful

using HeatExchange: HeatCoupling, SharedCore, ConductiveCoupling,
    CompartmentGraph, compartment_graph, num_compartments,
    compartment_part_names, compartment_of, parts_in_compartment

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
