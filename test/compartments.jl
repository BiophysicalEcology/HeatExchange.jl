using HeatExchange
using Test
using Unitful

using HeatExchange: HeatCoupling, SharedCore, ConductiveCoupling,
    CompartmentGraph, compartment_graph, num_compartments,
    compartment_part_names, compartment_of, parts_in_compartment,
    contribution_to_conductance, contribution_to_heat_load, build_conductance_matrix,
    solve_core_temperatures

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

@testset "solve_core_temperatures — single compartment reduces to zbrent criterion" begin
    # K = 1: core = skin + (metabolic − respiration) / G_flesh, and
    # net_metabolic = G_flesh·(core − skin) = metabolic − respiration exactly.
    graph = compartment_graph((:body,), ())
    net_generation = (10.0u"W",)          # metabolic − respiration
    flesh_conductance = (2.0u"W/K",)      # G_flesh
    skin = 305.0u"K"
    flesh_weighted_skin = (flesh_conductance[1] * skin,)
    cores = solve_core_temperatures(graph, (), net_generation, flesh_conductance, flesh_weighted_skin)
    @test length(cores) == 1
    @test cores[1] ≈ 310.0u"K"                       # 305 + 10/2
    # the closure criterion holds to machine precision
    net_metabolic = flesh_conductance[1] * (cores[1] - skin)
    @test net_metabolic ≈ net_generation[1]
end

@testset "solve_core_temperatures — two independent compartments" begin
    # No coupling entries → each compartment solves standalone.
    graph = compartment_graph((:a, :b), ())
    net_generation = (10.0u"W", 4.0u"W")
    flesh_conductance = (2.0u"W/K", 1.0u"W/K")
    skins = (305.0u"K", 300.0u"K")
    flesh_weighted_skin = (flesh_conductance[1]*skins[1], flesh_conductance[2]*skins[2])
    cores = solve_core_temperatures(graph, (), net_generation, flesh_conductance, flesh_weighted_skin)
    @test cores[1] ≈ 310.0u"K"    # 305 + 10/2
    @test cores[2] ≈ 304.0u"K"    # 300 + 4/1
end

@testset "solve_core_temperatures — two conductively coupled compartments" begin
    # Hand-verified: G_flesh = 1 W/K each, skins = 300 K, coupling G12 = 0.5 W/K,
    # generations 5 W and 1 W → cores 304 K and 302 K, with 1 W flowing 1→2.
    graph = compartment_graph((:a, :b), ())
    G12 = 0.5u"W/K"
    net_generation = (5.0u"W", 1.0u"W")
    flesh_conductance = (1.0u"W/K", 1.0u"W/K")
    skin = 300.0u"K"
    flesh_weighted_skin = (flesh_conductance[1]*skin, flesh_conductance[2]*skin)
    cores = solve_core_temperatures(graph, ((1, 2, G12),), net_generation, flesh_conductance, flesh_weighted_skin)
    @test cores[1] ≈ 304.0u"K"
    @test cores[2] ≈ 302.0u"K"
    # Energy conservation per compartment core node
    flesh_loss_1 = flesh_conductance[1] * (cores[1] - skin)
    coupling_flow = G12 * (cores[1] - cores[2])
    @test flesh_loss_1 + coupling_flow ≈ net_generation[1]           # comp 1: 4 + 1 = 5
    flesh_loss_2 = flesh_conductance[2] * (cores[2] - skin)
    @test flesh_loss_2 ≈ net_generation[2] + coupling_flow           # comp 2: 2 = 1 + 1
end

@testset "solve_core_temperatures — all-SharedCore composite is one node" begin
    # Two half-cylinders sharing a core → single compartment, generations sum,
    # flesh conductances sum, skin contributions flesh-weighted.
    graph = compartment_graph((:dorsal, :ventral), ((:dorsal, :ventral),))
    @test num_compartments(graph) == 1
    gen = (3.0u"W", 2.0u"W")
    gflesh = (1.5u"W/K", 0.5u"W/K")
    skins = (301.0u"K", 299.0u"K")
    net_generation = (gen[1] + gen[2],)                              # 5 W
    flesh_conductance = (gflesh[1] + gflesh[2],)                     # 2 W/K
    flesh_weighted_skin = (gflesh[1]*skins[1] + gflesh[2]*skins[2],) # 1.5·301 + 0.5·299
    cores = solve_core_temperatures(graph, (), net_generation, flesh_conductance, flesh_weighted_skin)
    # core = (net_gen + Σ G_flesh·skin) / Σ G_flesh
    expected = (5.0u"W" + (1.5u"W/K"*301.0u"K" + 0.5u"W/K"*299.0u"K")) / 2.0u"W/K"
    @test cores[1] ≈ expected
end

@testset "solve_core_temperatures — type stability" begin
    graph = compartment_graph((:a, :b), ())
    cores = @inferred solve_core_temperatures(graph, ((1, 2, 0.5u"W/K"),),
        (5.0u"W", 1.0u"W"), (1.0u"W/K", 1.0u"W/K"), (300.0u"W", 300.0u"W"))
    @test length(cores) == 2
end
