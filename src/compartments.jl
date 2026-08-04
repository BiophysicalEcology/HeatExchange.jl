# Heat coupling between parts + compartment partitioning.
#
# This is physics *topology* — how heat moves between joined parts — not geometry
# and not behaviour, so it lives in HeatExchange. `BiophysicalGeometry` stays pure
# geometry; `BiophysicalBehaviour` merely stores a coupling per join and hands it
# here at problem-construction time.
#
# The compartment partition is compile-time known (part names and coupling types
# are all in type parameters), so it is built once at `HeatBalanceProblem`
# construction and frozen into a `CompartmentGraph{...,NumCompartments}` whose
# `NumCompartments` type parameter lets the downstream linear solve size an
# `SMatrix{NumCompartments,NumCompartments}` statically.

"""
    HeatCoupling

Abstract supertype describing *how* heat exchanges between two joined parts.
Stored per join (keyed by join name) in `HeatBalanceProblem.couplings`. A join
with no coupling entry is insulated by default — zero conductance, zero heat load —
so no explicit `InsulatedJoin` sentinel is needed.

New coupling physics (e.g. blood perfusion, radiative enclosure) is added as new
subtypes plus `contribution_to_conductance_matrix` / `contribution_to_heat_load`
methods, never by rewriting the compartment solve.
"""
abstract type HeatCoupling end

"""
    SharedCore() <: HeatCoupling

The two parts equilibrate into a single thermal node: the join is contracted and
the parts share one `core_temperature`. Contributes zero to the conductance matrix
(the edge is already contracted away by the compartment partition) and zero heat
load. The dorsal/ventral split is the canonical example — one blood pool across
two body halves.
"""
struct SharedCore <: HeatCoupling end

"""
    ConductiveCoupling(interface_conductivity = nothing) <: HeatCoupling

Series-resistance conduction across a join.

- `ConductiveCoupling()` constructs `ConductiveCoupling{Nothing}`: the interface
  conductance is *derived* from each part's own `flesh_conductivity` and
  `internal_distance` to the join —
  `total_resistance = distance_parent / (conductivity_parent · area) +
   distance_child / (conductivity_child · area)`.
- `ConductiveCoupling(k)` constructs `ConductiveCoupling{typeof(k)}` and uses the
  explicit override `k`.

The two modes are distinguished at the type-parameter level, so
`contribution_to_conductance_matrix` dispatches without a runtime branch in the
hot loop.
"""
struct ConductiveCoupling{InterfaceConductivity} <: HeatCoupling
    interface_conductivity::InterfaceConductivity
end
ConductiveCoupling() = ConductiveCoupling(nothing)

"""
    CompartmentGraph{PartNames,NumParts,NumCompartments}

Partition of an organism's parts into thermal compartments. Parts joined by
`SharedCore` couplings equilibrate into one compartment (sharing a single
`core_temperature`); every other part is its own compartment.

Built once at `HeatBalanceProblem` construction via union-find over the
`SharedCore` edges, then frozen. `NumCompartments` is a type parameter so the
per-iteration compartment solve can size an `SMatrix{NumCompartments,NumCompartments}`
statically with no heap allocation.

`part_compartment[i]` is the compartment index (`1:NumCompartments`) of the `i`-th
part in `PartNames`. Compartments are labelled in first-appearance order over
`PartNames`.

Degenerate cases fall out for free:
- a single-part `Body` → one compartment → a 1×1 solve → today's numerics;
- an all-`SharedCore` composite → one compartment → one shared core.
"""
struct CompartmentGraph{PartNames,NumParts,NumCompartments}
    part_compartment::NTuple{NumParts,Int}
    function CompartmentGraph{PartNames,NumParts,NumCompartments}(
        part_compartment::NTuple{NumParts,Int},
    ) where {PartNames,NumParts,NumCompartments}
        new{PartNames,NumParts,NumCompartments}(part_compartment)
    end
end

"""
    compartment_graph(part_names::NTuple{N,Symbol}, shared_core_pairs) -> CompartmentGraph

Build the compartment partition from the part names and the list of `SharedCore`
edges (each a `(parent_name, child_name)` pair of part names). Union-find over the
edges groups equilibrating parts; remaining parts are singleton compartments.

This is the core primitive. The `HeatBalanceProblem` constructor derives
`shared_core_pairs` from a body's joins and its `couplings` NamedTuple and calls
this.
"""
function compartment_graph(part_names::NTuple{N,Symbol}, shared_core_pairs) where {N}
    # Union-find with path halving over part indices. Runs once at construction —
    # a small mutable scratch array here is fine; the *result* is a frozen tuple.
    parent = collect(1:N)
    index_of = NamedTuple{part_names}(ntuple(identity, Val(N)))
    function find(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    for (a, b) in shared_core_pairs
        root_a = find(index_of[a])
        root_b = find(index_of[b])
        root_a == root_b || (parent[root_a] = root_b)
    end
    # Relabel roots to 1:NumCompartments in first-appearance order over the parts.
    labels = zeros(Int, N)
    num_compartments = 0
    for i in 1:N
        root = find(i)
        if labels[root] == 0
            num_compartments += 1
            labels[root] = num_compartments
        end
        labels[i] = labels[root]
    end
    part_compartment = ntuple(i -> labels[i], Val(N))
    return CompartmentGraph{part_names,N,num_compartments}(part_compartment)
end

"""
    num_compartments(graph::CompartmentGraph) -> Int

Number of distinct thermal compartments (the dimension of the compartment solve).
"""
num_compartments(::CompartmentGraph{P,N,K}) where {P,N,K} = K

"""
    compartment_part_names(graph::CompartmentGraph) -> NTuple{N,Symbol}

The part names, in the order matching `graph.part_compartment`.
"""
compartment_part_names(::CompartmentGraph{P}) where {P} = P

"""
    compartment_of(graph::CompartmentGraph, part_name::Symbol) -> Int

Compartment index (`1:num_compartments`) that `part_name` belongs to.
"""
function compartment_of(graph::CompartmentGraph{P}, part_name::Symbol) where {P}
    i = findfirst(==(part_name), P)
    i === nothing && throw(ArgumentError("part $part_name is not in this compartment graph"))
    return graph.part_compartment[i]
end

"""
    parts_in_compartment(graph::CompartmentGraph, compartment::Integer) -> Tuple{Vararg{Symbol}}

The part names belonging to `compartment`.
"""
function parts_in_compartment(graph::CompartmentGraph{P}, compartment::Integer) where {P}
    return Tuple(P[i] for i in eachindex(P) if graph.part_compartment[i] == compartment)
end

# ---------------------------------------------------------------------------
# Coupling contributions to the compartment linear system.
#
# Two-method interface (plan §3.4): new coupling physics is a new pair of
# methods, never a rewrite of the compartment solve.
#   contribution_to_conductance(coupling, ...) -> W/K conductance across the join
#                                                 (nothing when the edge is contracted)
#   contribution_to_heat_load(coupling, ...)   -> W right-hand-side shift
#                                                 (zero except for perfusion-like couplings)
# ---------------------------------------------------------------------------

"""
    contribution_to_conductance(coupling, join_area, distance_parent, distance_child,
                                conductivity_parent, conductivity_child) -> conductance | nothing

Thermal conductance (W/K) that a join contributes across the two compartments it
connects, or `nothing` when the join contributes no matrix entry.

- `SharedCore` → `nothing`: the edge is already contracted away by the compartment
  partition, so it never appears in the conductance matrix.
- `ConductiveCoupling{Nothing}` → derived series resistance through each part's
  flesh: `1 / (distance_parent/(conductivity_parent·area) +
  distance_child/(conductivity_child·area))`.
- `ConductiveCoupling` with an explicit value → the override is an interface
  conductance coefficient (W/m²/K); total conductance is `value · join_area`.
"""
contribution_to_conductance(::SharedCore, args...) = nothing

function contribution_to_conductance(::ConductiveCoupling{Nothing}, join_area,
        distance_parent, distance_child, conductivity_parent, conductivity_child)
    resistance = distance_parent / (conductivity_parent * join_area) +
                 distance_child  / (conductivity_child  * join_area)
    # Canonicalise to W/K so every join contributes the same element type into the
    # conductance matrix regardless of the length units its inputs carried.
    return Unitful.uconvert(u"W/K", 1 / resistance)
end

# Explicit override (any non-Nothing parameter): interface conductance coefficient × area.
function contribution_to_conductance(coupling::ConductiveCoupling, join_area,
        distance_parent, distance_child, conductivity_parent, conductivity_child)
    return Unitful.uconvert(u"W/K", coupling.interface_conductivity * join_area)
end

"""
    contribution_to_heat_load(coupling) -> Power

Right-hand-side heat-load shift (W) a coupling adds to the compartment balance.
Zero for the conduction/contraction couplings here; a future
`BloodPerfusionCoupling` (plan §8.4) overrides this with a nonzero advective term.
"""
contribution_to_heat_load(::HeatCoupling) = 0.0u"W"

"""
    build_conductance_matrix(graph, entries) -> SMatrix{K,K}

Assemble the `K×K` compartment conductance (graph-Laplacian) matrix, where
`K == num_compartments(graph)`. `entries` is an iterable of
`(compartment_i, compartment_j, conductance)` triples — one per `ConductiveCoupling`
join whose endpoints lie in *different* compartments (a join internal to a
compartment, e.g. between two `SharedCore`-linked parts, contributes nothing).

Each entry adds `+conductance` to the two diagonal slots and `−conductance` to the
two off-diagonal slots, the standard thermal-network stiffness assembly. The result
is a stack-allocated `SMatrix` so the per-iteration compartment solve
`core_temperatures = matrix \\ heat_load` never touches the heap.
"""
function build_conductance_matrix(::CompartmentGraph{P,N,K}, entries) where {P,N,K}
    conductance_type = _conductance_eltype(entries)
    matrix = zeros(MMatrix{K,K,conductance_type})
    for (i, j, conductance) in entries
        matrix[i, i] += conductance
        matrix[j, j] += conductance
        matrix[i, j] -= conductance
        matrix[j, i] -= conductance
    end
    return SMatrix(matrix)
end

# Element type of the conductance entries (third slot of each triple). Falls back
# to W/K when there are no inter-compartment conductive joins.
function _conductance_eltype(entries)
    for e in entries
        return typeof(e[3])
    end
    return typeof(0.0u"W/K")
end

# ---------------------------------------------------------------------------
# Compartment core-temperature solve.
#
# The heat balance at compartment c's core node, treating its parts' skins and the
# neighbouring compartment cores as boundaries at their current estimates:
#
#   metabolic_c − respiration_c
#       = Σ_{parts p in c} G_flesh_p · (core_c − skin_p)      [conducted to own skins]
#       + Σ_{neighbours j} G_cj       · (core_c − core_j)     [conducted to neighbours]
#
# Collecting the unknown cores on the left gives the linear system
#
#   (L + diag(Σ_p G_flesh_p)) · core = (metabolic − respiration) + Σ_p G_flesh_p · skin_p
#
# where L is the inter-compartment conductance Laplacian (`build_conductance_matrix`).
# Adding the flesh conductance to the diagonal is what makes the system non-singular
# in the degenerate single-compartment case: for K = 1 with one part it reduces to
# core = skin + (metabolic − respiration) / G_flesh, i.e. the exact
# `net_metabolic = G_flesh·(core − skin) = metabolic − respiration` criterion the
# dorsal/ventral root-find enforces today.
# ---------------------------------------------------------------------------

"""
    solve_core_temperatures(graph, conductance_entries, net_generation,
                            flesh_conductance, flesh_weighted_skin) -> SVector{K}

Solve for the per-compartment core temperatures in one linear step.

`graph` fixes the compartment count `K`. `conductance_entries` are the
inter-compartment `(i, j, conductance)` triples (`build_conductance_matrix` input).
The three remaining arguments are per-compartment aggregates over the parts in each
compartment:

- `net_generation[c]` — `Σ (metabolic − respiration)` (W)
- `flesh_conductance[c]` — `Σ G_flesh_p` (W/K), the total core→skin conductance
- `flesh_weighted_skin[c]` — `Σ G_flesh_p · skin_p` (W)

Returns the `SVector{K}` of compartment core temperatures — a stack-allocated,
heap-free solve.
"""
function solve_core_temperatures(graph::CompartmentGraph{P,N,K}, conductance_entries,
        net_generation::NTuple{K}, flesh_conductance::NTuple{K},
        flesh_weighted_skin::NTuple{K}) where {P,N,K}
    laplacian = build_conductance_matrix(graph, conductance_entries)
    # Assemble and solve in unit-stripped SI space (W/K, W → K) so the StaticArrays
    # linear solve sees plain Float64 — Unitful matrix division is fragile.
    system = ustrip.(u"W/K", laplacian) + Diagonal(SVector(ustrip.(u"W/K", flesh_conductance)))
    load = SVector(ustrip.(u"W", net_generation)) .+ SVector(ustrip.(u"W", flesh_weighted_skin))
    return (system \ load) .* u"K"
end
