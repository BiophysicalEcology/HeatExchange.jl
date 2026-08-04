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
