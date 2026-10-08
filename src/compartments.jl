# How heat crosses the joins of a `CompositeBody`, and the core temperatures that follow.

"""
    HeatCoupling

How heat crosses a join. A body's couplings are a tuple with one for each of its
`joins`, in the same order.
"""
abstract type HeatCoupling end

"""
    Insulated() <: HeatCoupling

No heat crosses the join.
"""
struct Insulated <: HeatCoupling end

"""
    SharedCore() <: HeatCoupling

The two parts share one core temperature, as the halves of a torso share one blood pool.
"""
struct SharedCore <: HeatCoupling end

"""
    ConductiveCoupling(interface_conductance = nothing) <: HeatCoupling

Heat conducts across the join. With `nothing` it conducts through the two parts' flesh in
series, from each part's centroid to the join; otherwise `interface_conductance` (W/m²/K)
over the join area sets it.
"""
struct ConductiveCoupling{C} <: HeatCoupling
    interface_conductance::C
end
ConductiveCoupling() = ConductiveCoupling(nothing)

"""
    conductance(coupling, area, distances, conductivities)

The conductance (W/K) across a join of `area`, where `distances` and `conductivities` are
the `(parent, child)` centroid-to-join distances and flesh conductivities.
"""
conductance(::HeatCoupling, area, distances, conductivities) = 0.0u"W/K"
conductance(::ConductiveCoupling{Nothing}, area, distances, conductivities) =
    uconvert(u"W/K", area / (distances[1] / conductivities[1] + distances[2] / conductivities[2]))
conductance(coupling::ConductiveCoupling, area, distances, conductivities) =
    uconvert(u"W/K", coupling.interface_conductance * area)

shares_core(::HeatCoupling) = false
shares_core(::SharedCore) = true

"""
    ThermalNetwork(body::CompositeBody, couplings, flesh_conductivities)

The heat paths between the parts of `body`, given one `HeatCoupling` for each of its joins
and the flesh conductivity of each part, keyed like `body.parts`.

Parts joined by `SharedCore` form a compartment with one core temperature. `compartment`
holds the index of the first part of each part's compartment, and `conductances` a
`(parent, child, conductance)` triple for each join.
"""
struct ThermalNetwork{K,N,C}
    compartment::NTuple{N,Int}
    conductances::C
end
function ThermalNetwork(body::CompositeBody, couplings::Tuple, flesh_conductivities::NamedTuple)
    links = map((j, c) -> _link(body, j, c, flesh_conductivities), body.joins, couplings)
    compartment = _groups(body, map(l -> l.shared, links))
    conductances = map(l -> (l.parent, l.child, l.conductance), links)
    ThermalNetwork{keys(body.parts),length(compartment),typeof(conductances)}(compartment, conductances)
end

_indices(parts::NamedTuple{K}) where {K} = NamedTuple{K}(ntuple(identity, Val(length(K))))

function _link(body, j::Join{P,C}, coupling, flesh_conductivities) where {P,C}
    (; _ends(body, j)..., shared = shares_core(coupling),
       conductance = _conductance(coupling, j, body, (getfield(flesh_conductivities, P), getfield(flesh_conductivities, C))))
end

_ends(body, ::Join{P,C}) where {P,C} = (; parent = getfield(_indices(body.parts), P), child = getfield(_indices(body.parts), C))

# The parts brought together by the joins for which `linked` holds, as the index of the
# first part of each part's group.
_groups(body::CompositeBody, linked::Tuple) =
    foldl(_join_compartments, map((j, shared) -> (; _ends(body, j)..., shared), body.joins, linked);
        init = values(_indices(body.parts)))

_conductance(coupling::HeatCoupling, j, body, conductivities) = conductance(coupling, nothing, nothing, conductivities)
function _conductance(coupling::ConductiveCoupling, j, body, conductivities)
    distance = internal_distance(j, body)
    conductance(coupling, join_area(j, body), (distance.parent, distance.child), conductivities)
end

function _join_compartments(compartment, link)
    link.shared || return compartment
    low, high = minmax(compartment[link.parent], compartment[link.child])
    map(c -> c == high ? low : c, compartment)
end

"""
    core_temperatures(network, generation, flesh_conductance, skin_temperature; regulated = (;))

The core temperature of each part, from the heat each generates net of respiration (W),
its core-to-skin `flesh_conductance` (W/K) and its `skin_temperature`, all keyed like the
body's parts. Each compartment's core balances the heat its parts generate against the heat
they conduct to their own skins and across joins to other compartments.

`regulated`, a single entry like `(; torso = 310.0u"K")`, holds that part's compartment at
the temperature, and its generation is not used.
"""
function core_temperatures(network::ThermalNetwork{K,N}, generation, flesh_conductance, skin_temperature;
        regulated = (;)) where {K,N}
    g = map(x -> ustrip(u"W/K", x), values(flesh_conductance))
    A = zeros(MMatrix{N,N,Float64})
    b = MVector(map((q, c, t) -> ustrip(u"W", q) + c * ustrip(u"K", t), values(generation), g, values(skin_temperature)))
    for i in 1:N
        A[i, i] = g[i]
    end
    for (i, j, conductance) in network.conductances
        G = ustrip(u"W/K", conductance)
        A[i, i] += G
        A[j, j] += G
        A[i, j] -= G
        A[j, i] -= G
    end
    # The parts of a compartment share its first part's core: their balances add there.
    for i in 1:N
        root = network.compartment[i]
        root == i && continue
        for k in 1:N
            A[root, k] += A[i, k]
            A[i, k] = 0.0
        end
        b[root] += b[i]
        A[i, i] = 1.0
        A[i, root] = -1.0
        b[i] = 0.0
    end
    _regulate!(A, b, network, regulated)
    NamedTuple{K}(Tuple(SMatrix(A) \ SVector(b)) .* u"K")
end

_regulate!(A, b, network, ::NamedTuple{()}) = nothing
function _regulate!(A, b, network::ThermalNetwork{K,N}, regulated::NamedTuple{P,<:Tuple{Any}}) where {K,N,P}
    root = network.compartment[getfield(NamedTuple{K}(ntuple(identity, Val(N))), only(P))]
    for k in 1:N
        A[root, k] = 0.0
    end
    A[root, root] = 1.0
    b[root] = ustrip(u"K", only(values(regulated)))
    nothing
end
