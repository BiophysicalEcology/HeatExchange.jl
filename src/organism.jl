# We have a generic abstract type for all organisms
abstract type AbstractOrganism end

# With some generic methods to get the params and body
body(o::AbstractOrganism) = o.body # gets the body from an object of type AbstractOrganism
shape(o::AbstractOrganism) = body(o).shape # gets the shape from an object of type AbstractOrganism
shape(b::AbstractBody) = b.shape # gets the shape from an object of type AbstractBody
insulation(o::AbstractOrganism) = body(o).insulation # gets the insulation from an object of type AbstractOrganism
traits(o::AbstractOrganism) = o.traits # gets the traits from an object of type AbstractOrganism

Base.@kwdef struct FunctionalTraits{F,K,M,B}
    α_org_dorsal::F = Param(0.85, bounds=(0.2, 1.0))
    α_org_ventral::F = Param(0.85, bounds=(0.2, 1.0))
    ϵ_org_dorsal::F = Param(0.95, bounds=(0.1, 1.0))
    ϵ_org_ventral::F = Param(0.95, bounds=(0.1, 1.0))
    F_sky::F = Param(0.4, bounds=(0.3, 0.5))
    F_sub::F = Param(0.4, bounds=(0.3, 0.5))
    p_eyes::F = Param(0.0, bounds=(0.0, 4e-4))
    fO2_extract::F = Param(0.20, bounds=(0.10, 0.30))
    k_body::K = Param(0.5, bounds=(0.412, 2.8), units=u"W/m/K")
    rq::F = Param(0.8, bounds=(0.7, 0.9))
    M1::M = Param(0.013, bounds=(0.01, 0.02))
    M2::M = Param(0.8, bounds=(0.7, 0.9))
    M3::M = Param(0.038, bounds=(0.02, 0.04))
    p_wet::F = Param(0.001, bounds=(0.0, 1.0))
    p_cond::F = Param(0.1, bounds=(0.0, 1.0))
    pant::B = Param(1.0, bounds=(1.0, 10.0))
end

Base.@kwdef mutable struct OrganismalVars{T,P}
    T_core::T = K(20°C)
    T_surf::T = K(20°C)
    T_lung::T = K(20°C)
    ψ_org::P = -707J/kg
end

# Then define a concrete organism struct  
struct Organism{B<:Body,T<:FunctionalTraits} <: AbstractOrganism
    body::B
    traits::T
end