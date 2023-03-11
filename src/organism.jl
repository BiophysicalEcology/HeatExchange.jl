# We have a generic abstract type for all organisms
abstract type AbstractOrganism end

# With some generic methods to get the params and body
body(o::AbstractOrganism) = o.body # gets the body from an object of type AbstractOrganism
shape(o::AbstractOrganism) = body(o).shape # gets the shape from an object of type AbstractOrganism
shape(b::AbstractBody) = b.shape # gets the shape from an object of type AbstractBody
insulation(o::AbstractOrganism) = body(o).insulation # gets the insulation from an object of type AbstractOrganism
params(o::AbstractOrganism) = o.params # gets the parameters from an object of type AbstractOrganism

Base.@kwdef struct OrganismParams{F,K}
    α_org_dorsal::F = Param(0.8, bounds=(0.2, 1.0))
    α_org_ventral::F = Param(0.8, bounds=(0.2, 1.0))
    ϵ_org_dorsal::F = Param(0.95, bounds=(0.1, 0.0))
    ϵ_org_ventral::F = Param(0.95, bounds=(0.1, 0.0))
    F_sky::F = Param(0.4, bounds=(0.3, 0.5))
    F_sub::F = Param(0.4, bounds=(0.3, 0.5))
    p_eyes::F = Param(3e-4, bounds=(2e-4, 4e-4))
    fO2_extract::F = Param(0.20, bounds=(0.10, 0.30))
    k_body::K = Param(0.5, bounds=(0.412, 2.8))#, units=u"W/m/K")
    rq::F = Param(0.8, bounds=(0.7, 0.9))
    M1::F = Param(0.013, bounds=(0.01, 0.02))
    M2::F = Param(0.8, bounds=(0.7, 0.9))
    M3::F = Param(0.038, bounds=(0.02, 0.04))
end

Base.@kwdef struct OrganismalVars{T,P,F,B}
    T_core::T = K(20°C)
    T_surf::T = K(20°C)
    ψ_org::P = -707J/kg
    p_wet::F = 0.001
    p_cond::F = 0.1
    pant::B = 1
end

# Then define a concrete organism struct  
struct Organism{B<:Body,P<:OrganismParams} <: AbstractOrganism
    body::B
    params::P
end