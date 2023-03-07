

Base.@kwdef struct OrganismalPars{F}
    α_org_dorsal::F = Param(0.8, bounds=(0.2, 1.0))
    α_org_ventral::F = Param(0.8, bounds=(0.2, 1.0))
    ϵ_org::F = Param(0.95, bounds=(0.1, 0.0))
    F_sky::F = Param(0.4, bounds=(0.3, 0.5))
    F_sub::F = Param(0.4, bounds=(0.3, 0.5))
    p_eyes::F = Param(3e-4, bounds=(2e-4, 4e-4))
end

Base.@kwdef struct OrganismalVars{T,P,F}
    T_core::T = K(20.0°C)
    T_surf::T = K(20.0°C)
    ψ_org::P = -707J/kg
    p_wet::F = 0.001
    p_cond::F = 0.1
end