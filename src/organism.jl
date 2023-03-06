

Base.@kwdef struct OrganismalPars{F}
    α_org_dorsal::F = Param(0.8)
    α_org_ventral::F = Param(0.8)
    ϵ_org::F = Param(0.95)
    F_sky::F = Param(0.4)
    F_sub::F = Param(0.4)
    p_eyes::F = Param(0.03 / 100)
end

Base.@kwdef struct OrganismalVars{T,P,F}
    T_core::T = K(20.0°C)
    T_surf::T = K(20.0°C)
    ψ_org::P = -707J/kg
    p_wet::F = 0.001
    p_cond::F = 0.1
end