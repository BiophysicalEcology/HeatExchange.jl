abstract type AbstractEnvironmentalPars end

"""
    EnvironmentalParamsMorphoPars

    EnvironmentalParams(α_ground, ϵ_ground, ϵ_sky, elevation, fluid, fN2, fO2, fCO2)
    EnvironmentalParams(; kw...)

Environmental parameters for an organism model.
"""
Base.@kwdef struct EnvironmentalPars{A,S,E,L,F} <: AbstractEnvironmentalPars
    α_ground::A = Param(0.2, bounds=(0.0, 1.0))
    ϵ_ground::A = Param(1.0, bounds=(0.0, 1.0))
    ϵ_sky::A = Param(1.0, bounds=(0.0, 1.0))
    elevation::E = Param(0.0, units=u"m")
    fluid::L = Param(0)
    fN2::F = Param(0.7902)
    fO2::F = Param(0.2095)
    fCO2::F = Param(0.0003)
end

abstract type AbstractEnvironmentalVars end

"""
    EnvironmentalVars <: AbstractEnvironmentalVars

    EnvironmentalVars(T_air, T_sky, T_sub, rh, wind_speed, P_atmos, zenith_angle, k_substrate, Q_solar, Q_direct, Q_diffuse)
    EnvironmentalVars(; kw...)

Environmental variables for an organism model.
"""
Base.@kwdef struct EnvironmentalVars{TA,TU,TG,TS,RH,WS,PA,ZA,KS,RS,RB,RD} <: AbstractEnvironmentalVars
    T_air::TA
    T_sky::TU
    T_ground::TG
    T_substrate::TS
    rh::RH
    wind_speed::WS
    P_atmos::PA
    zenith_angle::ZA
    k_substrate::KS
    solar_radiation::RS
    direct_radiation::RB
    diffuse_radiation::RD
end

Base.@kwdef struct EnvironmentalVarsVec{T,R,V,P,Z,K,Q} <: AbstractEnvironmentalVars
    T_air::Matrix{T}
    T_sky::Vector{T}
    T_sub::Matrix{T}
    rh::Matrix{R}
    wind_speed::Matrix{V}
    P_atmos::Vector{P}
    zenith_angle::Vector{Z}
    k_substrate::Matrix{K}
    solar_radiation::Vector{Q}
    direct_radiation::Vector{Q}
    diffuse_radiation::Vector{Q}
end
