abstract type AbstractEnvironmentalPars end

"""
    EnvironmentalParamsMorphoPars

    EnvironmentalParams(α_ground, ϵ_ground, ϵ_sky, elevation, fluid, fN2, fO2, fCO2)
    EnvironmentalParams(; kw...)

Environmental parameters for an organism model.
"""
Base.@kwdef struct EnvironmentalPars{AG,EG,ES,EL,FL,G,CE,CD} <: AbstractEnvironmentalPars
    α_ground::AG = Param(0.2, bounds=(0.0, 1.0))
    ϵ_ground::EG = Param(1.0, bounds=(0.0, 1.0))
    ϵ_sky::ES = Param(1.0, bounds=(0.0, 1.0))
    elevation::EL = Param(0.0, units=u"m")
    fluid::FL = Param(0)
    gas::G = GasFractions()
    convection_enhancement::CE = Param(1.0)
    conduction_depth::CD = Param(2.5u"cm", bounds=(0.0, 200.0))
end

abstract type AbstractEnvironmentalVars end

"""
    EnvironmentalVars <: AbstractEnvironmentalVars

    EnvironmentalVars(T_air, T_sky, T_sub, rh, wind_speed, P_atmos, zenith_angle, k_substrate, global_solar_radiation, diffuse_fraction)
    EnvironmentalVars(; kw...)

Environmental variables for an organism model.
"""
Base.@kwdef struct EnvironmentalVars{TA,TR,TU,TD,TS,TB,TV,RH,
        WS,PA,ZA,KS,GR,FD,SD} <: AbstractEnvironmentalVars
    T_air::TA
    T_air_reference::TR = T_air
    T_sky::TU
    T_ground::TD
    T_substrate::TS
    T_bush::TB = T_air
    T_vegetation::TV = T_air
    rh::RH
    wind_speed::WS
    P_atmos::PA
    zenith_angle::ZA
    k_substrate::KS
    global_radiation::GR
    diffuse_fraction::FD
    shade::SD
end

Base.@kwdef struct EnvironmentalVarsVec{TA,TR,TU,TD,TS,TB,TV,RH,
        WS,PA,ZA,KS,GR,FD,SD} <: AbstractEnvironmentalVars
    T_air::Vector{TA}
    T_air_reference::Vector{TR} = T_air
    T_sky::Vector{TU}
    T_ground::Vector{TD}
    T_substrate::Vector{TS}
    T_bush::Vector{TB} = T_air
    T_vegetation::Vector{TV} = T_air
    rh::Vector{RH}
    wind_speed::Vector{WS}
    P_atmos::Vector{PA}
    zenith_angle::Vector{ZA}
    k_substrate::Vector{KS}
    global_radiation::Vector{GR}
    diffuse_fraction::Vector{FD}
    shade::Vector{SD}
end
