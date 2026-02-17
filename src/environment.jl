abstract type AbstractEnvironmentalPars end

"""
    EnvironmentalPars

    EnvironmentalPars(ground_albedo, ground_emissivity, sky_emissivity, elevation, fluid, gas_fractions, ...)
    EnvironmentalPars(; kw...)

Environmental parameters for an organism model.
"""
Base.@kwdef struct EnvironmentalPars{AG,EG,ES,EL,FL,G,CE,CD} <: AbstractEnvironmentalPars
    ground_albedo::AG = Param(0.2, bounds=(0.0, 1.0))
    ground_emissivity::EG = Param(1.0, bounds=(0.0, 1.0))
    sky_emissivity::ES = Param(1.0, bounds=(0.0, 1.0))
    elevation::EL = Param(0.0, units=u"m")
    fluid::FL = Param(0)
    gas_fractions::G = GasFractions()
    convection_enhancement::CE = Param(1.0)
    conduction_depth::CD = Param(2.5u"cm", bounds=(0.0, 200.0))
end

abstract type AbstractEnvironmentalVars end

"""
    EnvironmentalVars <: AbstractEnvironmentalVars

    EnvironmentalVars(air_temperature, sky_temperature, ground_temperature, substrate_temperature, relative_humidity, wind_speed, atmospheric_pressure, zenith_angle, k_substrate, global_radiation, diffuse_fraction)
    EnvironmentalVars(; kw...)

Environmental variables for an organism model.
"""
Base.@kwdef struct EnvironmentalVars{TA,TR,TU,TD,TS,TB,TV,RH,WS,PA,ZA,KS,GR,FD,SD} <:
                   AbstractEnvironmentalVars
    air_temperature::TA
    reference_air_temperature::TR = air_temperature
    sky_temperature::TU
    ground_temperature::TD
    substrate_temperature::TS
    bush_temperature::TB = air_temperature
    vegetation_temperature::TV = air_temperature
    relative_humidity::RH
    wind_speed::WS
    atmospheric_pressure::PA
    zenith_angle::ZA
    k_substrate::KS
    global_radiation::GR
    diffuse_fraction::FD
    shade::SD
end

Base.@kwdef struct EnvironmentalVarsVec{TA,TR,TU,TD,TS,TB,TV,RH,WS,PA,ZA,KS,GR,FD,SD} <:
                   AbstractEnvironmentalVars
    air_temperature::Vector{TA}
    reference_air_temperature::Vector{TR} = air_temperature
    sky_temperature::Vector{TU}
    ground_temperature::Vector{TD}
    substrate_temperature::Vector{TS}
    bush_temperature::Vector{TB} = air_temperature
    vegetation_temperature::Vector{TV} = air_temperature
    relative_humidity::Vector{RH}
    wind_speed::Vector{WS}
    atmospheric_pressure::Vector{PA}
    zenith_angle::Vector{ZA}
    k_substrate::Vector{KS}
    global_radiation::Vector{GR}
    diffuse_fraction::Vector{FD}
    shade::Vector{SD}
end
