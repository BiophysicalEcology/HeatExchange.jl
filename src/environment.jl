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
    fluid::FL = Air()
    gas_fractions::G = GasFractions()
    convection_enhancement::CE = Param(1.0)
    conduction_depth::CD = Param(2.5u"cm", bounds=(0.0, 200.0))
end

abstract type AbstractEnvironmentalVars end

"""
    EnvironmentalVars <: AbstractEnvironmentalVars

    EnvironmentalVars(air_temperature, sky_temperature, ground_temperature, substrate_temperature, relative_humidity, wind_speed, atmospheric_pressure, zenith_angle, substrate_conductivity, global_radiation, diffuse_fraction)
    EnvironmentalVars(; kw...)

Environmental variables for an organism model.

`shade` is the share of the sky hidden by the canopy, which is at `canopy_temperature`, and
`understory_cover` the share of the ground hidden by the understory, which is at
`understory_temperature`.
"""
Base.@kwdef struct EnvironmentalVars{TA,TR,TU,TD,TS,TB,TV,RH,WS,PA,ZA,KS,GR,FD,SD,UC} <:
                   AbstractEnvironmentalVars
    air_temperature::TA
    reference_air_temperature::TR = air_temperature
    sky_temperature::TU
    ground_temperature::TD
    substrate_temperature::TS
    understory_temperature::TB = air_temperature
    canopy_temperature::TV = reference_air_temperature
    relative_humidity::RH
    wind_speed::WS
    atmospheric_pressure::PA
    zenith_angle::ZA
    substrate_conductivity::KS
    global_radiation::GR
    diffuse_fraction::FD
    shade::SD
    understory_cover::UC = 0.0
end

Base.@kwdef struct EnvironmentalVarsVec{TA,TR,TU,TD,TS,TB,TV,RH,WS,PA,ZA,KS,GR,FD,SD,UC} <:
                   AbstractEnvironmentalVars
    air_temperature::Vector{TA}
    reference_air_temperature::Vector{TR} = air_temperature
    sky_temperature::Vector{TU}
    ground_temperature::Vector{TD}
    substrate_temperature::Vector{TS}
    understory_temperature::Vector{TB} = air_temperature
    canopy_temperature::Vector{TV} = reference_air_temperature
    relative_humidity::Vector{RH}
    wind_speed::Vector{WS}
    atmospheric_pressure::Vector{PA}
    zenith_angle::Vector{ZA}
    substrate_conductivity::Vector{KS}
    global_radiation::Vector{GR}
    diffuse_fraction::Vector{FD}
    shade::Vector{SD}
    understory_cover::Vector{UC} = zero(shade)
end
