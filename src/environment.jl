abstract type AbstractEnvironmentalPars end

"""
    EnvironmentalParamsMorphoPars

    EnvironmentalParams(α_sub, ϵ_sub, ϵ_sky, elev, fluid, fN2, fO2, fCO2)
    EnvironmentalParams(; kw...)

Environmental parameters for an organism model.
"""
Base.@kwdef struct EnvironmentalPars{A,E,L,F} <: AbstractEnvironmentalPars
    α_sub::A = Param(0.2, bounds=(0.0, 1.0))
    ϵ_sub::A = Param(1.0, bounds=(0.2, 1.0))
    ϵ_sky::A = Param(1.0, bounds=(0.2, 1.0))
    elev::E = Param(0.0, units=u"m")
    fluid::L = Param(0)
    fN2::F = Param(0.79)
    fO2::F = Param(0.2095)
    fCO2::F = Param(0.0003)
end

abstract type AbstractEnvironmentalVars end

"""
    EnvironmentalVars <: AbstractEnvironmentalVars

    EnvironmentalVars(T_air, T_sky, T_sub, rh, vel, P_atmos, zen, k_sub, Q_sol, Q_dir, Q_dif)
    EnvironmentalVars(; kw...)

Environmental variables for an organism model.
"""
Base.@kwdef struct EnvironmentalVars{T,R,V,P,Z,K,Q} <: AbstractEnvironmentalVars
    T_air::T = K(20.0°C)
    T_sky::T = K(-5.0°C)
    T_sub::T = K(30.0°C)
    rh::R = 5.0
    vel::V = 1.0m/s
    P_atmos::P = 101325.0Pa
    zen::Z = 20.0°
    k_sub::K = 0.5W/m/K
    Q_sol::Q = 1000.0W/m^2
    Q_dir::Q = 964.177772475912W/m^2
    Q_dif::Q = 100.0W/m^2
end
