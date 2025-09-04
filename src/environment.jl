abstract type AbstractEnvironmentalPars end

"""
    EnvironmentalParamsMorphoPars

    EnvironmentalParams(α_sub, ϵ_sub, ϵ_sky, elev, fluid, fN2, fO2, fCO2)
    EnvironmentalParams(; kw...)

Environmental parameters for an organism model.
"""
Base.@kwdef struct EnvironmentalPars{A,E,L,F} <: AbstractEnvironmentalPars
    α_sub::A = Param(0.8, bounds=(0.0, 1.0))
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
    T_air::T = (273.15+20.0)K
    T_sky::T = (273.15-5.0)K
    T_sub::T = (273.15+30.0)K
    rh::R = 5.0
    vel::V = 1.0m/s
    P_atmos::P = 101325.0Pa
    zen::Z = 20.0°
    k_sub::K = 0.5W/m/K
    Q_sol::Q = 1000.0W/m^2
    Q_dir::Q = 964.177772475912W/m^2
    Q_dif::Q = 100.0W/m^2
end

Base.@kwdef struct EnvironmentalVarsVec{T,R,V,P,Z,K,Q} <: AbstractEnvironmentalVars
    T_air::Vector{T} = (collect(15.0:5.0:35.0).+273.15).*1.0K
    T_sky::Vector{T} = fill((273.15-5.0)K, length(T_air))
    T_sub::Vector{T} = fill((273.15+30.0)K, length(T_air))
    rh::Vector{R} = fill(5.0, length(T_air))
    vel::Vector{V} = fill(1.0m/s, length(T_air))
    P_atmos::Vector{P} = fill(101325.0Pa, length(T_air))
    zen::Vector{Z} = fill(20.0°, length(T_air))
    k_sub::Vector{K} = fill(0.5W/m/K, length(T_air))
    Q_sol::Vector{Q} = fill(1000.0W/m^2, length(T_air))
    Q_dir::Vector{Q} = fill(964.177772475912W/m^2, length(T_air))
    Q_dif::Vector{Q} = fill(100.0W/m^2, length(T_air))
end
