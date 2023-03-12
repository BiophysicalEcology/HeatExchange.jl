
Base.@kwdef struct EnvironmentalParams{A,E,L,F}
    α_sub::A = Param(0.2, bounds=(0.0, 1.0))
    ϵ_sub::A = Param(1.0, bounds=(0.2, 1.0))
    ϵ_sky::A = Param(1.0, bounds=(0.2, 1.0))
    elev::E = Param(0.0, units=u"m")
    fluid::L = Param(0)
    fN2::F = Param(0.79)
    fO2::F = Param(0.2095)
    fCO2::F = Param(0.0003)
end


Base.@kwdef struct EnvironmentalVars{T,R,V,P,Z,K,Q}
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
