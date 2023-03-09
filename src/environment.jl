
Base.@kwdef struct EnvironmentalParams{A,P,E,L,F}
    α_sub::A = Param(0.2, bounds=(0.0, 1.0))
    ϵ_sub::A = Param(1.0, bounds=(0.2, 1.0))
    ϵ_sky::A = Param(1.0, bounds=(0.2, 1.0))
    P_atmos::P = Param(101325, units=u"Pa")
    elev::E = Param(0, units=u"m")
    fluid::L = Param(0)
    fN2::F = Param(0.79)
    fO2::F = Param(0.2095)
    fCO2::F = Param(0.00042)
end


Base.@kwdef struct EnvironmentalVars{T,R,V,Z,K,Q}
    Ta::T = K(20°C)
    Tsky::T = K(-5°C)
    Tsub::T = K(30°C)
    rh::R = 5
    vel::V = 1m/s
    zen::Z = 20°
    k_sub::K = 0.5W/m/K
    Q_dir::Q = 1000W/m^2
    Q_dif::Q = 200W/m^2
end
