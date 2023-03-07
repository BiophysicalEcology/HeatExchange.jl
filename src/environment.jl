
Base.@kwdef struct EnvironmentalPars{F,P,E,L}
    α_sub::F = Param(0.85, bounds=(0.2, 1.0))
    ϵ_sub::F = Param(1.0, bounds=(0.2, 1.0))
    ϵ_sky::F = Param(1.0, bounds=(0.2, 1.0))
    P_atmos::P = Param(101325, units=u"Pa")
    elev::E = Param(0, units=u"m")
    fluid::L = Param(0)
    fN2::F = Param(0.7902)
    fO2::F = Param(0.2095)
    fCO2::F = Param(0.0003)
end


Base.@kwdef struct EnvironmentalVars{T,R,V,Z,K,Q}
    Ta::T = K(20.0°C)
    Tsky::T = K(10.0°C)
    Tsub::T = K(30.0°C)
    rh::R = 50
    vel::V = 1m/s
    zen::Z = 0°
    k_sub::K = 2.5W/m/K
    Q_dir::Q = 700W/m^2
    Q_dif::Q = 150W/m^2
end
