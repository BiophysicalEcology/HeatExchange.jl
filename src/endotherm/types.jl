"""
    ConductanceCoeffs{T1,T2,T3}

Thermal conductance coefficients through insulation layers.

Computed by `radiant_temperature()` and used by downstream functions.

# Fields
- `cd1::T1` — Total composite conductance (compressed + uncompressed pathways)
- `cd2::T2` — Conductance through compressed insulation layer
- `cd3::T3` — Conductance through uncompressed insulation layer
"""
struct ConductanceCoeffs{T1,T2,T3}
    cd1::T1
    cd2::T2
    cd3::T3
end

"""
    DivisorCoeffs{T1,T2,T3,T4}

Intermediate divisor values for heat balance calculations.

Computed by `radiant_temperature()` and used by downstream functions.

# Fields
- `dv1::T1` — Geometric/thermal divisor
- `dv2::T2` — Evaporative heat contribution term
- `dv3::T3` — Temperature solution numerator
- `dv4::T4` — Radiative conductance divisor
"""
struct DivisorCoeffs{T1,T2,T3,T4}
    dv1::T1
    dv2::T2
    dv3::T3
    dv4::T4
end

"""
    RadiationCoeffs{T1,T2,T3,T4}

Linearized radiation exchange coefficients to environmental surfaces.

# Fields
- `Q_rad1::T1` — Radiation coefficient to sky
- `Q_rad2::T2` — Radiation coefficient to bush/shrub layer
- `Q_rad3::T3` — Radiation coefficient to vegetation canopy
- `Q_rad4::T4` — Radiation coefficient to ground surface
"""
struct RadiationCoeffs{T1,T2,T3,T4}
    Q_rad1::T1
    Q_rad2::T2
    Q_rad3::T3
    Q_rad4::T4
end

"""
    EnvironmentTemperatures{T1,T2,T3,T4,T5,T6}

Temperatures of environmental surfaces for radiation exchange.

# Fields
- `T_air::T1` — Air temperature
- `T_sky::T2` — Effective sky temperature
- `T_ground::T3` — Ground surface temperature
- `T_vegetation::T4` — Vegetation canopy temperature
- `T_bush::T5` — Bush/shrub layer temperature
- `T_substrate::T6` — Substrate temperature (for conduction)
"""
struct EnvironmentTemperatures{T1,T2,T3,T4,T5,T6}
    T_air::T1
    T_sky::T2
    T_ground::T3
    T_vegetation::T4
    T_bush::T5
    T_substrate::T6
end

"""
    OrganismTemperatures{T1,T2,T3}

Temperatures of organism body layers.

# Fields
- `T_core::T1` — Core body temperature
- `T_skin::T2` — Skin temperature
- `T_insulation::T3` — Insulation/surface temperature
"""
struct OrganismTemperatures{T1,T2,T3}
    T_core::T1
    T_skin::T2
    T_insulation::T3
end

"""
    ViewFactors{T1,T2,T3,T4}

View factors to environmental surfaces for radiation exchange.

# Fields
- `F_sky::T1` — View factor to sky
- `F_ground::T2` — View factor to ground
- `F_bush::T3` — View factor to bush/shrub layer
- `F_vegetation::T4` — View factor to vegetation canopy
"""
struct ViewFactors{T1,T2,T3,T4}
    F_sky::T1
    F_ground::T2
    F_bush::T3
    F_vegetation::T4
end

"""
    AtmosphericConditions{T1,T2,T3}

Atmospheric conditions for heat exchange calculations.

# Fields
- `rh::T1` — Relative humidity (fraction 0-1)
- `wind_speed::T2` — Wind speed
- `P_atmos::T3` — Atmospheric pressure
"""
struct AtmosphericConditions{T1,T2,T3}
    rh::T1
    wind_speed::T2
    P_atmos::T3
end

"""
    ThermalConductivities{F,FA,I}

Thermal conductivities of organism tissues and insulation.

# Fields
- `k_flesh::F` — Thermal conductivity of lean tissue (W/m/K)
- `k_fat::FA` — Thermal conductivity of fat tissue (W/m/K)
- `k_insulation::I` — Effective thermal conductivity of insulation (W/m/K), can be `Nothing`
"""
struct ThermalConductivities{F,FA,I}
    k_flesh::F
    k_fat::FA
    k_insulation::I
end

"""
    MolarFluxes

Molar fluxes for respiratory gas exchange (mol/s).

# Fields
- `J_air_in` — Molar flux of air inhaled
- `J_air_out` — Molar flux of air exhaled
- `J_H2O_in` — Molar flux of water vapor inhaled
- `J_H2O_out` — Molar flux of water vapor exhaled
- `J_O2_in` — Molar flux of oxygen inhaled
- `J_O2_out` — Molar flux of oxygen exhaled
- `J_CO2_in` — Molar flux of CO2 inhaled
- `J_CO2_out` — Molar flux of CO2 exhaled
- `J_N2_in` — Molar flux of nitrogen inhaled
- `J_N2_out` — Molar flux of nitrogen exhaled
"""
struct MolarFluxes{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
    J_air_in::T1
    J_air_out::T2
    J_H2O_in::T3
    J_H2O_out::T4
    J_O2_in::T5
    J_O2_out::T6
    J_CO2_in::T7
    J_CO2_out::T8
    J_N2_in::T9
    J_N2_out::T10
end

"""
    HeatFluxes{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}

Heat flux components from the heat balance solution.

# Fields
- `Q_convection` — Convective heat loss to air
- `Q_conduction` — Conductive heat loss to substrate
- `Q_gen_net` — Net metabolic heat generation
- `Q_evap_skin` — Evaporative heat loss from skin
- `Q_evap_insulation` — Evaporative heat loss from insulation surface
- `Q_longwave` — Net longwave radiation exchange
- `Q_solar` — Absorbed solar radiation
- `Q_rad_sky` — Radiation exchange with sky
- `Q_rad_bush` — Radiation exchange with bush layer
- `Q_rad_vegetation` — Radiation exchange with vegetation
- `Q_rad_ground` — Radiation exchange with ground
"""
struct HeatFluxes{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}
    Q_convection::T1
    Q_conduction::T2
    Q_gen_net::T3
    Q_evap_skin::T4
    Q_evap_insulation::T5
    Q_longwave::T6
    Q_solar::T7
    Q_rad_sky::T8
    Q_rad_bush::T9
    Q_rad_vegetation::T10
    Q_rad_ground::T11
end
