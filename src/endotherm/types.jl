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
    BodyRegionValues{T}

Container for average, dorsal, and ventral surface values.

# Fields
- `average::T` — Weighted average value across body surface
- `dorsal::T` — Dorsal (upper/back) surface value
- `ventral::T` — Ventral (lower/belly) surface value
"""
struct BodyRegionValues{T}
    average::T
    dorsal::T
    ventral::T
end

"""
    FibreProperties{D,L,N,I,R,C}

Physical properties of insulation fibres (fur/feathers) for a body region.

# Fields
- `diameter::D` — Fibre diameter
- `length::L` — Fibre length
- `density::N` — Fibre density (count per area)
- `depth::I` — Insulation depth
- `reflectance::R` — Solar reflectance (0-1)
- `conductivity::C` — Thermal conductivity of fibre material (W/m/K)
"""
Base.@kwdef struct FibreProperties{D,L,N,I,R,C}
    diameter::D
    length::L
    density::N
    depth::I
    reflectance::R
    conductivity::C
end

"""
    EnvironmentTemperatures{T1,T2,T3,T4,T5,T6}

Temperatures of environmental surfaces for radiation exchange.

# Fields
- `air::T1` — Air temperature
- `sky::T2` — Effective sky temperature
- `ground::T3` — Ground surface temperature
- `vegetation::T4` — Vegetation canopy temperature
- `bush::T5` — Bush/shrub layer temperature
- `substrate::T6` — Substrate temperature (for conduction)
"""
struct EnvironmentTemperatures{T1,T2,T3,T4,T5,T6}
    air::T1
    sky::T2
    ground::T3
    vegetation::T4
    bush::T5
    substrate::T6
end
function EnvironmentTemperatures(e_vars::AbstractEnvironmentalVars)
    EnvironmentTemperatures(
        e_vars.T_air,
        e_vars.T_sky,
        e_vars.T_ground,
        e_vars.T_vegetation,
        e_vars.T_bush,
        e_vars.T_substrate,
    )
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
- `sky::T1` — View factor to sky
- `ground::T2` — View factor to ground
- `bush::T3` — View factor to bush/shrub layer
- `vegetation::T4` — View factor to vegetation canopy
"""
struct ViewFactors{T1,T2,T3,T4}
    sky::T1
    ground::T2
    bush::T3
    vegetation::T4
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
function AtmosphericConditions(e_vars::AbstractEnvironmentalVars)
    AtmosphericConditions(e_vars.rh, e_vars.wind_speed, e_vars.P_atmos)
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
    MolarFlows

Molar flow rates for respiratory gas exchange (mol/s).

# Fields
- `J_air_in` — Molar flow rate of air inhaled
- `J_air_out` — Molar flow rate of air exhaled
- `J_H2O_in` — Molar flow rate of water vapor inhaled
- `J_H2O_out` — Molar flow rate of water vapor exhaled
- `J_O2_in` — Molar flow rate of oxygen inhaled
- `J_O2_out` — Molar flow rate of oxygen exhaled
- `J_CO2_in` — Molar flow rate of CO2 inhaled
- `J_CO2_out` — Molar flow rate of CO2 exhaled
- `J_N2_in` — Molar flow rate of nitrogen inhaled
- `J_N2_out` — Molar flow rate of nitrogen exhaled
"""
struct MolarFlows{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
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
    HeatFlows{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}

Heat flow components from the heat balance solution.

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
struct HeatFlows{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}
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

"""
    SolarConditions

Solar radiation conditions for heat exchange calculations.

# Fields
- `zenith_angle` — Solar zenith angle
- `global_radiation` — Total incoming solar radiation (W/m²)
- `diffuse_fraction` — Fraction of radiation that is diffuse (0-1)
- `shade` — Fraction of organism in shade (0-1)
"""
Base.@kwdef struct SolarConditions{Z,G,D,S}
    zenith_angle::Z
    global_radiation::G
    diffuse_fraction::D
    shade::S
end
function SolarConditions(e_vars::AbstractEnvironmentalVars)
    SolarConditions(;
        zenith_angle=e_vars.zenith_angle,
        global_radiation=e_vars.global_radiation,
        diffuse_fraction=e_vars.diffuse_fraction,
        shade=e_vars.shade,
    )
end

"""
    TransferCoefficients

Heat and mass transfer coefficients from convection calculations.

# Fields
- `heat` — Heat transfer coefficient (W/m²/K)
- `mass` — Mass transfer coefficient, combined free + forced (m/s)
- `mass_free` — Mass transfer coefficient, free convection only (m/s)
"""
Base.@kwdef struct TransferCoefficients{HC,HD,HDF}
    heat::HC
    mass::HD
    mass_free::HDF = 0.0u"m/s"
end

"""
    MetabolicRates

Metabolic heat generation rates for respiration calculations.

# Fields
- `metabolic` — Current metabolic heat rate (W)
- `sum` — Sum of heat flows for balance (W), defaults to metabolic
- `minimum` — Minimum metabolic rate (W), defaults to metabolic
"""
Base.@kwdef struct MetabolicRates{M,S,N}
    metabolic::M
    sum::S = metabolic
    minimum::N = metabolic
end

"""
    Emissivities

Longwave emissivity values for radiation exchange.

# Fields
- `body_dorsal` — Dorsal body surface emissivity (0-1)
- `body_ventral` — Ventral body surface emissivity (0-1)
- `ground` — Ground surface emissivity (0-1)
- `sky` — Effective sky emissivity (0-1)
"""
Base.@kwdef struct Emissivities{BD,BV,G,S}
    body_dorsal::BD
    body_ventral::BV
    ground::G
    sky::S
end
function Emissivities(rad::RadiationParameters, env::AbstractEnvironmentalPars)
    Emissivities(;
        body_dorsal=rad.ϵ_body_dorsal,
        body_ventral=rad.ϵ_body_ventral,
        ground=env.ϵ_ground,
        sky=env.ϵ_sky,
    )
end

"""
    Absorptivities

Shortwave absorptivity values for solar radiation.

# Fields
- `body_dorsal` — Dorsal body surface absorptivity (0-1)
- `body_ventral` — Ventral body surface absorptivity (0-1)
- `ground` — Ground surface absorptivity (0-1)
"""
Base.@kwdef struct Absorptivities{BD,BV,G}
    body_dorsal::BD
    body_ventral::BV
    ground::G
end
function Absorptivities(rad::RadiationParameters, env::AbstractEnvironmentalPars)
    Absorptivities(;
        body_dorsal=rad.α_body_dorsal,
        body_ventral=rad.α_body_ventral,
        ground=env.α_ground,
    )
end

"""
    InsulationProperties{F,C,A,O,IT,CC}

Computed thermal properties of insulation layers, returned by `insulation_properties()`.

# Fields
- `fibres::BodyRegionValues{<:FibreProperties}` — Fibre properties for average/dorsal/ventral.
- `conductivities::BodyRegionValues` — Effective thermal conductivities (W/m/K).
- `absorption_coefficients::BodyRegionValues` — Absorption coefficients (m⁻¹).
- `optical_thickness::BodyRegionValues` — Optical thickness factors (dimensionless).
- `insulation_test` — Bare-skin test parameter (m⁴); zero indicates no insulation.
- `conductivity_compressed` — Conductivity of compressed ventral insulation (W/m/K).
"""
struct InsulationProperties{F<:BodyRegionValues,C<:BodyRegionValues,A<:BodyRegionValues,O<:BodyRegionValues,IT,CC}
    fibres::F
    conductivities::C
    absorption_coefficients::A
    optical_thickness::O
    insulation_test::IT
    conductivity_compressed::CC
end

"""
    GeometryVariables

Geometric and thermal parameters for heat exchange calculations on a body side.

# Fields
- `side` — Body side (`:dorsal` or `:ventral`)
- `substrate_conductance` — Thermal conductance to substrate (W/K), Q_cond = substrate_conductance × ΔT
- `ventral_fraction` — Fraction of body surface that is ventral (0-1)
- `conduction_fraction` — Fraction of surface area in contact with substrate (0-1)
- `longwave_depth_fraction` — Fraction of insulation depth for longwave radiation exchange (0-1)
"""
Base.@kwdef struct GeometryVariables{S,SC,VF,CF,LDF}
    side::S
    substrate_conductance::SC
    ventral_fraction::VF
    conduction_fraction::CF
    longwave_depth_fraction::LDF
end
