"""
    ConductanceCoeffs{T1,T2,T3}

Thermal conductance coefficients through insulation layers.

Computed by `radiant_temperature()` and used by downstream functions.

# Fields
- `total::T1` — Total composite conductance (compressed + uncompressed pathways)
- `compressed::T2` — Conductance through compressed insulation layer
- `uncompressed::T3` — Conductance through uncompressed insulation layer
"""
struct ConductanceCoeffs{T1,T2,T3}
    total::T1
    compressed::T2
    uncompressed::T3
end

"""
    DivisorCoeffs{T1,T2,T3,T4}

Intermediate divisor values for heat balance calculations.

Computed by `radiant_temperature()` and used by downstream functions.

# Fields
- `geometric::T1` — Geometric/thermal divisor
- `evaporative::T2` — Evaporative heat contribution term
- `numerator::T3` — Temperature solution numerator
- `radiative::T4` — Radiative conductance divisor
"""
struct DivisorCoeffs{T1,T2,T3,T4}
    geometric::T1
    evaporative::T2
    numerator::T3
    radiative::T4
end

"""
    RadiationCoeffs{T1,T2,T3,T4}

Linearized radiation exchange coefficients to environmental surfaces.

# Fields
- `sky::T1` — Radiation coefficient to sky
- `bush::T2` — Radiation coefficient to bush/shrub layer
- `vegetation::T3` — Radiation coefficient to vegetation canopy
- `ground::T4` — Radiation coefficient to ground surface
"""
struct RadiationCoeffs{T1,T2,T3,T4}
    sky::T1
    bush::T2
    vegetation::T3
    ground::T4
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
    DorsalVentral{D,V}

Container for dorsal and ventral surface values without an average.
Used for input parameters where only dorsal/ventral are specified.

# Fields
- `dorsal::D` — Dorsal (upper/back) surface value
- `ventral::V` — Ventral (lower/belly) surface value
"""
struct DorsalVentral{D,V}
    dorsal::D
    ventral::V
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
function EnvironmentTemperatures(environment_vars::AbstractEnvironmentalVars)
    EnvironmentTemperatures(
        environment_vars.air_temperature,
        environment_vars.sky_temperature,
        environment_vars.ground_temperature,
        environment_vars.vegetation_temperature,
        environment_vars.bush_temperature,
        environment_vars.substrate_temperature,
    )
end

"""
    OrganismTemperatures{T1,T2,T3}

Temperatures of organism body layers.

# Fields
- `core_temperature::T1` — Core body temperature
- `skin_temperature::T2` — Skin temperature
- `insulation_temperature::T3` — Insulation/surface temperature
"""
struct OrganismTemperatures{T1,T2,T3}
    core_temperature::T1
    skin_temperature::T2
    insulation_temperature::T3
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
- `relative_humidity::T1` — Relative humidity (fraction 0-1)
- `wind_speed::T2` — Wind speed
- `atmospheric_pressure::T3` — Atmospheric pressure
"""
struct AtmosphericConditions{T1,T2,T3}
    relative_humidity::T1
    wind_speed::T2
    atmospheric_pressure::T3
end
function AtmosphericConditions(environment_vars::AbstractEnvironmentalVars)
    AtmosphericConditions(environment_vars.relative_humidity, environment_vars.wind_speed, environment_vars.atmospheric_pressure)
end

"""
    ThermalConductivities{F,FA,I}

Thermal conductivities of organism tissues and insulation.

# Fields
- `flesh::F` — Thermal conductivity of lean tissue (W/m/K)
- `fat::FA` — Thermal conductivity of fat tissue (W/m/K)
- `insulation::I` — Effective thermal conductivity of insulation (W/m/K), can be `Nothing`
"""
struct ThermalConductivities{F,FA,I}
    flesh::F
    fat::FA
    insulation::I
end

"""
    MolarFluxes

Molar fluxes for respiratory gas exchange (mol/s).

# Fields
- `air` — Molar flux of air
- `water` — Molar flux of water vapor
- `oxygen` — Molar flux of oxygen
- `carbon_dioxide` — Molar flux of CO2
- `nitrogen` — Molar flux of nitrogen
"""
struct MolarFluxes{A,W,O,C,N}
    air::A
    water::W
    oxygen::O
    carbon_dioxide::C
    nitrogen::N
end

"""
    HeatFluxes{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}

Heat flux components from the heat balance solution.

# Fields
- `convection` — Convective heat loss to air
- `conduction` — Conductive heat loss to substrate
- `net_generated` — Net metabolic heat generation
- `skin_evaporation` — Evaporative heat loss from skin
- `insulation_evaporation` — Evaporative heat loss from insulation surface
- `longwave` — Net longwave radiation exchange
- `solar` — Absorbed solar radiation
- `sky_radiation` — Radiation exchange with sky
- `bush_radiation` — Radiation exchange with bush layer
- `vegetation_radiation` — Radiation exchange with vegetation
- `ground_radiation` — Radiation exchange with ground
"""
struct HeatFluxes{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}
    convection::T1
    conduction::T2
    net_generated::T3
    skin_evaporation::T4
    insulation_evaporation::T5
    longwave::T6
    solar::T7
    sky_radiation::T8
    bush_radiation::T9
    vegetation_radiation::T10
    ground_radiation::T11
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
function SolarConditions(environment_vars::AbstractEnvironmentalVars)
    SolarConditions(;
        zenith_angle=environment_vars.zenith_angle,
        global_radiation=environment_vars.global_radiation,
        diffuse_fraction=environment_vars.diffuse_fraction,
        shade=environment_vars.shade,
    )
end

"""
    TransferCoefficients

Transfer coefficients for convection (heat or mass).

# Fields
- `combined` — Combined free + forced transfer coefficient
- `free` — Free convection transfer coefficient
- `forced` — Forced convection transfer coefficient
"""
Base.@kwdef struct TransferCoefficients{C,F,FO}
    combined::C
    free::F
    forced::FO
end

"""
    MetabolicRates

Metabolic heat generation rates for respiration calculations.

# Fields
- `metabolic` — Current metabolic heat rate (W)
- `sum` — Sum of heat fluxes for balance (W), defaults to metabolic
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
- `body` — Body surface emissivities (DorsalVentral, 0-1)
- `ground` — Ground surface emissivity (0-1)
- `sky` — Effective sky emissivity (0-1)
"""
Base.@kwdef struct Emissivities{B,G,S}
    body::B
    ground::G
    sky::S
end
function Emissivities(rad_pars::RadiationParameters, env::AbstractEnvironmentalPars)
    Emissivities(;
        body=DorsalVentral(rad_pars.body_emissivity_dorsal, rad_pars.body_emissivity_ventral),
        ground=env.ground_emissivity,
        sky=env.sky_emissivity,
    )
end

"""
    Absorptivities

Shortwave absorptivity values for solar radiation.

# Fields
- `body` — Body surface absorptivities (DorsalVentral, 0-1)
- `ground` — Ground surface absorptivity (0-1)
"""
Base.@kwdef struct Absorptivities{B,G}
    body::B
    ground::G
end
function Absorptivities(rad_pars::RadiationParameters, env::AbstractEnvironmentalPars)
    Absorptivities(;
        body=DorsalVentral(rad_pars.body_absorptivity_dorsal, rad_pars.body_absorptivity_ventral),
        ground=env.ground_albedo,
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
- `substrate_conductance` — Thermal conductance to substrate (W/K), conduction_flux = substrate_conductance × ΔT
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
