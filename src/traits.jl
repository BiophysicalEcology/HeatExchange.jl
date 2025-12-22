
"""
    ExternalConductionParameters <: AbstractMorphoParameters

Morphological parameters relating to conductive heat exchange with the external environment.

# Parameters
- `conduction_fraction::F` — Fraction of total surface area that is contacting the substrate (0–1).

"""
Base.@kwdef struct ExternalConductionParameters{F} <: AbstractMorphoParameters
    conduction_fraction::F =  Param(0.0, bounds=(0.0, 1.0))
end

"""
    InternalConductionParameters <: AbstractPhysioParameters

Morphological parameters relating to conductive heat flow within the organism.

# Parameters
- `fat_fraction` — Fraction of body mass that is fat (0–1).
- `k_flesh::K` — Thermal conductivity of lean tissue (W/m/K).
- `k_fat::K` — Thermal conductivity of fat tissue (W/m/K).

"""
Base.@kwdef struct InternalConductionParameters{FF,FL,FA,DF} <: AbstractPhysioParameters
    fat_fraction::FF = Param(0.0, bounds=(0.0, 1.0))
    k_flesh::FL      = Param(0.9u"W/m/K")
    k_fat::FA        = Param(0.230u"W/m/K")
    ρ_fat::DF        = Param(901.0u"kg/m^3")
end

"""
    ConvectionParameters <: AbstractMorphoParameters

Morphological parameters relating to convective heat exchange.

# Parameters
- `convection_area` — surface area involved in convection.

"""
Base.@kwdef struct ConvectionParameters{A} <: AbstractMorphoParameters
    convection_area::A = Param(0.0u"m^2")
end

"""
    RadiationParameters <: AbstractRadiationParameters

Morphological parameters relating to radiation exchange.

# Parameters
- `α_body_dorsal::F` — Shortwave absorptivity of dorsal body surface (0–1).
- `α_body_ventral::F` — Shortwave absorptivity of ventral body surface (0–1).
- `ϵ_body_dorsal::F` — Longwave emissivity of dorsal body surface (0–1).
- `ϵ_body_ventral::F` — Longwave emissivity of ventral body surface (0–1).
- `F_sky::F` — Radiative configuration factor to sky (0–1).
- `F_ground::F` — Radiative configuration factor to substrate/ground (0–1).
- `F_vegetation::F` — Configuration factor to surrounding vegetation (0–1).
- `F_bush::F` — Configuration factor to vegetation at animal height (0–1).
- `ventral_fraction::F` — Fraction of total surface area that is ventral (0–1).

"""
Base.@kwdef struct RadiationParameters{AD,AV,ED,EV,SA,DA,VA,FS,FG,FV,FB,
  VF} <: AbstractMorphoParameters
    α_body_dorsal::AD     = Param(0.85, bounds=(0.0, 1.0))
    α_body_ventral::AV    = Param(0.85, bounds=(0.0, 1.0))
    ϵ_body_dorsal::ED     = Param(0.95, bounds=(0.0, 1.0))
    ϵ_body_ventral::EV    = Param(0.95, bounds=(0.0, 1.0))
    A_silhouette::SA      = Param(0.0u"m^2")
    A_dorsal::DA          = Param(0.0u"m^2")
    A_ventral::VA         = Param(0.0u"m^2")
    F_sky::FS             = Param(0.5, bounds=(0.0, 1.0))
    F_ground::FG          = Param(0.5, bounds=(0.0, 1.0))
    F_vegetation::FV      = Param(0.0, bounds=(0.0, 1.0))
    F_bush::FB            = Param(0.0, bounds=(0.0, 1.0))
    ventral_fraction::VF  = Param(0.5, bounds=(0.0, 1.0))
    solar_orientation     = Intermediate() 
end

"""
    EvaporationParameters <: AbstractMorphoParameters

Morphological parameters relating to cutaneous evaporation of the organism.

# Parameters
- `skin_wetness::F` — Fraction of skin surface that is wet (0–1).
- `insulation_wetness::F` — Fraction of insulation surface that is wet (0–1).
- `eye_fraction::F` — Fraction of surface area that is eye (0–1).
- `bare_skin_fraction::F` — Fraction of surface available for evaporation that is bare skin (0–1).
- `insulation_fraction::F` — Fraction of surface area covered by insulation (0–1).

"""
Base.@kwdef struct EvaporationParameters{SW,IW,EF,BF,IF} <: AbstractMorphoParameters
    skin_wetness::SW        = Param(0.0, bounds=(0.0, 1.0))
    insulation_wetness::IW  = Param(1, bounds=(0.0, 1.0))  
    eye_fraction::EF        = Param(0.0, bounds=(0.0, 1.0))
    bare_skin_fraction::BF  = Param(1.0, bounds=(0.0, 1.0))
    insulation_fraction::IF = Param(0.0, bounds=(0.0, 1.0))
end

"""
    HydraulicParameters <: AbstractPhysioParameters

Morphological parameters relating to radiation exchange.

# Parameters
- `water_potential` — Body water potential (ψ_org). Determines humidity at skin surface 
    and liquid water exchange. (J/kg)
- `hydraulic_conductance` — Hydraulic conductance of skin (K_skin) (kg/(m2 s (J/kg))
- `specific_hydration` — Specific hydration of body (m3 / (m3 (J/kg))). Drives liquid water 
  exchange with substrate if K_skin > 0

"""
Base.@kwdef struct HydraulicParameters{WP,HC,SH} <: AbstractPhysioParameters
    water_potential::WP       = Param(0.0u"J/kg", bounds=(-Inf, 0.0))
    hydraulic_conductance::HC = Param(0.0u"kg / (m^2 * s * (J/kg))", bounds=(0.0, Inf))
    specific_hydration::SH    = Param(0.000304u"m^3 / (m^3 * (J/kg))", bounds=(0.0, Inf))
end

"""
    RespirationParameters <: AbstractPhysioParameters

Physiological parameters relating to respiration of the organism.

# Parameters
- `fO2_extract::F` — Fraction of inspired oxygen extracted per breath (0–1).
- `pant` — Multiplier on breathing rate to simulate panting.
- `rq::F` — Respiratory quotient relating CO₂ produced to O₂ consumed (0–1).
- `Δ_breath::B` — Temperature offset between ambient air and exhaled air.
- `rh_exit::F` — Relative humidity of exhaled air (fraction 0–1).

These parameters influence radiative and evaporative exchange.
"""
Base.@kwdef struct RespirationParameters{OE,PA,RQ,DB,RE} <: AbstractPhysioParameters
    fO2_extract::OE = Param(0.2, bounds=(0.0, 1.0))
    pant::PA        = Param(1.0, bounds=(1.0, Inf))
    rq::RQ          = Param(0.8, bounds=(0.5, 1.2))
    Δ_breath::DB    = Param(0.0u"K")
    rh_exit::RE     = Param(1.0, bounds=(0.0, 1.0))
end

"""
    InsulationPars(; kwargs...) <: AbstractMorphoParameters

Parameters describing the physical and optical properties of fur or other
insulative pelage layers on an animal. These values control heat retention,
radiation exchange, and conduction through the insulation.

All parameters may be overridden at construction. Defaults are based on
typical small-mammal pelage properties.

# Parameters

## Thermal Conductivity
- `insulation_conductivity_dorsal` — User-specified dorsal fur conductivity
  (W/m/K). If `nothing`, conductivity is computed from fibre properties.
- `insulation_conductivity_ventral` — Same as above, but for ventral fur.

## Fibre Geometry
- `fibre_diameter_dorsal` — Mean dorsal hair diameter (μm).
- `fibre_diameter_ventral` — Mean ventral hair diameter (μm).
- `fibre_length_dorsal` — Dorsal hair length (mm). Also defines maximum
  achievable dorsal fur depth.
- `fibre_length_ventral` — Ventral hair length (mm).

## Fur Depth
- `insulation_depth_dorsal` — Actual dorsal fur depth (mm).
- `insulation_depth_ventral` — Actual ventral fur depth (mm).
- `max_insulation_depth_dorsal` — Maximum dorsal fur depth (mm), typically equal
  to fibre length.
- `max_insulation_depth_ventral` — Maximum ventral fur depth (mm).

## Fur Density
- `fibre_density_dorsal` — Dorsal hair density (cm⁻²).
- `fibre_density_ventral` — Ventral hair density (cm⁻²).

## Optical/Reflective Properties
- `insulation_reflectance_dorsal` — Dorsal fur reflectance (fraction, 0–1).
- `insulation_reflectance_ventral` — Ventral fur reflectance (fraction, 0–1).

## Compressed Fur (Conduction Layer)
- `insulation_depth_compressed` — Effective depth of compressed fur during
  conduction contact (mm).

## Fibre Material Properties
- `fibre_conductivity` — Thermal conductivity of keratin fibres (W/m/K).

## Radiative Exchange
- `longwave_depth_fraction` — Fraction (0–1) of the fur depth at which
  longwave radiation is exchanged. A value of 1 means radiation interacts at
  the full fur depth; lower values represent radiation being absorbed higher
  in the pelage.

# Notes
- All defaults use `Param` wrappers with units (Unitful.jl).
- Users may specify custom conductivity values to override internally computed
  pelage conductivities.
- This struct represents both dorsal and ventral insulation separately,
  allowing asymmetric fur properties.
"""
Base.@kwdef struct InsulationParameters{ICD,ICV,FDD,FDV,FLD,FLV,IDD,IDV,
    FRD,FRV,IRD,IRV,IDC,FCN,LDF} <: AbstractMorphoParameters
    insulation_conductivity_dorsal::ICD     = nothing
    insulation_conductivity_ventral::ICV    = nothing
    fibre_diameter_dorsal::FDD              = Param(30.0u"μm")
    fibre_diameter_ventral::FDV             = Param(30.0u"μm")
    fibre_length_dorsal::FLD                = Param(23.9u"mm")
    fibre_length_ventral::FLV               = Param(23.9u"mm")
    insulation_depth_dorsal::IDD             = Param(2.0u"mm")
    insulation_depth_ventral::IDV            = Param(2.0u"mm")    
    fibre_density_dorsal::FRD               = Param(3000.0u"cm^-2")
    fibre_density_ventral::FRV              = Param(3000.0u"cm^-2") 
    insulation_reflectance_dorsal::IRD      = Param(0.301, bounds=(0.0, 1.0))
    insulation_reflectance_ventral::IRV     = Param(0.301, bounds=(0.0, 1.0))
    insulation_depth_compressed::IDC        = Param(2.0u"mm")
    fibre_conductivity::FCN                 = Param(0.209u"W/m/K")
    longwave_depth_fraction::LDF            = Param(1, bounds=(0.0, 1.0))
end

"""
    MetabolicParameters <: AbstractPhysioParameters

A collection of physiological parameters relating to metabolic rate.

# Parameters
- `T_core::B` — Core body temperature (K)
- `Q_metabolism::B` — Metabolic heat generation rate (W)
- `q10::F` — Q10 factor describing metabolic rate sensitivity to core temperature.

"""
Base.@kwdef struct MetabolismParameters{TC,QM,QT} <: AbstractPhysioParameters
    T_core::TC       = Param(37u"°C" |> u"K")
    Q_metabolism::QM = Param(0.0u"W")
    q10::QT          = Param(2.0)
    model            = AndrewsPough2()
end