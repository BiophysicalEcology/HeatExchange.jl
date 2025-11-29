
"""
    BodyPars <: AbstractMorphoParameters

Morphological parameters describing an animal’s body size, shape
and tissue composition.

- `mass` — Body mass (kg).
- `ρ_flesh` — Density of flesh (kg/m³).
- `ρ_fat` — Density of fat (kg/m³).
- `fat_fraction` — Fraction of body mass that is fat (0–1).
- `shape_b` — Current long–short axis ratio of ellipsoid body (>1).
- `shape_b_max` — Maximum allowable `shape_b` value (>1).
- `shape_c` — Plate-like length:height ratio (often equal to `shape_b`).

These parameters are used to construct an object of type Body.
"""
Base.@kwdef struct BodyPars{BM,FL,FA,FF,SB,SC,MB} <: AbstractMorphoParameters
    mass::BM =         Param(65.0u"kg")
    ρ_flesh::FL =      Param(1000.0u"kg/m^3")
    ρ_fat::FA =        Param(901.0u"kg/m^3")
    fat_fraction::FF = Param(0.0, bounds=(0.0, 1.0))
    shape_b::SB =      Param(1.1, bounds=(1.0, Inf))
    shape_c::SC =      Param(1.1, bounds=(1.0, Inf))
    shape_b_max::MB =  Param(5.0, bounds=(1.0, Inf))
end

"""
    IntegumentPars <: AbstractMorphoParameters

Morphological parameters relating to the integument (boundary) of the organism.

# Parameters
- `α_body_dorsal::F` — Shortwave absorptivity of dorsal body surface (0–1).
- `α_body_ventral::F` — Shortwave absorptivity of ventral body surface (0–1).
- `ϵ_body_dorsal::F` — Longwave emissivity of dorsal body surface (0–1).
- `ϵ_body_ventral::F` — Longwave emissivity of ventral body surface (0–1).
- `F_sky::F` — Radiative configuration factor to sky (0–1).
- `F_ground::F` — Radiative configuration factor to substrate/ground (0–1).
- `F_vegetation::F` — Configuration factor to surrounding vegetation (0–1).
- `F_bush::F` — Configuration factor to vegetation at animal height (0–1).
- `eye_fraction::F` — Fraction of surface area that is eye (0–1).
- `skin_wetness::F` — Fraction of skin surface that is wet (0–1).
- `bare_skin_fraction::F` — Fraction of surface available for evaporation that is bare skin (0–1).
- `insulation_fraction::F` — Fraction of surface area covered by insulation (0–1).
- `conduction_fraction::F` — Fraction of ventral surface directly contacting substrate (0–1).
- `ventral_fraction::F` — Fraction of total surface area that is ventral (0–1).

These parameters influence radiative and evaporative exchange.
"""
Base.@kwdef struct IntegumentPars{AD,AV,ED,EV,FS,FG,FV,FB,EF,SW,BF,IF,CF,VF} <: AbstractMorphoParameters
    α_body_dorsal::AD =       Param(0.85, bounds=(0.0, 1.0))
    α_body_ventral::AV =      Param(0.85, bounds=(0.0, 1.0))
    ϵ_body_dorsal::ED =       Param(0.95, bounds=(0.0, 1.0))
    ϵ_body_ventral::EV =      Param(0.95, bounds=(0.0, 1.0))
    F_sky::FS =               Param(0.5, bounds=(0.0, 1.0))
    F_ground::FG =            Param(0.5, bounds=(0.0, 1.0))
    F_vegetation::FV =        Param(0.0, bounds=(0.0, 1.0))
    F_bush::FB =              Param(0.0, bounds=(0.0, 1.0))
    eye_fraction::EF =        Param(0.0, bounds=(0.0, 1.0))
    skin_wetness::SW =        Param(0.0, bounds=(0.0, 1.0))
    bare_skin_fraction::BF =  Param(1.0, bounds=(0.0, 1.0))
    insulation_fraction::IF = Param(0.0, bounds=(0.0, 1.0))
    conduction_fraction::CF = Param(0.0, bounds=(0.0, 1.0))
    ventral_fraction::VF =    Param(0.5, bounds=(0.0, 1.0))
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
Base.@kwdef struct InsulationPars{ICD,ICV,FDD,FDV,FLD,FLV,IDD,IDV,MID,
    MIV,FRD,FRV,IRD,IRV,IDC,FC,LDF,INW} <: AbstractMorphoParameters
    insulation_conductivity_dorsal::ICD     = nothing
    insulation_conductivity_ventral::ICV    = nothing
    fibre_diameter_dorsal::FDD              = Param(30.0u"μm")
    fibre_diameter_ventral::FDV             = Param(30.0u"μm")
    fibre_length_dorsal::FLD                = Param(23.9u"mm")
    fibre_length_ventral::FLV               = Param(23.9u"mm")
    insulation_depth_dorsal::IDD            = Param(2.0u"mm")
    insulation_depth_ventral::IDV           = Param(2.0u"mm")
    max_insulation_depth_dorsal::MID        = Param(2.0u"mm")
    max_insulation_depth_ventral::MIV       = Param(2.0u"mm")
    fibre_density_dorsal::FRD               = Param(3000.0u"cm^-2")
    fibre_density_ventral::FRV              = Param(3000.0u"cm^-2") 
    insulation_reflectance_dorsal::IRD      = Param(0.301, bounds=(0.0, 1.0))
    insulation_reflectance_ventral::IRV     = Param(0.301, bounds=(0.0, 1.0))
    insulation_depth_compressed::IDC        = Param(2.0u"mm")
    fibre_conductivity::FC                  = Param(0.209u"W/m/K")
    longwave_depth_fraction::LDF            = Param(1, bounds=(0.0, 1.0))
    insulation_wetness::INW                 = Param(1, bounds=(0.0, 1.0))
end

"""
    PhysioPars <: AbstractPhysioParameters

A collection of physiological parameters describing metabolic,
respiratory, and thermal tissue properties of an organism.

# Parameters
- `Q_minimum::B` — Minimum allowed metabolic heat generation rate (W), 
    e.g., resting, active.
- `q10::F` — Q10 factor describing metabolic rate sensitivity to core temperature.
- `k_flesh::K` — Thermal conductivity of lean tissue (W/m/K).
- `k_fat::K` — Thermal conductivity of fat tissue (W/m/K).
- `fO2_extract::F` — Fraction of inspired oxygen extracted per breath (0–1).
- `rq::F` — Respiratory quotient relating CO₂ produced to O₂ consumed (0–1).
- `Δ_breath::B` — Temperature offset between ambient air and exhaled air.
- `rh_exit::F` — Relative humidity of exhaled air (fraction 0–1).

All parameters may be given using `Param(...)` wrappers for bounds,
units, or documentation, and support Unitful values.
"""
Base.@kwdef struct PhysioPars{QM,QT,KF,KA,FO,RQ,DB,RE} <: AbstractPhysioParameters
    Q_minimum::QM =              Param(0.0u"W")
    q10::QT =                    Param(2.0)
    k_flesh::KF =                Param(0.9u"W/m/K")
    k_fat::KA =                  Param(0.230u"W/m/K")
    fO2_extract::FO =            Param(0.20, bounds=(0.0, 1.0))
    rq::RQ =                     Param(0.8, bounds=(0.0, 1.2))
    Δ_breath::DB =               Param(0.0u"K")
    rh_exit::RE =                Param(1.0, bounds=(0.0, 1.0))
end

"""
    ThermoregulationPars <: AbstractBehavParameters

Behavioural and physiological parameters governing thermoregulatory
responses such as panting, sweating, core-temperature adjustment, and
changes in body shape or tissue conductivity.

# Fields

- `insulation_step` — Incremental reduction in insulation depth from the
  piloerect (maximum fur depth) state.
- `shape_b_step` — Increment by which the body‐shape parameter `shape_b`
  increases per iteration, allowing the animal to uncurl toward
  `shape_b_max`.
- `T_core_target` — Target (normothermic) core temperature (K).
- `T_core_max` — Maximum core temperature (K).
- `T_core_min` — Minimum core temperature during torpor (K).
- `T_core_step` — Increment by which core temperature is elevated per
  iteration (K).
- `k_flesh_step` — Increment in flesh thermal conductivity (W/m/K).
- `k_flesh_max` — Maximum flesh thermal conductivity (W/m/K).
- `pant` — Multiplier on breathing rate to simulate panting.
- `pant_step` — Increment in panting multiplier per iteration.
- `pant_multiplier` — Multiplier applied to basal metabolic rate at
  maximum panting effort.
- `pant_max` — Maximum panting multiplier.
- `skin_wetness_step` — Increment in surface wetness fraction used to
  model sweating behaviour.
- `skin_wetness_max` — Maximum fraction of body surface area that can be
  wetted (0–1).

All parameters use `Param` wrappers where appropriate for unit support
and bounds checking.
"""
Base.@kwdef struct ThermoregulationPars{
    IS, SBS, TCT, TCMX, TCMN, TCS, KFS, KFM,
    P, PS, PM, PMX, SWS, SWM
} <: AbstractBehavParameters
    insulation_step::IS    = Param(1.0)
    shape_b_step::SBS      = Param(0.1)
    T_core_target::TCT     = Param(37u"°C" |> u"K")
    T_core_max::TCMX       = Param(39u"°C" |> u"K")
    T_core_min::TCMN       = Param(19u"°C" |> u"K")
    T_core_step::TCS       = Param(0.1u"K")
    k_flesh_step::KFS      = Param(0.1u"W/m/K")
    k_flesh_max::KFM       = Param(2.8u"W/m/K")
    pant::P                = Param(1.0)
    pant_step::PS          = Param(0.1)
    pant_multiplier::PM    = Param(1.05)
    pant_max::PMX          = Param(5.0)
    skin_wetness_step::SWS = Param(0.001)
    skin_wetness_max::SWM  = Param(1.0, bounds=(0.0, 1.0))
end
