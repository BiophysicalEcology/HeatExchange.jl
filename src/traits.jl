
"""
    ExternalConductionParameters <: AbstractMorphologyParameters

Morphological parameters relating to conductive heat exchange with the external environment.

# Parameters
- `conduction_fraction::F` — Fraction of total surface area that is contacting the substrate (0–1).

"""
Base.@kwdef struct ExternalConductionParameters{F} <: AbstractMorphologyParameters
    conduction_fraction::F = Param(0.0, bounds=(0.0, 1.0))
end

"""
    InternalConductionParameters <: AbstractPhysiologyParameters

Morphological parameters relating to conductive heat flow within the organism.

# Parameters
- `fat_fraction` — Fraction of body mass that is fat (0–1).
- `k_flesh::K` — Thermal conductivity of lean tissue (W/m/K).
- `k_fat::K` — Thermal conductivity of fat tissue (W/m/K).

"""
Base.@kwdef struct InternalConductionParameters{FF,FL,FA,DF} <: AbstractPhysiologyParameters
    fat_fraction::FF = Param(0.0, bounds=(0.0, 1.0))
    k_flesh::FL = Param(0.9u"W/m/K")
    k_fat::FA = Param(0.230u"W/m/K")
    ρ_fat::DF = Param(901.0u"kg/m^3")
end

"""
    ConvectionParameters <: AbstractMorphologyParameters

Morphological parameters relating to convective heat exchange.

# Parameters
- `convection_area` — surface area involved in convection.

"""
Base.@kwdef struct ConvectionParameters{A} <: AbstractMorphologyParameters
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
Base.@kwdef struct RadiationParameters{AD,AV,ED,EV,AS,AT,AC,FS,FG,FV,FB,VF} <:
                   AbstractMorphologyParameters
    α_body_dorsal::AD = Param(0.85, bounds=(0.0, 1.0))
    α_body_ventral::AV = Param(0.85, bounds=(0.0, 1.0))
    ϵ_body_dorsal::ED = Param(0.95, bounds=(0.0, 1.0))
    ϵ_body_ventral::EV = Param(0.95, bounds=(0.0, 1.0))
    A_silhouette::AS = Param(0.0u"m^2")
    A_total::AT = Param(0.0u"m^2")
    A_conduction::AC = Param(0.0u"m^2")
    F_sky::FS = Param(0.5, bounds=(0.0, 1.0))
    F_ground::FG = Param(0.5, bounds=(0.0, 1.0))
    F_vegetation::FV = Param(0.0, bounds=(0.0, 1.0))
    F_bush::FB = Param(0.0, bounds=(0.0, 1.0))
    ventral_fraction::VF = Param(0.5, bounds=(0.0, 1.0))
    solar_orientation = Intermediate()
end

"""
    EvaporationParameters <: AbstractMorphologyParameters

Morphological parameters relating to cutaneous evaporation of the organism.

# Parameters
- `skin_wetness::F` — Fraction of skin surface that is wet (0–1).
- `insulation_wetness::F` — Fraction of insulation surface that is wet (0–1).
- `eye_fraction::F` — Fraction of surface area that is eye (0–1).
- `bare_skin_fraction::F` — Fraction of surface available for evaporation that is bare skin (0–1).
- `insulation_fraction::F` — Fraction of surface area covered by insulation (0–1).

"""
Base.@kwdef struct EvaporationParameters{SW,IW,EF,BF,IF} <: AbstractMorphologyParameters
    skin_wetness::SW = Param(0.0, bounds=(0.0, 1.0))
    insulation_wetness::IW = Param(1, bounds=(0.0, 1.0))
    eye_fraction::EF = Param(0.0, bounds=(0.0, 1.0))
    bare_skin_fraction::BF = Param(1.0, bounds=(0.0, 1.0))
    insulation_fraction::IF = Param(0.0, bounds=(0.0, 1.0))
end

"""
    HydraulicParameters <: AbstractPhysiologyParameters

Morphological parameters relating to radiation exchange.

# Parameters
- `water_potential` — Body water potential (ψ_org). Determines humidity at skin surface 
    and liquid water exchange. (J/kg)
- `hydraulic_conductance` — Hydraulic conductance of skin (K_skin) (kg/(m2 s (J/kg))
- `specific_hydration` — Specific hydration of body (m3 / (m3 (J/kg))). Drives liquid water 
  exchange with substrate if K_skin > 0

"""
Base.@kwdef struct HydraulicParameters{WP,HC,SH} <: AbstractPhysiologyParameters
    water_potential::WP = Param(0.0u"J/kg", bounds=(-Inf, 0.0))
    hydraulic_conductance::HC = Param(0.0u"kg / (m^2 * s * (J/kg))", bounds=(0.0, Inf))
    specific_hydration::SH = Param(0.000304u"m^3 / (m^3 * (J/kg))", bounds=(0.0, Inf))
end

"""
    RespirationParameters <: AbstractPhysiologyParameters

Physiological parameters relating to respiration of the organism.

# Parameters
- `fO2_extract::F` — Fraction of inspired oxygen extracted per breath (0–1).
- `pant` — Multiplier on breathing rate to simulate panting.
- `rq::F` — Respiratory quotient relating CO₂ produced to O₂ consumed (0–1).
- `Δ_breath::B` — Temperature offset between ambient air and exhaled air.
- `rh_exit::F` — Relative humidity of exhaled air (fraction 0–1).

These parameters influence radiative and evaporative exchange.
"""
Base.@kwdef struct RespirationParameters{OE,PA,RQ,DB,RE} <: AbstractPhysiologyParameters
    fO2_extract::OE = Param(0.2, bounds=(0.0, 1.0))
    pant::PA = Param(1.0, bounds=(1.0, Inf))
    rq::RQ = Param(0.8, bounds=(0.5, 1.2))
    Δ_breath::DB = Param(0.0u"K")
    rh_exit::RE = Param(1.0, bounds=(0.0, 1.0))
end

"""
    InsulationParameters(; kwargs...) <: AbstractMorphologyParameters

Parameters describing the physical and optical properties of fur or other
insulative pelage layers on an animal.

# Parameters
- `dorsal::FibreProperties` — Fibre properties for dorsal (back) surface
- `ventral::FibreProperties` — Fibre properties for ventral (belly) surface
- `conductivity_dorsal` — User-specified dorsal conductivity (W/m/K), or `nothing` to compute
- `conductivity_ventral` — User-specified ventral conductivity (W/m/K), or `nothing` to compute
- `depth_compressed` — Depth of compressed insulation during ground contact (mm)
- `fibre_conductivity` — Thermal conductivity of keratin fibres (W/m/K)
- `longwave_depth_fraction` — Fraction of insulation depth for longwave exchange (0-1)
"""
Base.@kwdef struct InsulationParameters{D,V,CD,CV,DC,FC,LDF} <: AbstractMorphologyParameters
    dorsal::D = FibreProperties(;
        diameter=Param(30.0u"μm"),
        length=Param(23.9u"mm"),
        density=Param(3000.0u"cm^-2"),
        depth=Param(2.0u"mm"),
        reflectance=Param(0.301, bounds=(0.0, 1.0)),
    )
    ventral::V = FibreProperties(;
        diameter=Param(30.0u"μm"),
        length=Param(23.9u"mm"),
        density=Param(3000.0u"cm^-2"),
        depth=Param(2.0u"mm"),
        reflectance=Param(0.301, bounds=(0.0, 1.0)),
    )
    conductivity_dorsal::CD = nothing
    conductivity_ventral::CV = nothing
    depth_compressed::DC = Param(2.0u"mm")
    fibre_conductivity::FC = Param(0.209u"W/m/K")
    longwave_depth_fraction::LDF = Param(1, bounds=(0.0, 1.0))
end

"""
    MetabolicParameters <: AbstractPhysiologyParameters

A collection of physiological parameters relating to metabolic rate.

# Keywords
- `T_core::B` — Core body temperature (K)
- `Q_metabolism::B` — Metabolic heat generation rate (W)
- `q10::F` — Q10 factor describing metabolic rate sensitivity to core temperature.
"""
Base.@kwdef struct MetabolismParameters{TC,QM,QT} <: AbstractPhysiologyParameters
    T_core::TC = Param(u"K"(37u"°C"))
    Q_metabolism::QM = Param(0.0u"W")
    q10::QT = Param(2.0)
    model = AndrewsPough2()
end
