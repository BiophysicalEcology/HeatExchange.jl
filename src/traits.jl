
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
- `flesh_conductivity::K` — Thermal conductivity of lean tissue (W/m/K).
- `fat_conductivity::K` — Thermal conductivity of fat tissue (W/m/K).
- `fat_density` — Density of fat tissue (kg/m³).

"""
Base.@kwdef struct InternalConductionParameters{FF,FL,FA,DF} <: AbstractPhysiologyParameters
    fat_fraction::FF = Param(0.0, bounds=(0.0, 1.0))
    flesh_conductivity::FL = Param(0.9u"W/m/K")
    fat_conductivity::FA = Param(0.230u"W/m/K")
    fat_density::DF = Param(901.0u"kg/m^3")
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
- `body_absorptivity_dorsal::F` — Shortwave absorptivity of dorsal body surface (0–1).
- `body_absorptivity_ventral::F` — Shortwave absorptivity of ventral body surface (0–1).
- `body_emissivity_dorsal::F` — Longwave emissivity of dorsal body surface (0–1).
- `body_emissivity_ventral::F` — Longwave emissivity of ventral body surface (0–1).
- `sky_view_factor::F` — Radiative view factor to sky (0–1).
- `ground_view_factor::F` — Radiative view factor to substrate/ground (0–1).
- `vegetation_view_factor::F` — View factor to surrounding vegetation (0–1).
- `bush_view_factor::F` — View factor to vegetation at animal height (0–1).
- `ventral_fraction::F` — Fraction of total surface area that is ventral (0–1).
- `silhouette_area` — Projected silhouette area (m²).
- `total_area` — Total surface area (m²).
- `conduction_area` — Area in contact with substrate (m²).

"""
Base.@kwdef struct RadiationParameters{AD,AV,ED,EV,AS,AT,AC,FS,FG,FV,FB,VF} <:
                   AbstractMorphologyParameters
    body_absorptivity_dorsal::AD = Param(0.85, bounds=(0.0, 1.0))
    body_absorptivity_ventral::AV = Param(0.85, bounds=(0.0, 1.0))
    body_emissivity_dorsal::ED = Param(0.95, bounds=(0.0, 1.0))
    body_emissivity_ventral::EV = Param(0.95, bounds=(0.0, 1.0))
    silhouette_area::AS = Param(0.0u"m^2")
    total_area::AT = Param(0.0u"m^2")
    conduction_area::AC = Param(0.0u"m^2")
    sky_view_factor::FS = Param(0.5, bounds=(0.0, 1.0))
    ground_view_factor::FG = Param(0.5, bounds=(0.0, 1.0))
    vegetation_view_factor::FV = Param(0.0, bounds=(0.0, 1.0))
    bush_view_factor::FB = Param(0.0, bounds=(0.0, 1.0))
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
    LeafEvaporationParameters <: AbstractMorphologyParameters

Evaporation parameters for a leaf using stomatal vapour conductances (mol/m²/s).
The boundary layer conductance is supplied via `mass::TransferCoefficients` from
a prior `convection()` call, ensuring the same body geometry and convective
enhancement factor are used for both heat and mass transfer.

# Parameters
- `abaxial_vapour_conductance` — Stomatal conductance of the abaxial (lower) surface (mol/m²/s); default 0.3.
- `adaxial_vapour_conductance` — Stomatal conductance of the adaxial (upper) surface (mol/m²/s); default 0.0.
- `cuticular_conductance` — Baseline conductance when stomata are closed (mol/m²/s); added equally
  to both surfaces (`cuticular_conductance / 2` per side). Equivalent to ~0.1 µmol H₂O/(m²·s·Pa); default 0.01.
"""
Base.@kwdef struct LeafEvaporationParameters{AB,AD,CC} <: AbstractMorphologyParameters
    abaxial_vapour_conductance::AB = Param(0.3u"mol/m^2/s",  bounds=(0.0, Inf))
    adaxial_vapour_conductance::AD = Param(0.0u"mol/m^2/s",  bounds=(0.0, Inf))
    cuticular_conductance::CC      = Param(0.01u"mol/m^2/s", bounds=(0.0, Inf))
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
- `oxygen_extraction_efficiency` — Fraction of inspired oxygen extracted per breath (0–1).
- `pant` — Multiplier on breathing rate to simulate panting.
- `respiratory_quotient` — Respiratory quotient relating CO₂ produced to O₂ consumed (0–1).
- `exhaled_temperature_offset` — Temperature offset between ambient air and exhaled air (K).
- `exhaled_relative_humidity` — Relative humidity of exhaled air (fraction 0–1).
- `mouth_fraction` — Fractional addition to `skin_wetness` when mouth is open for panting
  (`pant > 1`). Replicates NicheMapR PMOUTH parameter. Default 0.0 (adds 0% of surface area).

These parameters influence radiative and evaporative exchange.
"""
Base.@kwdef struct RespirationParameters{OE,PA,RQ,DB,RE,MF} <: AbstractPhysiologyParameters
    oxygen_extraction_efficiency::OE = Param(0.2, bounds=(0.0, 1.0))
    pant::PA = Param(1.0, bounds=(1.0, Inf))
    respiratory_quotient::RQ = Param(0.8, bounds=(0.5, 1.2))
    exhaled_temperature_offset::DB = Param(0.0u"K")
    exhaled_relative_humidity::RE = Param(1.0, bounds=(0.0, 1.0))
    mouth_fraction::MF = Param(0.0, bounds=(0.0, 1.0))
end

"""
    InsulationParameters(; kwargs...) <: AbstractMorphologyParameters

Parameters describing the physical and optical properties of fur or other
insulative pelage layers on an animal.

# Parameters
- `dorsal::FibreProperties` — Fibre properties for dorsal (back) surface
- `ventral::FibreProperties` — Fibre properties for ventral (belly) surface
- `depth_compressed` — Depth of compressed insulation during ground contact (mm)
- `longwave_depth_fraction` — Fraction of insulation depth for longwave exchange (0-1)
"""
Base.@kwdef struct InsulationParameters{D,V,DC,LDF} <: AbstractMorphologyParameters
    dorsal::D = FibreProperties(;
        diameter=Param(30.0u"μm"),
        length=Param(23.9u"mm"),
        density=Param(3000.0u"cm^-2"),
        depth=Param(2.0u"mm"),
        reflectance=Param(0.301, bounds=(0.0, 1.0)),
        conductivity=Param(0.209u"W/m/K"),
    )
    ventral::V = FibreProperties(;
        diameter=Param(30.0u"μm"),
        length=Param(23.9u"mm"),
        density=Param(3000.0u"cm^-2"),
        depth=Param(2.0u"mm"),
        reflectance=Param(0.301, bounds=(0.0, 1.0)),
        conductivity=Param(0.209u"W/m/K"),
    )
    depth_compressed::DC = Param(2.0u"mm")
    longwave_depth_fraction::LDF = Param(1, bounds=(0.0, 1.0))
end

"""
    MetabolicParameters <: AbstractPhysiologyParameters

A collection of physiological parameters relating to metabolic rate.

# Keywords
- `core_temperature::B` — Core body temperature (K)
- `metabolic_heat_flow::B` — Metabolic heat generation rate (W)
- `q10::F` — Q10 factor describing metabolic rate sensitivity to core temperature.
"""
Base.@kwdef struct MetabolismParameters{TC,QM,QT} <: AbstractPhysiologyParameters
    core_temperature::TC = Param(u"K"(37u"°C"))
    metabolic_heat_flow::QM = Param(0.0u"W")
    q10::QT = Param(2.0)
    model = AndrewsPough2()
end
