

"""
    EndoThermoregParameters <: AbstractBehavParameters

Behavioural and physiological parameters governing thermoregulatory
responses such as panting, sweating, core-temperature adjustment, and
changes in body shape or tissue conductivity.

# Fields

- `insulation_step` — Incremental reduction in insulation depth from the
  piloerect (maximum fur depth) state.
- `shape_b_step` — Increment by which the body‐shape parameter `shape_b`
  increases per iteration, allowing the animal to uncurl toward
  `shape_b_max`.
- `shape_b_max` — Maximum allowed value of shape_b when adjusting posture.
- `T_core_max` — Maximum core temperature (K).
- `T_core_min` — Minimum core temperature during torpor (K).
- `T_core_step` — Increment by which core temperature is elevated per
  iteration (K).
- `k_flesh_step` — Increment in flesh thermal conductivity (W/m/K).
- `k_flesh_max` — Maximum flesh thermal conductivity (W/m/K).
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
Base.@kwdef struct EndoThermoregParameters{
    INS, SBS, SBM, TMX, TMN, TCS, KFS, KFM,
    PTS, PTM, PMX, SWS, SWM
} <: AbstractBehavParameters
    insulation_step::INS   = Param(1.0)
    shape_b_step::SBS      = Param(0.1)
    shape_b_max::SBM       = Param(5.0)
    T_core_max::TMX        = Param(39u"°C" |> u"K")
    T_core_min::TMN        = Param(19u"°C" |> u"K")
    T_core_step::TCS       = Param(0.1u"K")
    k_flesh_step::KFS      = Param(0.1u"W/m/K")
    k_flesh_max::KFM       = Param(2.8u"W/m/K")
    pant_step::PTS         = Param(0.1)
    pant_multiplier::PTM   = Param(1.05)
    pant_max::PMX          = Param(5.0)
    skin_wetness_step::SWS = Param(0.001)
    skin_wetness_max::SWM  = Param(1.0, bounds=(0.0, 1.0))
end

# # Behavioral rules here
# abstract type AbstractBehavior end

# abstract type AbstractBehavModel <: AbstractBehavior end 
        
# abstract type AbstractBehavParameters <: AbstractBehavior end
        
# abstract type AbstractBehavThresholds <: AbstractBehavior end

#abstract type Diurnal <: AbstractBehavModel end
#abstract type Nocturnal <: AbstractBehavModel end
#abstract type Crepuscular <: AbstractBehavModel end
#abstract type ShadeSeeker <: AbstractBehavModel end
#abstract type Burrower <: AbstractBehavModel end
#abstract type Climber <: AbstractBehavModel end

# abstract type ActivityPeriod end

# struct Diurnal <: ActivityPeriod end
# struct Norturnal <: ActivityPeriod end
# struct Crepuscular <: ActivityPeriod end

# """
#     ResponsiveActivity <: ActivityPeriod

#     ResponsiveActivity(isactive)

# # Arguments

# - `isactive` a `Function` or functor that recieves a `ModelParEnvironment`
#     object with the current system state, and decides whether to be "active"
#     or "innactive".
# """
# struct ResponsiveActivity{F} <: ActivityPeriod 
#     isactive::F
# end

# struct NullBehavior <: AbstractBehavModel end

# initialise_state(::NullBehavior) = ()

# abstract type BodyTemperatureThresholds end

# struct EctoThermoregParameters{B<:BodyTemperatureThresholds} <: AbstractBehavParameters
#     T_F_min::T
#     T_F_max::T
#     T_B_min::T
#     T_RB_min::T
#     T_pref::T
# end


# function photo_phase(phase::Diurnal, zenith_angle)
#     Body(shape, insulation, geometry(shape, insulation))
# end

# function photo_phase()

# end

#' @param T_F_min = 24, Minimum foraging temperature, °C (also affects burrow depth selection)
#' @param T_F_max = 34, Maximum foraging temperature, °C
#' @param T_B_min = 17.5, Minimum basking temperature, °C
#' @param T_RB_min = 17.5, Minimum temperature at which animal will move from retreat to basking site, °C
#' @param T_pref = 30, Preferred body temperature, °C
#' @param CT_max = 40, Critical thermal maximum, °C (affects burrow depth selection and may be used to impose death from heat stress)
#' @param CT_min = 6, Critical thermal minimum, °C (affects burrow depth selection and may be used to impose death from cold stress)
