
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

# struct Thermoregulate{B<:BodyTemperatureThresholds} <: AbstractBehavParameters
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
