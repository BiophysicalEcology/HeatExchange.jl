
# Behavioral rules here
abstract type AbstractBehavior end

abstract type AbstractBehavModel <: AbstractBehavior end 
        
abstract type AbstractBehavParameters <: AbstractBehavior end
        
abstract type AbstractBehavThresholds <: AbstractBehavior end

abstract type Diurnal <: AbstractBehavModel end
abstract type Nocturnal <: AbstractBehavModel end
abstract type Crepuscular <: AbstractBehavModel end
abstract type ShadeSeeker <: AbstractBehavModel end
abstract type Burrower <: AbstractBehavModel end
abstract type Climber <: AbstractBehavModel end

struct NullBehavior <: AbstractBehavModel end

initialise_state(::NullBehavior) = ()

abstract type BodyTemperatureThresholds end

struct Thermoregulate{B<:BodyTemperatureThresholds} <: AbstractBehavParameters
    T_F_min::T
    T_F_max::T
    T_B_min::T
    T_RB_min::T
    T_pref::T
end


#' @param T_F_min = 24, Minimum foraging temperature, °C (also affects burrow depth selection)
#' @param T_F_max = 34, Maximum foraging temperature, °C
#' @param T_B_min = 17.5, Minimum basking temperature, °C
#' @param T_RB_min = 17.5, Minimum temperature at which animal will move from retreat to basking site, °C
#' @param T_pref = 30, Preferred body temperature, °C
#' @param CT_max = 40, Critical thermal maximum, °C (affects burrow depth selection and may be used to impose death from heat stress)
#' @param CT_min = 6, Critical thermal minimum, °C (affects burrow depth selection and may be used to impose death from cold stress)
