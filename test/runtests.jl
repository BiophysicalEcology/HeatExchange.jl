using Ectotherm
using Test
using Plots
using UnitfulRecipes
using Unitful: °, rad, °F, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R
#import Base: promote_rule, convert

# abstract type Temperature end

# struct Celsius <: Temperature
#     value::Float64
# end

# struct Kelvin <: Temperature
#    value::Float64
# end

# promote_rule(::Type{Kelvin}, ::Type{Celsius})     = Kelvin

# convert(::Type{Kelvin},  t::Celsius)     = Kelvin(t.value + 273.15)
# convert(::Type{Celsius}, t::Kelvin)      = Celsius(t.value - 273.15)
# promote(Kelvin(1), Celsius(1))
# promote(Celsius(1), Kelvin(1))

#include("biophysics.jl") end
#include("geometry.jl") end
#include("ectotherm.jl") end
