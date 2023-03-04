module Ectotherm

using Unitful

using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R


include("GEOM.jl")
include("biophysics.jl")
include("heat_balance.jl")

end
