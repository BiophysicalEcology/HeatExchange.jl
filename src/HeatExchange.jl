module HeatExchange

using Unitful, ModelParameters

using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

export Shape, Cylinder, Ellipsoid, Plate, LeopardFrog, DesertIguana 

export Body

export Insulation, Naked, Fur

export Organism, FunctionalTraits, OrganismalVars

export EnvironmentalParams, EnvironmentalVars

export calc_area, calc_silhouette_area

export geometry, shape, insulation

export heat_balance

export metabolism, respiration, solar, radin, radout, evaporation, conduction, convection, vapour_pressure, 
    wet_air, dry_air, water_prop, get_nusselt_free, get_nusselt_forced, get_Tsurf_Tlung

include("geometry.jl")
include("environment.jl")
include("organism.jl")
include("biophysics.jl")
include("heat_balance.jl")

end
