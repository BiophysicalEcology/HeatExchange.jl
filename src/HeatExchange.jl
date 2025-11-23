module HeatExchange

function __init__()\
    Unitful.register(HeatExchange)
end

using FluidProperties: wet_air_properties, dry_air_properties, vapour_pressure, enthalpy_of_vaporisation, water_properties, atmospheric_pressure

using Unitful, UnitfulMoles, ModelParameters

using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

export Shape, Cylinder, Sphere, Ellipsoid, Plate, LeopardFrog, DesertIguana

export bird_skin_area, bird_plumage_area, mammal_skin_area, mammal_fur_area

export Body

export Insulation, CompositeInsulation, Naked, Fur, Fat

export Organism, MorphoPars, PhysioPars, OrganismalVars

export EnvironmentalPars, EnvironmentalVars, EnvironmentalVarsVec

export calc_area, calc_silhouette_area

export geometry, shape, insulation, geometry_with_layers

export heat_balance, get_Tb, flip2vectors

export metabolism, respiration, solar, radin, radout, evaporation, conduction, convection, get_nusselt_free, get_nusselt_forced, get_Tsurf_Tlung

export ellipsoid_endotherm

export insulation_thermal_conductivity, insulation_properties

include("geometry.jl")
include("insulation.jl")
include("environment.jl")
include("organism.jl")
include("biophysics.jl")
include("ectotherm.jl")
include("endotherm.jl")

@compound H2O
@compound O2
@compound CO2
@compound N2

end
