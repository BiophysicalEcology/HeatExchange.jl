module HeatExchange

function __init__()\
    Unitful.register(HeatExchange)
end

using FluidProperties: wet_air_properties, dry_air_properties, vapour_pressure, enthalpy_of_vaporisation, water_properties, atmospheric_pressure

using Unitful, UnitfulMoles, ModelParameters

using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

export Shape, Cylinder, Sphere, Ellipsoid, Plate, LeopardFrog, DesertIguana

export bird_skin_area, bird_plumage_area, mammal_skin_area, mammal_fur_area

export Body, get_total_area, get_skin_area, get_convection_area, get_r_skin, get_r_insulation, get_r_flesh

export Insulation, CompositeInsulation, Naked, Fur, Fat

export Organism, MorphoPars, PhysioPars, OrganismalVars

export EnvironmentalPars, EnvironmentalVars, EnvironmentalVarsVec, InsulationPars

export surface_area, silhouette_area

export geometry, shape, insulation, geometry_with_layers

export heat_balance, get_Tb, flip2vectors

export metabolism, respiration, solar, radin, radout, evaporation, conduction, convection

export nusselt_free, nusselt_forced, Tsurf_and_Tlung, radiant_temperature

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
