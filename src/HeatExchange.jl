module HeatExchange

function __init__()\
    Unitful.register(HeatExchange)
end

using Unitful, UnitfulMoles, ModelParameters, Roots

using FluidProperties
using FluidProperties: wet_air_properties, dry_air_properties, vapour_pressure, 
    enthalpy_of_vaporisation, water_properties, atmospheric_pressure

# to go in GeometryOfOrganisms.jl
export AbstractGeometryModel, AbstractGeometryPars, AbstractShape, AbstractInsulation, Body
export Cylinder, Sphere, Ellipsoid, Plate, LeopardFrog, DesertIguana
export CompositeInsulation, Naked, Fur, Fat
export geometry, shape, insulation
export bird_skin_area, bird_plumage_area, mammal_skin_area, mammal_fur_area
export get_total_area, get_skin_area, get_evaporation_area, get_r_skin, get_r_insulation, get_r_flesh
export surface_area, silhouette_area, SolarOrientation, Intermediate, ParallelToSun, NormalToSun

export Organism, EndoModelPars, BodyPars, InsulationPars, IntegumentPars, PhysioPars, 
    ThermoregulationPars, ThermoregulationVars, OrganismalVars

export EnvironmentalPars, EnvironmentalVars, EnvironmentalVarsVec

export heat_balance, get_Tb, flip2vectors

export solar, radin, radout, evaporation, conduction, convection, nusselt_free, nusselt_forced 

export ectotherm, Tsurf_and_Tlung, respiration_ectotherm

export radiant_temperature, insulation_radiant_temperature, compressed_radiant_temperature

export endotherm, ellipsoid_endotherm, update_T_insulation!, solve_with_insulation!, solve_without_insulation!

export insulation_thermal_conductivity, insulation_properties, net_metabolic_heat

export simulsol, respiration_endotherm, mean_skin_temperature, respiration

export MetabolicRateEquation, metabolic_rate, AndrewsPough2, Kleiber, McKechnieWolf

export OxygenJoulesConversion, O2_to_Joules, Joules_to_O2, Typical, Kleiber1961

include("geometry/geometry.jl")
include("geometry/plate.jl")
include("geometry/cylinder.jl")
include("geometry/sphere.jl")
include("geometry/ellipsoid.jl")
include("geometry/desert_iguana.jl")
include("geometry/leopard_frog.jl")

include("organism.jl")
include("traits.jl")
include("environment.jl")

include("biophysics.jl")
include("metabolism.jl")
include("respiration.jl")
include("insulation.jl")

include("ectotherm/lung_and_surface_temperature.jl")
include("ectotherm/ectotherm.jl")

include("endotherm/ellipsoid_model.jl")
include("endotherm/radiant_temperature.jl")
include("endotherm/insulation_radiant_temperature.jl")
include("endotherm/compressed_radiant_temperature.jl")
include("endotherm/mean_skin_temperature.jl")
include("endotherm/net_metabolic_heat.jl")
include("endotherm/skin_and_insulation_temperature.jl")
include("endotherm/endotherm.jl")

end
