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
export total_area, skin_area, evaporation_area, skin_radius, insulation_radius, flesh_radius
export surface_area, silhouette_area, SolarOrientation, Intermediate, ParallelToSun, NormalToSun

export Organism, Traits, EndoModelPars, BodyPars, InsulationParameters, ExternalConductionParameters,
  InternalConductionParameters, RadiationParameters, ConvectionParameters, EvaporationParameters, 
  HydraulicParameters, RespirationParameters, MetabolismParameters

export EnvironmentalPars, EnvironmentalVars, EnvironmentalVarsVec

export body, traits, shape, insulation, insulationpars, conductionpars_external, conductionpars_internal, 
  convectionpars, radiationpars, evaporationpars, hydraulicpars, respirationpars, metabolismpars

export heat_balance, get_Tb

export solar, radin, radout, evaporation, conduction, convection, nusselt_free, nusselt_forced 

export ectotherm, Tsurf_and_Tlung, respiration_ectotherm

export radiant_temperature, insulation_radiant_temperature, compressed_radiant_temperature

export endotherm, solve_metabolic_rate, ellipsoid_endotherm, update_T_insulation!, solve_with_insulation!, solve_without_insulation!

export insulation_thermal_conductivity, insulation_properties, net_metabolic_heat

export simulsol, respiration_endotherm, mean_skin_temperature, respiration

export MetabolicRateEquation, metabolic_rate, AndrewsPough2, Kleiber, McKechnieWolf

export OxygenJoulesConversion, O2_to_Joules, Joules_to_O2, Typical, Kleiber1961

include("geometry/geometry.jl")
include("geometry/shapes/plate.jl")
include("geometry/shapes/cylinder.jl")
include("geometry/shapes/sphere.jl")
include("geometry/shapes/ellipsoid.jl")
include("geometry/shapes/desert_iguana.jl")
include("geometry/shapes/leopard_frog.jl")

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
