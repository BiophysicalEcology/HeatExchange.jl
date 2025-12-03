module HeatExchange

function __init__()\
    Unitful.register(HeatExchange)
end

using FluidProperties
using FluidProperties: wet_air_properties, dry_air_properties, vapour_pressure, 
    enthalpy_of_vaporisation, water_properties, atmospheric_pressure

using Unitful, UnitfulMoles, ModelParameters, Roots

export Shape, Cylinder, Sphere, Ellipsoid, Plate, LeopardFrog, DesertIguana

export bird_skin_area, bird_plumage_area, mammal_skin_area, mammal_fur_area

export Body, get_total_area, get_skin_area, get_evaporation_area, get_r_skin, get_r_insulation, get_r_flesh

export Insulation, CompositeInsulation, Naked, Fur, Fat

export Organism, EndoModelPars, BodyPars, InsulationPars, IntegumentPars, PhysioPars, 
    ThermoregulationPars, ThermoregulationVars, OrganismalVars

export EnvironmentalPars, EnvironmentalVars, EnvironmentalVarsVec

export surface_area, silhouette_area, SolarOrientation, Intermediate, ParallelToSun, NormalToSun

export geometry, shape, insulation, geometry_with_layers

export heat_balance, get_Tb, flip2vectors

export solar, radin, radout, evaporation, conduction, convection, nusselt_free, nusselt_forced 

export ectotherm, Tsurf_and_Tlung, respiration_ectotherm

export radiant_temperature, insulation_radiant_temperature, compressed_radiant_temperature

export endotherm, ellipsoid_endotherm, update_T_insulation!, solve_with_insulation!, solve_without_insulation!

export insulation_thermal_conductivity, insulation_properties, net_metabolic_heat

export simulsol, respiration_endotherm, mean_skin_temperature

export metabolic_rate, AndrewsPough2, Kleiber, McKechnieWolf, MetabolicRateEquation

include("organism.jl")
include("traits.jl")
include("insulation.jl")
include("geometry/geometry.jl")
include("geometry/plate.jl")
include("geometry/cylinder.jl")
include("geometry/sphere.jl")
include("geometry/ellipsoid.jl")
include("geometry/desert_iguana.jl")
include("geometry/leopard_frog.jl")
include("environment.jl")
include("biophysics.jl")
include("metabolism.jl")
include("ectotherm.jl")
include("endotherm.jl")

end
