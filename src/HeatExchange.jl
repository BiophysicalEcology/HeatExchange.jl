module HeatExchange

using Unitful, UnitfulMoles, ModelParameters, Roots

using FluidProperties
using FluidProperties: wet_air_properties, dry_air_properties, vapour_pressure, 
    enthalpy_of_vaporisation, water_properties, atmospheric_pressure
using FluidProperties: wet_air_properties, dry_air_properties, vapour_pressure, 
    enthalpy_of_vaporisation, water_properties, atmospheric_pressure

using FluidProperties: molCO₂, molH₂O, molO₂, molN₂

using BiophysicalGeometry
using BiophysicalGeometry: AbstractBody, shape

export Organism, HeatExchangeTraits, EndoModelPars, InsulationParameters, ExternalConductionParameters,
  InternalConductionParameters, RadiationParameters, ConvectionParameters, EvaporationParameters, 
  HydraulicParameters, RespirationParameters, MetabolismParameters

export EnvironmentalPars, EnvironmentalVars, EnvironmentalVarsVec

export body, traits, insulationpars, conductionpars_external, conductionpars_internal, 
  convectionpars, radiationpars, evaporationpars, hydraulicpars, respirationpars, metabolismpars

export get_Tb

export solar, radin, radout, evaporation, conduction, convection, nusselt_free, nusselt_forced 

export ectotherm, Tsurf_and_Tlung

export radiant_temperature, insulation_radiant_temperature, compressed_radiant_temperature

export solve_metabolic_rate, ellipsoid_endotherm, update_T_insulation!, solve_with_insulation!, solve_without_insulation!

export insulation_thermal_conductivity, insulation_properties, net_metabolic_heat

export simulsol, mean_skin_temperature, respiration

export MetabolicRateEquation, metabolic_rate, AndrewsPough2, Kleiber, McKechnieWolf

export OxygenJoulesConversion, O2_to_Joules, Joules_to_O2, Typical, Kleiber1961

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

# function __init__()\
#     Unitful.register(HeatExchange)
# end

end
