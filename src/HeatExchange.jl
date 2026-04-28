module HeatExchange

using BiophysicalGeometry
using ConstructionBase: getproperties, setproperties
using FluidProperties
using ModelParameters
using Roots
using Unitful
using UnitfulMoles

using FluidProperties:
    wet_air_properties,
    dry_air_properties,
    vapour_pressure,
    enthalpy_of_vaporisation,
    water_properties,
    atmospheric_pressure,
    GasFractions,
    molCO₂,
    molH₂O,
    molO₂,
    molN₂

using BiophysicalGeometry: AbstractBody, shape, outer_insulation

export Organism,
    HeatExchangeTraits,
    SolveMetabolicRateOptions,
    InsulationParameters,
    ExternalConductionParameters,
    InternalConductionParameters,
    RadiationParameters,
    ConvectionParameters,
    AnimalEvaporationParameters,
    LeafEvaporationParameters,
    HydraulicParameters,
    RespirationParameters,
    MetabolismParameters

export EnvironmentalPars, EnvironmentalVars, EnvironmentalVarsVec

export body,
    traits,
    shape_pars,
    insulation_pars,
    conduction_pars_external,
    conduction_pars_internal,
    convection_pars,
    radiation_pars,
    evaporation_pars,
    hydraulic_pars,
    respiration_pars,
    metabolism_pars,
    options

export solar,
    radiation_in, radiation_out, evaporation, conduction, convection, nusselt_free, nusselt_forced

export heat_balance, solve_temperature, surface_and_lung_temperature

export radiant_temperature, insulation_radiant_temperature, compressed_radiant_temperature

export EvaluationStrategy, SingleBody, MultiSided, evaluation_strategy

export solve_metabolic_rate,
    ellipsoid_endotherm,
    solve_with_insulation!,
    solve_without_insulation!

export insulation_thermal_conductivity, insulation_properties, net_metabolic_heat

export solve_temperatures, mean_skin_temperature, respiration

export ConductanceCoeffs,
    DivisorCoeffs,
    RadiationCoeffs,
    BodyRegionValues,
    DorsalVentral,
    FibreProperties,
    EnvironmentTemperatures,
    OrganismTemperatures,
    GasFractions,
    ViewFactors,
    AtmosphericConditions,
    ThermalConductivities,
    MolarFluxes,
    HeatFlows,
    SolarConditions,
    TransferCoefficients,
    MetabolicRates,
    Emissivities,
    Absorptivities,
    InsulationProperties,
    GeometryVariables

export CharacteristicDimFormula, VolumeCubeRoot, ScaledDimension, characteristic_dimension

export MetabolicRateEquation, metabolic_rate, AndrewsPough2, Kleiber, McKechnieWolf, PlantDarkRespiration

export OxygenJoulesConversion, O2_to_Joules, Joules_to_O2, Typical, Kleiber1961

export zbrac, zbrent

export NLPStrategy, WeightedMeanNLP, MultiSidedNLP
export WeightedMeanNLPPacked, MultiSidedNLPPacked
export nlp_pack, nlp_residuals, nlp_assemble_output

include("rootfinding.jl")
include("organism.jl")
include("traits.jl")
include("environment.jl")

include("insulated/types.jl")
include("insulated/utils.jl")

include("biophysics.jl")
include("metabolism.jl")
include("respiration.jl")
include("insulated/insulation.jl")

include("internal_temperature.jl")
include("heat_balance.jl")
include("solve_temperature.jl")

include("ellipsoid_model.jl")
include("insulated/radiant_temperature.jl")
include("insulated/insulation_radiant_temperature.jl")
include("insulated/compressed_radiant_temperature.jl")
include("insulated/mean_skin_temperature.jl")
include("insulated/net_metabolic_heat.jl")
include("insulated/skin_and_insulation_temperature.jl")
include("solve_metabolic_rate.jl")
include("nlp_interface.jl")

end
