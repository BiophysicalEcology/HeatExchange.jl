# API

## Organisms and traits

```@docs
Organism
HeatExchangeTraits
body
traits
shape_pars
insulation_pars
conduction_pars_external
conduction_pars_internal
radiation_pars
convection_pars
evaporation_pars
hydraulic_pars
respiration_pars
metabolism_pars
options
```

## Parameters

```@docs
RadiationParameters
ConvectionParameters
ExternalConductionParameters
InternalConductionParameters
AnimalEvaporationParameters
LeafEvaporationParameters
HydraulicParameters
RespirationParameters
MetabolismParameters
InsulationParameters
FibreProperties
SolveMetabolicRateOptions
```

## Environment

```@docs
EnvironmentalPars
EnvironmentalVars
EnvironmentalVarsVec
```

## Solvers

```@docs
solve_temperature
solve_metabolic_rate
heat_balance
ThermoregulationOutput
ThermoregulationState
MorphologyState
EnergyFlowState
MassFlowState
EvaluationStrategy
SingleBody
MultiSided
evaluation_strategy
```

### CommonSolve interface

`init`, `solve!` and `solve` are those of [CommonSolve.jl](https://github.com/SciML/CommonSolve.jl), with methods
for a [`HeatBalanceProblem`](@ref).

```@docs
HeatBalanceProblem
HeatBalanceSolver
reinit!
```

## Heat flows

```@docs
solar
radiation_in
radiation_out
convection
nusselt_free
nusselt_forced
conduction
evaporation
respiration
surface_and_lung_temperature
CharacteristicDimFormula
VolumeCubeRoot
ScaledDimension
characteristic_dimension
```

## Metabolism

```@docs
metabolic_rate
MetabolicRateEquation
AndrewsPough2
Kleiber
McKechnieWolf
PlantDarkRespiration
OxygenJoulesConversion
Typical
Kleiber1961
O2_to_Joules
Joules_to_O2
```

## Insulation

```@docs
insulation_properties
insulation_thermal_conductivity
InsulationProperties
solve_temperatures
radiant_temperature
insulation_radiant_temperature
compressed_radiant_temperature
mean_skin_temperature
net_metabolic_heat
```

## Radial layers

```@docs
AbstractRadialLayer
GeneratingCore
ConductiveShell
stack_resistance
core_to_skin_stack
radial_net_metabolic_heat
```

## Bodies of many parts

```@docs
solve_part_heat_balance
solve_part_surface
part_surface_residuals
solve_coupled_metabolic_rate
HeatCoupling
SharedCore
ConductiveCoupling
CompartmentGraph
compartment_graph
num_compartments
compartment_part_names
compartment_of
parts_in_compartment
contribution_to_conductance
contribution_to_heat_load
build_conductance_matrix
solve_core_temperatures
solve_regulated_core_temperatures
```

## Smoothing and the NLP interface

```@docs
SmoothingStrategy
HardBound
SmoothBound
safe_abs
safe_relu
safe_step
safe_max
safe_min
safe_clamp
NLPStrategy
nlp_pack
```

## The ellipsoid model

```@docs
ellipsoid_endotherm
```

## Root finding

```@docs
zbrent
zbrac
```

## Containers

These hold the inputs and outputs of the heat-flow functions.

```@docs
BodySide
Dorsal
Ventral
Fluid
Air
Water
DorsalVentral
BodyRegionValues
ViewFactors
Absorptivities
Emissivities
SolarConditions
AtmosphericConditions
EnvironmentTemperatures
OrganismTemperatures
ThermalConductivities
TransferCoefficients
MetabolicRates
MolarFluxes
HeatFlows
GeometryVariables
ConductanceCoeffs
DivisorCoeffs
RadiationCoeffs
```

## Example parameters

```@docs
example_heat_exchange_traits
example_shape_pars
example_ellipsoid_shape_pars
example_insulation_pars
example_conduction_pars_external
example_conduction_pars_internal
example_radiation_pars
example_convection_pars
example_evaporation_pars
example_leaf_evaporation_pars
example_hydraulic_pars
example_respiration_pars
example_metabolism_pars
example_metabolic_rate_options
example_environment_vars
example_environment_pars
example_ectotherm_heat_exchange_traits
example_ectotherm_conduction_pars_external
example_ectotherm_conduction_pars_internal
example_ectotherm_radiation_pars
example_ectotherm_evaporation_pars
example_ectotherm_hydraulic_pars
example_ectotherm_respiration_pars
example_ectotherm_metabolism_pars
```
