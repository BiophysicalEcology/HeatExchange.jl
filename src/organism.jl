abstract type AbstractPhysiologyModel end

"""
    MetabolicRateEquation

Abstract supertype for equations that give metabolic rate from mass and temperature, used with
[`metabolic_rate`](@ref). Subtypes are [`AndrewsPough2`](@ref), [`Kleiber`](@ref), [`McKechnieWolf`](@ref) and
[`PlantDarkRespiration`](@ref).
"""
abstract type MetabolicRateEquation <: AbstractPhysiologyModel end
"""
    OxygenJoulesConversion

Abstract supertype for conversions between oxygen consumption and heat production, used with
[`O2_to_Joules`](@ref) and [`Joules_to_O2`](@ref). Subtypes are [`Typical`](@ref) and [`Kleiber1961`](@ref).
"""
abstract type OxygenJoulesConversion <: AbstractPhysiologyModel end

abstract type AbstractPhysiologyParameters end

abstract type AbstractMorphologyParameters end

abstract type AbstractModelParameters end

abstract type AbstractFunctionalTraits end

"""
    shape_pars(organism)
    shape_pars(traits)

The shape of an [`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
shape_pars(t::AbstractFunctionalTraits) = stripparams(t.shape_pars)
"""
    insulation_pars(organism)
    insulation_pars(traits)

The insulation parameters ([`InsulationParameters`](@ref)) of an [`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
insulation_pars(t::AbstractFunctionalTraits) = stripparams(t.insulation_pars)
"""
    conduction_pars_external(organism)
    conduction_pars_external(traits)

The parameters for conduction to the substrate ([`ExternalConductionParameters`](@ref)) of an
[`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
function conduction_pars_external(t::AbstractFunctionalTraits)
    stripparams(t.conduction_pars_external)
end
"""
    conduction_pars_internal(organism)
    conduction_pars_internal(traits)

The parameters for conduction within the body ([`InternalConductionParameters`](@ref)) of an
[`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
function conduction_pars_internal(t::AbstractFunctionalTraits)
    stripparams(t.conduction_pars_internal)
end
"""
    convection_pars(organism)
    convection_pars(traits)

The convection parameters ([`ConvectionParameters`](@ref)) of an [`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
convection_pars(t::AbstractFunctionalTraits) = stripparams(t.convection_pars)
"""
    radiation_pars(organism)
    radiation_pars(traits)

The radiation parameters ([`RadiationParameters`](@ref)) of an [`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
radiation_pars(t::AbstractFunctionalTraits) = stripparams(t.radiation_pars)
"""
    evaporation_pars(organism)
    evaporation_pars(traits)

The evaporation parameters ([`AnimalEvaporationParameters`](@ref) or [`LeafEvaporationParameters`](@ref)) of an [`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
evaporation_pars(t::AbstractFunctionalTraits) = stripparams(t.evaporation_pars)
"""
    hydraulic_pars(organism)
    hydraulic_pars(traits)

The hydraulic parameters ([`HydraulicParameters`](@ref)) of an [`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
hydraulic_pars(t::AbstractFunctionalTraits) = stripparams(t.hydraulic_pars)
"""
    respiration_pars(organism)
    respiration_pars(traits)

The respiration parameters ([`RespirationParameters`](@ref)) of an [`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
respiration_pars(t::AbstractFunctionalTraits) = stripparams(t.respiration_pars)
"""
    metabolism_pars(organism)
    metabolism_pars(traits)

The metabolism parameters ([`MetabolismParameters`](@ref)) of an [`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
metabolism_pars(t::AbstractFunctionalTraits) = stripparams(t.metabolism_pars)
"""
    options(organism)
    options(traits)

The solver options ([`SolveMetabolicRateOptions`](@ref)) of an [`Organism`](@ref) or a [`HeatExchangeTraits`](@ref), with any `Param` wrappers removed.
"""
options(t::AbstractFunctionalTraits) = stripparams(t.options)

# TODO more specific subtypes
"""
    HeatExchangeTraits(shape_pars, insulation_pars, conduction_pars_external, conduction_pars_internal,
                       radiation_pars, convection_pars, evaporation_pars, hydraulic_pars, respiration_pars,
                       metabolism_pars, options)

The functional traits of an organism that its heat budget needs, one set of parameters for each process.

# Fields
- `shape_pars` — the shape of the body, an `AbstractShape` of BiophysicalGeometry.jl
- `insulation_pars` — [`InsulationParameters`](@ref)
- `conduction_pars_external` — [`ExternalConductionParameters`](@ref)
- `conduction_pars_internal` — [`InternalConductionParameters`](@ref)
- `radiation_pars` — [`RadiationParameters`](@ref)
- `convection_pars` — [`ConvectionParameters`](@ref)
- `evaporation_pars` — [`AnimalEvaporationParameters`](@ref) or [`LeafEvaporationParameters`](@ref)
- `hydraulic_pars` — [`HydraulicParameters`](@ref)
- `respiration_pars` — [`RespirationParameters`](@ref)
- `metabolism_pars` — [`MetabolismParameters`](@ref)
- `options` — [`SolveMetabolicRateOptions`](@ref)

See [`example_heat_exchange_traits`](@ref) and [`example_ectotherm_heat_exchange_traits`](@ref) for ready-made sets.
"""
struct HeatExchangeTraits{
    SP<:AbstractShape,
    IN<:AbstractMorphologyParameters,
    CE<:AbstractMorphologyParameters,
    CI<:AbstractPhysiologyParameters,
    RA<:AbstractMorphologyParameters,
    CO<:AbstractMorphologyParameters,
    EV<:AbstractMorphologyParameters,
    HD<:AbstractPhysiologyParameters,
    RE<:AbstractPhysiologyParameters,
    ME<:AbstractPhysiologyParameters,
    OP<:AbstractModelParameters,
} <: AbstractFunctionalTraits
    shape_pars::SP
    insulation_pars::IN
    conduction_pars_external::CE
    conduction_pars_internal::CI
    radiation_pars::RA
    convection_pars::CO
    evaporation_pars::EV
    hydraulic_pars::HD
    respiration_pars::RE
    metabolism_pars::ME
    options::OP
end

"""
    AbstractOrganism

Abstract supertype for organisms.
"""
abstract type AbstractOrganism end

# With some generic methods to get the params and body
"""
    body(organism)

The body of an [`Organism`](@ref), a `Body` or `CompositeBody` of BiophysicalGeometry.jl.
"""
body(o::AbstractOrganism) = o.body # gets the body from an object of type AbstractOrganism
"""
    traits(organism)

The traits of an [`Organism`](@ref), a [`HeatExchangeTraits`](@ref).
"""
traits(o::AbstractOrganism) = o.traits
#shape(o::AbstractOrganism) = shape(body(o)) # gets the shape from an object of type AbstractOrganism
#insulation(o::AbstractOrganism) = insulation(body(o)) # gets the insulation from an object of type AbstractOrganism

# Forwarding methods from organism to traits
shape_pars(o::AbstractOrganism) = shape_pars(traits(o))
insulation_pars(o::AbstractOrganism) = insulation_pars(traits(o))
conduction_pars_external(o::AbstractOrganism) = conduction_pars_external(traits(o))
conduction_pars_internal(o::AbstractOrganism) = conduction_pars_internal(traits(o))
convection_pars(o::AbstractOrganism) = convection_pars(traits(o))
radiation_pars(o::AbstractOrganism) = radiation_pars(traits(o))
evaporation_pars(o::AbstractOrganism) = evaporation_pars(traits(o))
hydraulic_pars(o::AbstractOrganism) = hydraulic_pars(traits(o))
respiration_pars(o::AbstractOrganism) = respiration_pars(traits(o))
metabolism_pars(o::AbstractOrganism) = metabolism_pars(traits(o))
options(o::AbstractOrganism) = options(traits(o))

"""
    Organism <: AbstractOrganism

    Organism(body, traits)

A concrete implementation of `AbstractOrganism`. It accepts an `AbstractBody` of BiophysicalGeometry.jl and an
`AbstractFunctionalTraits` object, such as [`HeatExchangeTraits`](@ref).
"""
struct Organism{B<:AbstractBody,T<:AbstractFunctionalTraits} <: AbstractOrganism
    body::B
    traits::T
end

# TODO use this as a container for outputs? Or remove?
"""
    AbstractOrganismalVars

Abstract supertype for organismal variables.
"""
abstract type AbstractOrganismalVars end

"""
    OrganismalVars <: AbstractOrganismalVars

    - `water_potential` — Body water potential (determines humidity at skin surface
    and liquid water exchange) (J/kg).
Variables for an [`AbstractOrganism`](@ref) model.
"""
Base.@kwdef mutable struct OrganismalVars{TC,TS,TI,TL,P} <: AbstractOrganismalVars
    core_temperature::TC
    skin_temperature::TS
    insulation_temperature::TI = core_temperature
    lung_temperature::TL
    water_potential::P
end
