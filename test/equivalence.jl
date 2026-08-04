# Phase 7 equivalence gate: the multi-part coupled solve reproduces the
# dorsal/ventral MultiSided result for the one case dorsal/ventral can represent.
#
# Setup is maximally symmetric so the reformulation is isolated from confounders:
# sky_view_factor = 0.5 and shade = 0 give dmult = vmult = 0.5, insulation is
# identical dorsal/ventral, there is no solar load, and no ground contact. In this
# regime HalfCylinder(mass/2) is geometrically exactly half of Cylinder(mass) —
# same radius and length, half the volume — so its flesh conductance, surface area,
# and radius-based convection all scale by exactly 1/2. The plain sum of two such
# parts should therefore reproduce the 1/2 + 1/2 weighted mean the dorsal/ventral
# path computes.

using HeatExchange
using BiophysicalGeometry
using Unitful, UnitfulMoles
using FluidProperties
using Test

using HeatExchange: solve_coupled_metabolic_rate,
    EnvironmentTemperatures, ViewFactors, AtmosphericConditions,
    example_environment_vars, example_environment_pars,
    example_insulation_pars, example_conduction_pars_external,
    example_conduction_pars_internal, example_evaporation_pars,
    example_radiation_pars, example_respiration_pars, example_metabolism_pars,
    example_heat_exchange_traits

# --- shared physiology -----------------------------------------------------
mass = 1.0u"kg"
ρ    = 1000.0u"kg/m^3"
b    = 3.0

ins_pars = example_insulation_pars()          # dorsal == ventral by default
fibre    = ins_pars.dorsal
fur      = FibrousLayer(fibre.depth, fibre.diameter, fibre.density)
fat      = FatLayer(0.0, 901.0u"kg/m^3")       # fat_fraction = 0

internal  = example_conduction_pars_internal()
external  = example_conduction_pars_external(; conduction_fraction = 0.0)
evap      = example_evaporation_pars()
rad_pars  = example_radiation_pars()           # symmetric: sky = ground = 0.5, ε = 0.99
resp_pars = example_respiration_pars()
metab     = example_metabolism_pars()
core      = metab.core_temperature

env_vars = example_environment_vars()          # uniform temps, global_radiation = 0, shade = 0
env_pars = example_environment_pars()

skin0  = core - 5u"K"
insul0 = env_vars.air_temperature + 2u"K"

# --- baseline: dorsal/ventral MultiSided on the full cylinder --------------
full_shape = Cylinder(mass, ρ, b)
full_body  = Body(full_shape, CompositeInsulation(fur, fat))
traits = example_heat_exchange_traits(;
    shape_pars = full_shape,
    insulation_pars = ins_pars,
    conduction_pars_external = external,
    conduction_pars_internal = internal,
    radiation_pars = rad_pars,
    evaporation_pars = evap,
    respiration_pars = resp_pars,
    metabolism_pars = metab,
)
organism = Organism(full_body, traits)
environment = (; environment_pars = env_pars, environment_vars = env_vars)
baseline = solve_metabolic_rate(organism, environment, skin0, insul0)

# --- multi-part: two HalfCylinder(mass/2) parts joined by a shared core ----
packed(view_factors) = (;
    temperature = EnvironmentTemperatures(
        env_vars.air_temperature, env_vars.sky_temperature,
        env_vars.ground_temperature, env_vars.vegetation_temperature,
        env_vars.bush_temperature, env_vars.substrate_temperature,
    ),
    view_factors,
    atmos = AtmosphericConditions(env_vars),
    fluid = env_pars.fluid,
    solar_flow = 0.0u"W",
    gas_fractions = env_pars.gas_fractions,
    convection_enhancement = env_pars.convection_enhancement,
)

part_traits(ϵ_body) = (;
    core_temperature   = core,
    flesh_conductivity = internal.flesh_conductivity,
    fat_conductivity   = internal.fat_conductivity,
    ϵ_body,
    skin_wetness       = evap.skin_wetness,
    insulation_wetness = evap.insulation_wetness,
    bare_skin_fraction = evap.bare_skin_fraction,
    eye_fraction       = evap.eye_fraction,
)

# The two halves reassemble the full cylinder, so their external curvature — hence
# the convective characteristic dimension — is the full cylinder's, not the value
# VolumeCubeRoot derives from a half's reduced volume. (Per §3.6 this is a Tier-1
# cache input the caller supplies.)
full_characteristic_dim = characteristic_dimension(VolumeCubeRoot(), full_body)

function part_setup(view_factors, ϵ_body)
    part_body = Body(HalfCylinder(mass / 2, ρ, b), CompositeInsulation(fur, fat))
    # Covered join patch: the flat face that mates with the other half
    # (2·radius_skin·length_skin), matching CompositeBody's covered-area accounting.
    covered_area = 2 * part_body.geometry.length.radius_skin * part_body.geometry.length.length_skin
    return (;
        body = part_body,
        insulation_pars = ins_pars,
        traits = part_traits(ϵ_body),
        environment_vars = packed(view_factors),
        conduction_fraction = 0.0,
        conductance_coefficient = 0.0u"W/K",
        ventral_fraction = 0.5,
        longwave_depth_fraction = 1.0,
        covered_area,
        characteristic_dim = full_characteristic_dim,
    )
end

# Per-side view factors mirror _pack_sides: sky/vegetation-facing (dorsal) and
# ground/bush-facing (ventral), each doubled (half the body, twice the exposure).
sky_factor    = rad_pars.sky_view_factor
ground_factor = 1 - sky_factor
dorsal_setup  = part_setup(ViewFactors(sky_factor * 2, 0.0, 0.0, 0.0), rad_pars.body_emissivity_dorsal)
ventral_setup = part_setup(ViewFactors(0.0, ground_factor * 2, 0.0, 0.0), rad_pars.body_emissivity_ventral)

multipart = solve_coupled_metabolic_rate(;
    part_surface_setups = (dorsal_setup, ventral_setup),
    core_temperature = core,
    skin_temperature = skin0,
    insulation_temperature = insul0,
    temperature_tolerance = 1e-3u"K",
    respire = true,
    respiration_pars = resp_pars,
    lung_mass = mass,
    air_temperature = env_vars.air_temperature,
    atmos = AtmosphericConditions(env_vars),
    gas_fractions = env_pars.gas_fractions,
    metabolic_heat_flow_setpoint = metab.metabolic_heat_flow,
    resp_tolerance = 1e-5,
)

# --- strong gate: the whole body as a single part reproduces dorsal/ventral ---
# solve_part_surface on the full cylinder (seeing sky+ground with the un-doubled
# view factors) is the faithful per-part reformulation of one dorsal/ventral side;
# in the symmetric case the two sides average to the full body, so a single full
# part must reproduce the dorsal/ventral solve essentially exactly.
whole_setup = (;
    body = full_body,
    insulation_pars = ins_pars,
    traits = part_traits(rad_pars.body_emissivity_dorsal),
    environment_vars = packed(ViewFactors(sky_factor, ground_factor, 0.0, 0.0)),
    conduction_fraction = 0.0,
    conductance_coefficient = 0.0u"W/K",
    ventral_fraction = 0.5,
    longwave_depth_fraction = 1.0,
)
whole = solve_coupled_metabolic_rate(;
    part_surface_setups = (whole_setup,),
    core_temperature = core,
    skin_temperature = skin0,
    insulation_temperature = insul0,
    temperature_tolerance = 1e-3u"K",
    respire = true,
    respiration_pars = resp_pars,
    lung_mass = mass,
    air_temperature = env_vars.air_temperature,
    atmos = AtmosphericConditions(env_vars),
    gas_fractions = env_pars.gas_fractions,
    metabolic_heat_flow_setpoint = metab.metabolic_heat_flow,
    resp_tolerance = 1e-5,
)

@testset "Phase 7 gate — whole body as one part ≡ dorsal/ventral (exact)" begin
    base_metabolic = baseline.energy_flows.metabolic_heat_flow
    # The per-part surface solve reproduces the dorsal/ventral result to solver tolerance.
    @test whole.metabolic_heat_flow ≈ base_metabolic rtol = 1e-3
    @test whole.skin_temperature ≈ baseline.thermoregulation.skin_temperature rtol = 1e-4
    @test whole.parts[1].net_metabolic ≈ baseline.energy_flows.dorsal.net_metabolic rtol = 1e-4
end

@testset "Phase 7 gate — 2 HalfCylinder + SharedCore ≈ dorsal/ventral" begin
    base_metabolic = baseline.energy_flows.metabolic_heat_flow
    mp_metabolic   = multipart.metabolic_heat_flow

    # The two parts are identical and equal-mass in the uniform symmetric case,
    # and their surface temperatures track the dorsal/ventral mean closely.
    @test multipart.parts[1].net_metabolic ≈ multipart.parts[2].net_metabolic
    @test multipart.skin_temperature ≈ baseline.thermoregulation.skin_temperature rtol = 1e-3
    @test multipart.insulation_temperature ≈ baseline.thermoregulation.insulation_temperature rtol = 1e-3

    # Metabolic heat flow now matches the dorsal/ventral result closely. Closing the
    # earlier ~13% gap needed three geometric corrections for the joined half-shape:
    #  - convective/evaporative area = total − covered flat face (solve_part_surface);
    #  - characteristic dimension supplied by the caller (the full body's), not
    #    derived from the half's reduced volume (§3.6 Tier-1 cache);
    #  - the 2π (full-circumference) fur-conductance factor in the cylindrical
    #    radiant-temperature formulas scaled by the exposed shell fraction
    #    (_shell_angle_fraction = 0.5 for a HalfCylinder).
    @test multipart.metabolic_heat_flow ≈ base_metabolic rtol = 1e-3
end
