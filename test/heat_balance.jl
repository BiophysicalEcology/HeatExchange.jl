using HeatExchange
using BiophysicalGeometry
using BiologicalScaling: BasalMetabolicRate, EutherianMammal
using ModelParameters
using Unitful, UnitfulMoles
using FluidProperties
using Test

# Parameters from test/data/endoR_input_4_1.csv (Ellipsoid with fur)
# Shape 4 = Ellipsoid, furmult 1 = with insulation
shape_pars = Ellipsoid(mass = 65.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 5.0, axis_ratio_c = 5.0)
fat        = FatLayer(0.20, 901.0u"kg/m^3")

pven = 0.4
insulation_pars = InsulationParameters(;
    dorsal = FibreProperties(;
        diameter    = 3e-5u"m",
        length      = 0.0239u"m",
        density     = 3e7u"1/m^2",
        depth       = 0.0239u"m",
        reflectance = 0.2,
        conductivity = 0.209u"W/m/K",
    ),
    ventral = FibreProperties(;
        diameter    = 1e-5u"m",
        length      = 0.0139u"m",
        density     = 1e7u"1/m^2",
        depth       = 0.0139u"m",
        reflectance = 0.2,
        conductivity = 0.209u"W/m/K",
    ),
    depth_compressed       = 0.0139u"m",
    longwave_depth_fraction = 1.0,
)

# mean-weighted geometry (same as endotherm.jl test)
mean_ins_depth    = insulation_pars.dorsal.depth    * (1 - pven) + insulation_pars.ventral.depth    * pven
mean_fibre_diam   = insulation_pars.dorsal.diameter * (1 - pven) + insulation_pars.ventral.diameter * pven
mean_fibre_density = insulation_pars.dorsal.density * (1 - pven) + insulation_pars.ventral.density  * pven
fur      = FibrousLayer(mean_ins_depth, mean_fibre_diam, mean_fibre_density)
geometry = Body(shape_pars, CompositeInsulation(fur, fat))

environment_vars = EnvironmentalVars(;
    air_temperature           = u"K"(20.0u"°C"),
    reference_air_temperature = u"K"(20.0u"°C"),
    sky_temperature           = u"K"(10.0u"°C"),
    ground_temperature        = u"K"(40.0u"°C"),
    substrate_temperature     = u"K"(30.0u"°C"),
    understory_temperature    = u"K"(10.0u"°C"),
    canopy_temperature        = u"K"(20.0u"°C"),
    relative_humidity         = 0.05,
    wind_speed                = 1.0u"m/s",
    atmospheric_pressure      = 101325.0u"Pa",
    zenith_angle              = 20.0u"°",
    substrate_conductivity    = 2.79u"W/m/K",
    global_radiation          = 500.0u"W/m^2",
    diffuse_fraction          = 0.15,
    shade                     = 0.1,
)

environment_pars = EnvironmentalPars(;
    ground_albedo        = 0.2,
    ground_emissivity    = 1.0,
    sky_emissivity       = 1.0,
    elevation            = 0.0u"m",
    fluid                = Air(),
    gas_fractions        = GasFractions(0.2095, 0.000412, 0.7902),
    convection_enhancement = 1.0,
)

internal_conduction = InternalConductionParameters(;
    fat_fraction        = 0.20,
    flesh_conductivity  = 0.9u"W/m/K",
    fat_conductivity    = 0.23u"W/m/K",
    fat_density         = 901.0u"kg/m^3",
)
external_conduction = ExternalConductionParameters(; conduction_fraction = 0.3)

radiation_pars = RadiationParameters(;
    body_absorptivity_dorsal  = 0.8,
    body_absorptivity_ventral = 0.8,
    body_emissivity_dorsal    = 0.99,
    body_emissivity_ventral   = 0.99,
    sky_view_factor           = 0.5,
    ground_view_factor        = 0.5,
    ventral_fraction          = pven,
    solar_orientation         = NormalToSun(),
)

evaporation_pars = AnimalEvaporationParameters(;
    skin_wetness        = 0.0005,
    insulation_wetness  = 0.0005,
    eye_fraction        = 0.0003,
    bare_skin_fraction  = 0.0,
    insulation_fraction = 1.0,
)
respiration_pars = RespirationParameters(;
    oxygen_extraction_efficiency = 0.2,
    pant                         = 1.0,
    respiratory_quotient         = 0.8,
    exhaled_temperature_offset   = nothing,   # exhaled air at lung temperature
    exhaled_relative_humidity    = 1.0,
)
metabolism_pars = MetabolismParameters(;
    core_temperature    = u"K"(37.0u"°C"),
    metabolic_heat_flow = 77.6184u"W",
    q10                 = 2.0,
    model               = (BasalMetabolicRate(), EutherianMammal()),
)

traits = HeatExchangeTraits(
    shape_pars,
    insulation_pars,
    external_conduction,
    internal_conduction,
    radiation_pars,
    ConvectionParameters(),
    evaporation_pars,
    HydraulicParameters(),
    respiration_pars,
    metabolism_pars,
    SolveMetabolicRateOptions(; respire = true, temperature_error_tolerance = 1e-3u"K", resp_tolerance = 1e-5),
)
organism = Organism(geometry, traits)

# Solve iterative model
endo_out   = solve_metabolic_rate(organism, (; environment_pars, environment_vars), u"K"(34.0u"°C"), u"K"(29.0u"°C"))
thermoreg  = endo_out.thermoregulation
metabolic_heat_flow      = endo_out.energy_flows.metabolic_heat_flow
core_temperature     = thermoreg.core_temperature

e = (; environment_pars, environment_vars)

@testset "Each side's surface balances at the solution" begin
    for (part, side) in zip(sides(organism, e), (thermoreg.dorsal, thermoreg.ventral))
        temperatures = HeatExchange.OrganismTemperatures(core_temperature, side.skin_temperature, side.insulation_temperature)
        b = heat_balance(temperatures, metabolic_heat_flow, part; respiration_pars)
        @test abs(ustrip(u"W", b.surface_balance)) < 1e-6
        @test abs(ustrip(u"K", b.residual_skin_temperature)) < 1e-6
    end
end

# An insulated organism's temperature comes from solve_temperature.
@test_throws MethodError heat_balance(u"K"(30.0u"°C"), organism, e)

# solve_temperature returns a rich ThermoregulationOutput with whole-organism and per-side results
result = solve_temperature(organism, e)
@test result isa ThermoregulationOutput
@test :thermoregulation in propertynames(result)
@test :energy_flows in propertynames(result)
@test :dorsal in propertynames(result.energy_flows)
@test :ventral in propertynames(result.energy_flows)
@test :dorsal in propertynames(result.thermoregulation)
@test :ventral in propertynames(result.thermoregulation)
@test u"K"(0.0u"°C") < result.thermoregulation.core_temperature < u"K"(80.0u"°C")

# -------------------------------------------------------------------------
# CommonSolve interface (HeatBalanceProblem / init / solve! / reinit! / solve)
# wraps the same solve_metabolic_rate physics — results must be identical.
# -------------------------------------------------------------------------
let
    init_skin = u"K"(34.0u"°C")
    init_ins  = u"K"(29.0u"°C")
    direct = solve_metabolic_rate(organism, e, init_skin, init_ins)

    problem = HeatBalanceProblem(organism, e)
    initial_state = (; skin_temperature = init_skin, insulation_temperature = init_ins)

    # One-shot solve == direct call.
    oneshot = solve(problem; initial_state)
    @test oneshot isa ThermoregulationOutput
    @test oneshot.energy_flows.metabolic_heat_flow ≈ direct.energy_flows.metabolic_heat_flow
    @test oneshot.thermoregulation.skin_temperature ≈ direct.thermoregulation.skin_temperature
    @test oneshot.thermoregulation.core_temperature ≈ direct.thermoregulation.core_temperature

    # Stepwise: init → solve! gives the same, and threads state forward.
    solver = init(problem; initial_state)
    stepwise = solve!(solver)
    @test stepwise.energy_flows.metabolic_heat_flow ≈ direct.energy_flows.metabolic_heat_flow
    @test solver.state.skin_temperature ≈ direct.thermoregulation.skin_temperature

    # reinit! with the same problem, warm-started from the converged state,
    # re-solves consistently — but only to the solver's convergence tolerance
    # (1e-3 K), since the temperature seed differs from the cold-start `direct`.
    reinit!(solver, problem)
    resolved = solve!(solver)
    @test resolved.energy_flows.metabolic_heat_flow ≈ direct.energy_flows.metabolic_heat_flow rtol=1e-2

    # reinit! without warm start resets to the default initial state.
    reinit!(solver, problem; warm_start = false)
    @test solver.state.skin_temperature == HeatExchange.default_initial_state(problem).skin_temperature
end
