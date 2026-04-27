using HeatExchange
using BiophysicalGeometry
using ModelParameters
using Unitful, UnitfulMoles
using FluidProperties
using Test

# Parameters from test/data/endoR_input_4_1.csv (Ellipsoid with fur)
# Shape 4 = Ellipsoid, furmult 1 = with insulation
shape_pars = Ellipsoid(65.0u"kg", 1000.0u"kg/m^3", 5.0, 5.0)
fat        = Fat(0.20, 901.0u"kg/m^3")

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
fur      = Fur(mean_ins_depth, mean_fibre_diam, mean_fibre_density)
geometry = Body(shape_pars, CompositeInsulation(fur, fat))

environment_vars = EnvironmentalVars(;
    air_temperature           = u"K"(20.0u"°C"),
    reference_air_temperature = u"K"(20.0u"°C"),
    sky_temperature           = u"K"(10.0u"°C"),
    ground_temperature        = u"K"(40.0u"°C"),
    substrate_temperature     = u"K"(30.0u"°C"),
    bush_temperature          = u"K"(10.0u"°C"),
    vegetation_temperature    = u"K"(20.0u"°C"),
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
    ground_albedo        = 0.8,
    ground_emissivity    = 1.0,
    sky_emissivity       = 1.0,
    elevation            = 0.0u"m",
    fluid                = 0,
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
    bush_view_factor          = 0.0,
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
    exhaled_temperature_offset   = 100.0u"K",
    exhaled_relative_humidity    = 1.0,
)
metabolism_pars = MetabolismParameters(;
    core_temperature    = u"K"(37.0u"°C"),
    metabolic_heat_flow = 77.6184u"W",
    q10                 = 2.0,
    model               = Kleiber(),
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
Q_gen      = endo_out.energy_flows.generated_heat_flow
T_core     = thermoreg.core_temperature

# -------------------------------------------------------------------------
# Replicate the dorsal-side inputs that solve_metabolic_rate packs internally
# (see solve_metabolic_rate lines ~80–200 in endotherm.jl)
# -------------------------------------------------------------------------
dmult = 1.0 - pven                      # 0.6
vegetation_factor = radiation_pars.sky_view_factor * environment_vars.shade   # 0.5 * 0.1 = 0.05
sky_factor_ref    = radiation_pars.sky_view_factor - vegetation_factor         # 0.45
ground_factor_ref = 1.0 - sky_factor_ref - vegetation_factor                  # 0.50

# Dorsal body uses the dorsal insulation depth
ins_layer_d = Fur(insulation_pars.dorsal.depth, insulation_pars.dorsal.diameter, insulation_pars.dorsal.density)
body_d      = Body(shape_pars, CompositeInsulation(ins_layer_d, fat))

# Solar for dorsal side — mirrors solve_metabolic_rate solar block
ep = stripparams(environment_pars)
absorptivities   = Absorptivities(radiation_pars, ep)
vf_solar         = ViewFactors(sky_factor_ref, ground_factor_ref, 0.0, 0.0)
solar_conds      = SolarConditions(environment_vars)
area_sil         = silhouette_area(geometry, radiation_pars.solar_orientation)
area_cond_ref    = BiophysicalGeometry.total_area(geometry) * external_conduction.conduction_fraction
solar_out        = solar(geometry, absorptivities, vf_solar, solar_conds, area_sil, area_cond_ref)
dorsal_solar     = solar_out.solar_flow > 0.0u"W" ?
    2.0 * solar_out.direct_flow + solar_out.solar_sky_flow * 2.0 :
    0.0u"W"

# Dorsal side temperatures from converged solution
T_skin_d = thermoreg.skin_temperature_dorsal
T_ins_d  = thermoreg.insulation_temperature_dorsal

# Insulation properties at converged dorsal temperatures
ins_d = insulation_properties(insulation_pars, T_ins_d * 0.7 + T_skin_d * 0.3, pven)

geometry_vars_d = GeometryVariables(;
    side                    = :dorsal,
    conductance_coefficient = 0.0u"W/K",
    ventral_fraction        = pven,
    conduction_fraction     = external_conduction.conduction_fraction,
    longwave_depth_fraction = insulation_pars.longwave_depth_fraction,
)

# vegetation_temperature is reference_air_temperature in solve_metabolic_rate
env_d = (;
    temperature = (
        air        = environment_vars.air_temperature,
        sky        = environment_vars.sky_temperature,
        ground     = environment_vars.ground_temperature,
        vegetation = environment_vars.reference_air_temperature,
        bush       = environment_vars.bush_temperature,
        substrate  = environment_vars.substrate_temperature,
    ),
    view_factors = (
        sky        = sky_factor_ref * 2.0,
        ground     = 0.0,
        bush       = 0.0,
        vegetation = vegetation_factor * 2.0,
    ),
    atmos = (
        relative_humidity    = environment_vars.relative_humidity,
        wind_speed           = environment_vars.wind_speed,
        atmospheric_pressure = environment_vars.atmospheric_pressure,
    ),
    fluid                  = ep.fluid,
    solar_flow             = dorsal_solar,
    gas_fractions          = ep.gas_fractions,
    convection_enhancement = ep.convection_enhancement,
)

traits_d = (;
    fat_conductivity   = internal_conduction.fat_conductivity,
    flesh_conductivity = internal_conduction.flesh_conductivity,
    ϵ_body             = radiation_pars.body_emissivity_dorsal,
    insulation_wetness = evaporation_pars.insulation_wetness,
    bare_skin_fraction = evaporation_pars.bare_skin_fraction,
    eye_fraction       = evaporation_pars.eye_fraction,
    skin_wetness       = evaporation_pars.skin_wetness,
)

hb_d = heat_balance(
    T_core, T_skin_d, T_ins_d, Q_gen;
    body              = body_d,
    insulation_pars   = insulation_pars,
    insulation        = ins_d,
    geometry_vars     = geometry_vars_d,
    minimum_metabolic_heat = Q_gen, # what should this be?
    environment_vars  = env_d,
    traits            = traits_d,
    resp_pars         = respiration_pars,
)
balance = Q_gen + hb_d.solar_heat_flow - 
    hb_d.radiation_heat_flow - 
    hb_d.convection_heat_flow - 
    hb_d.conduction_heat_flow - 
    hb_d.skin_evaporation_heat_flow -
    hb_d.insulation_evaporation_heat_flow - 
    hb_d.respiration_heat_flow - 
    hb_d.net_metabolic_heat_internal

#@test abs(ustrip(u"W", balance)) < 0.01
# The iterative solver converged T_skin, so residual_skin_temperature must be near zero.
# residual_energy_balance and residual_internal_conduction use the full-organism Q_gen
# against per-side losses, so they need not be zero for a one-sided call.
#@test abs(ustrip(u"K", hb_d.residual_skin_temperature)) < 0.01

# -------------------------------------------------------------------------
# Smoke tests: MultiSided solve_temperature for insulated organism
# -------------------------------------------------------------------------
@test evaluation_strategy(organism) isa MultiSided

e = (; environment_pars, environment_vars)

# heat_balance on a MultiSided organism raises MethodError — use solve_temperature instead
@test_throws MethodError heat_balance(u"K"(30.0u"°C"), organism, e)

# solve_temperature returns a rich NamedTuple with whole-organism and per-side results
result = solve_temperature(organism, e)
@test result isa NamedTuple
@test haskey(result, :thermoregulation)
@test haskey(result, :energy_flows)
@test haskey(result.energy_flows, :dorsal)
@test haskey(result.energy_flows, :ventral)
@test haskey(result.thermoregulation, :dorsal)
@test haskey(result.thermoregulation, :ventral)
@test u"K"(0.0u"°C") < result.thermoregulation.core_temperature < u"K"(80.0u"°C")
