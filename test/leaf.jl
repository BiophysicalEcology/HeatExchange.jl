using HeatExchange
using BiophysicalGeometry
using ModelParameters
using Unitful, UnitfulMoles
using Test

# Environment
air_temperature       = (273.15 + 25.0)u"K"
sky_temperature       = (273.15 + 10.0)u"K"
ground_temperature    = (273.15 + 35.0)u"K"
substrate_temperature = air_temperature
relative_humidity     = 0.40
wind_speed            = 1.0u"m/s"
atmospheric_pressure  = 101325.0u"Pa"
zenith_angle          = 30.0u"°"
global_radiation      = 800.0u"W/m^2"
substrate_conductivity = 0.5u"W/m/K"

gas_fractions = GasFractions(0.2095, 0.0003, 0.79)

environment_pars = EnvironmentalPars(;
    ground_albedo=0.15,
    ground_emissivity=0.95,
    sky_emissivity=1.0,
    elevation=0.0u"m",
    fluid=0,
    gas_fractions,
)

environment_vars = EnvironmentalVars(;
    air_temperature,
    sky_temperature,
    ground_temperature,
    substrate_temperature,
    relative_humidity,
    wind_speed,
    atmospheric_pressure,
    substrate_conductivity,
    zenith_angle,
    global_radiation,
    diffuse_fraction=0.15,
    shade=0.0,
)

environment = (; environment_pars, environment_vars)

# Leaf body: thin Plate (5 cm × 5 cm × 0.5 mm), ShortestDimension(0.7) after Gates 1980
leaf_density  = 700.0u"kg/m^3"       # typical mesophyll
leaf_mass     = 0.2u"g" |> u"kg"     # a small leaf
leaf_shape    = Plate(leaf_mass, leaf_density, 1.0, 100.0)   # b=1 (square), c=100 (thin)
leaf_body     = Body(leaf_shape, Naked())

# Traits
leaf_traits = HeatExchangeTraits(
    leaf_shape,
    InsulationParameters(),                                    # no insulation
    ExternalConductionParameters(; conduction_fraction=0.0),   # leaf not on substrate
    InternalConductionParameters(; flesh_conductivity=100.0u"W/m/K"),  # near-isothermal
    RadiationParameters(;
        body_absorptivity_dorsal=0.5,
        body_absorptivity_ventral=0.5,
        body_emissivity_dorsal=0.97,
        body_emissivity_ventral=0.97,
        sky_view_factor=0.5,
        ground_view_factor=0.5,
        vegetation_view_factor=0.0,
        bush_view_factor=0.0,
        ventral_fraction=0.5,
        solar_orientation=Intermediate(),
    ),
    ConvectionParameters(; characteristic_dimension_formula=ScaledDimension(0.7, :width_skin)),
    LeafEvaporationParameters(;
        abaxial_vapour_conductance=0.3u"mol/m^2/s",    # stomata open
        adaxial_vapour_conductance=0.0u"mol/m^2/s",
        cuticular_conductance=0.01u"mol/m^2/s",
    ),
    HydraulicParameters(; water_potential=0.0u"J/kg"),
    RespirationParameters(),
    MetabolismParameters(; metabolic_heat_flow=0.0u"W", model=PlantDarkRespiration()),
    SolveMetabolicRateOptions(),
)

leaf = Organism(leaf_body, leaf_traits)

# --- Tests ---

leaf_temperature = (273.15 + 25.0)u"K"
out = heat_balance(leaf_temperature, leaf, environment)

@test isfinite(ustrip(out.heat_balance))        # heat_balance is finite
@test isfinite(ustrip(out.surface_temperature)) # surface temperature computed

# With stomata open, transpiration should be positive under subsaturated air
@test out.mass_balance.transpiration_mass > 0.0u"g/s"

# Closing stomata (both conductances → 0) should reduce transpiration to near-cuticular only
closed_pars = LeafEvaporationParameters(;
    abaxial_vapour_conductance=0.0u"mol/m^2/s",
    adaxial_vapour_conductance=0.0u"mol/m^2/s",
    cuticular_conductance=0.01u"mol/m^2/s",
)
closed_traits = HeatExchangeTraits(
    leaf_shape,
    InsulationParameters(),
    ExternalConductionParameters(; conduction_fraction=0.0),
    InternalConductionParameters(; flesh_conductivity=100.0u"W/m/K"),
    RadiationParameters(;
        body_absorptivity_dorsal=0.5,
        body_absorptivity_ventral=0.5,
        body_emissivity_dorsal=0.97,
        body_emissivity_ventral=0.97,
        sky_view_factor=0.5,
        ground_view_factor=0.5,
        vegetation_view_factor=0.0,
        bush_view_factor=0.0,
        ventral_fraction=0.5,
        solar_orientation=Intermediate(),
    ),
    ConvectionParameters(; 
    characteristic_dimension_formula=ScaledDimension(0.7, :width_skin)), # Campbell and Norman 1998 suggests 0.7 for leaf boundary layer
    closed_pars,
    HydraulicParameters(; water_potential=0.0u"J/kg"),
    RespirationParameters(),
    MetabolismParameters(; metabolic_heat_flow=0.0u"W", model=PlantDarkRespiration()),
    SolveMetabolicRateOptions(),
)
leaf_closed = Organism(leaf_body, closed_traits)
out_closed = heat_balance(leaf_temperature, leaf_closed, environment)

@test out_closed.mass_balance.transpiration_mass < out.mass_balance.transpiration_mass

# solve_temperature should converge to a steady state
eq_out = solve_temperature(leaf, environment)
equilibrium_temperature = eq_out.core_temperature
@test isfinite(ustrip(equilibrium_temperature))
@test 273.0u"K" < equilibrium_temperature < 370.0u"K"

# At steady state, heat_balance should be near zero
@test abs(ustrip(u"W", eq_out.heat_balance)) < 1e-2   # within 0.01 W of zero

# PlantDarkRespiration: metabolic rate at 25°C should be positive and small
dark_resp = PlantDarkRespiration()
resp_rate = metabolic_rate(dark_resp, leaf_mass, leaf_temperature)
@test resp_rate > 0.0u"W"
@test resp_rate < 1.0u"W"   # a tiny leaf
