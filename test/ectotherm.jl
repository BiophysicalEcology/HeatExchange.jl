using HeatExchange
using BiophysicalGeometry
using FluidProperties
using ModelParameters
using Unitful, UnitfulMoles
using Roots
using Test
using DataFrames, CSV

testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

ecto_input_vec = DataFrame(CSV.File("$testdir/data/ectoR_input.csv"))[:, 2]
ecto_output_vec = DataFrame(CSV.File("$testdir/data/ectoR_output.csv"))[:, 2]

names = [
    :Ww_g,
    :alpha,
    :epsilon,
    :rho_body,
    :fatosk,
    :fatosb,
    :shape,
    :shape_a,
    :shape_b,
    :shape_c,
    :conv_enhance,
    :custom_shape1,
    :custom_shape2,
    :custom_shape3,
    :custom_shape4,
    :custom_shape5,
    :custom_shape6,
    :custom_shape7,
    :custom_shape8,
    :pct_cond,
    :pct_touch,
    :postur,
    :orient,
    :k_flesh,
    :M_1,
    :M_2,
    :M_3,
    :M_4,
    :Q_act,
    :pct_wet,
    :pct_eyes,
    :pct_mouth,
    :psi_body,
    :pantmax,
    :F_O2,
    :RQ,
    :delta_air,
    :leaf,
    :g_vs_ab,
    :g_vs_ad,
    :elevation,
    :alpha_sub,
    :epsilon_sub,
    :epsilon_sky,
    :pres,
    :fluid,
    :O2gas,
    :CO2gas,
    :N2gas,
    :K_sub,
    :PDIF,
    :SHADE,
    :QSOLR,
    :Z,
    :TA,
    :TGRD,
    :TSUBST,
    :TSKY,
    :VEL,
    :RH,
]

ecto_input = (; zip(names, ecto_input_vec)...)

names = [
    :TC,
    :TSKIN,
    :TLUNG,
    :QSOL,
    :QIRIN,
    :QMET,
    :QRESP,
    :QEVAP,
    :QIROUT,
    :QCONV,
    :QCOND,
    :ENB,
    :O2_ml,
    :H2OResp_g,
    :H2OCut_g,
    :H2OEyes_g,
    :AREA,
    :AV,
    :AT,
    :AL,
    :VOL,
    :R,
    :R1,
    :ASEMAJR,
    :BSEMINR,
    :CSEMINR,
    :ASILN,
    :ASILP,
    :AEFF,
    :QSOLAR,
    :QSDIFF,
    :QSRSB,
    :QSSKY,
    :QSOBJ,
    :QRESP2,
    :GEVAP,
    :AIRML1,
    :AIRML2,
    :WMOL1,
    :WMOL2,
    :O2MOL1,
    :O2MOL2,
    :CO2MOL1,
    :CO2MOL2,
    :CONVAR,
    :QCONV2,
    :HC,
    :HD,
    :SH,
    :QFREE,
    :HCFREE,
    :HCFORC,
    :SHFREE,
    :SHFORC,
    :HDFREE,
    :HDFORC,
    :QSEVAP,
    :WEVAP,
    :WRESP,
    :WCUT,
    :WEYES,
    :G_V,
]

ecto_output = (; zip(names, ecto_output_vec)...)

# environment
air_temperature = u"K"((ecto_input.TA)u"°C")
sky_temperature = u"K"((ecto_input.TSKY)u"°C")
substrate_temperature = u"K"((ecto_input.TSUBST)u"°C")
ground_temperature = substrate_temperature
substrate_conductivity = (ecto_input.K_sub)u"W/m/K"
atmospheric_pressure = (ecto_input.pres)u"Pa"
relative_humidity = (ecto_input.RH/100)
elevation = (ecto_input.elevation)u"m"
wind_speed = (ecto_input.VEL)u"m/s"
fluid = ecto_input.fluid
zenith_angle = (ecto_input.Z)u"°"
global_radiation = (ecto_input.QSOLR)u"W/m^2"
ground_albedo = ecto_input.alpha_sub
shade = ecto_input.SHADE / 100
ground_emissivity = ecto_input.epsilon_sub
sky_emissivity = ecto_input.epsilon_sky
fO2 = ecto_input.O2gas / 100
fCO2 = ecto_input.CO2gas / 100
fN2 = ecto_input.N2gas / 100
gas_fractions = GasFractions(fO2, fCO2, fN2)
convection_enhancement = ecto_input.conv_enhance

# organism morphology
mass = u"kg"((ecto_input.Ww_g)u"g")
ρ_flesh = (ecto_input.rho_body)u"kg/m^3"
shape_b = ecto_input.shape_b
shape_c = ecto_input.shape_c
eye_fraction = ecto_input.pct_eyes / 100
conduction_fraction = ecto_input.pct_cond / 100

sky_view_factor = ecto_input.fatosk
ground_view_factor = ecto_input.fatosb
body_absorptivity_dorsal = ecto_input.alpha
body_absorptivity_ventral = ecto_input.alpha
body_emissivity_dorsal = ecto_input.epsilon
body_emissivity_ventral = ecto_input.epsilon
Le = 0.025u"m"
#conduction_depth::D =    Param(2.5u"cm") this should be an environmental variable

ventral_fraction = 0.5

# organism physiology
flesh_conductivity = (ecto_input.k_flesh)u"W/m/K"
ψ_org = (ecto_input.psi_body)u"J/kg"
skin_wetness = ecto_input.pct_wet / 100
oxygen_extraction_efficiency = ecto_input.F_O2 / 100
respiratory_quotient = ecto_input.RQ
mass_normalisation = ecto_input.M_1
mass_exponent = ecto_input.M_2
thermal_sensitivity = ecto_input.M_3
metabolic_state = ecto_input.M_4

# organism behaviour
if ecto_input.postur == 0.0
    solar_orientation = Intermediate()
elseif ecto_input.postur == 1.0
    solar_orientation = NormalToSun()
else
    solar_orientation = ParallelToSun()
end
pant = 1.0

# geometry
shape_pars = DesertIguana(mass, ρ_flesh)
#shape_pars = Cylinder(mass, ρ_flesh, shape_b)
geometry = Body(shape_pars, Naked())
A_total = total_area(geometry)
A_sil_normal = silhouette_area(shape_pars, NormalToSun())
A_sil_parallel = silhouette_area(shape_pars, ParallelToSun())
A_silhouette = silhouette_area(shape_pars, solar_orientation)
A_conduction = A_total * conduction_fraction
A_convection = A_total * (1 - conduction_fraction)
A_eff = A_convection * skin_wetness

# geometry test
@test A_total ≈ (ecto_output.AREA)u"m^2" rtol=1e-9
@test A_conduction ≈ (ecto_output.AV)u"m^2" rtol=1e-9
@test A_sil_normal ≈ (ecto_output.ASILN)u"m^2" rtol=1e-9
@test A_sil_parallel ≈ (ecto_output.ASILP)u"m^2" rtol=1e-9
@test A_eff ≈ (ecto_output.AEFF)u"m^2" rtol=1e-9
@test geometry.geometry.characteristic_dimension ≈ (ecto_output.AL)u"m" rtol=1e-9
@test geometry.geometry.length[1] ≈ (ecto_output.R1)u"m" rtol=1e-9
@test geometry.geometry.volume ≈ (ecto_output.VOL)u"m^3" rtol=1e-9

# organism state
core_temperature = u"K"((ecto_output.TC)u"°C")

# calculate heat fluxes

# metabolism
metabolic_heat_flow = metabolic_rate(
    AndrewsPough2(; mass_normalisation, mass_exponent, thermal_sensitivity, metabolic_state, O2conversion=Typical()), mass, core_temperature
)

# metabolism test
@test metabolic_heat_flow ≈ (ecto_output.QMET)u"W" rtol=1e-9

# respiration
resp_out = respiration(
    MetabolicRates(; metabolic=metabolic_heat_flow),
    RespirationParameters(; oxygen_extraction_efficiency, pant, respiratory_quotient),
    AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure),
    mass,
    core_temperature,  # lung_temperature
    air_temperature;
    gas_fractions,
)
respiration_heat_flow = resp_out.respiration_heat_flow

# respiration test
molar_fluxes_in = resp_out.molar_fluxes_in
molar_fluxes_out = resp_out.molar_fluxes_out
@test u"mol/s"(molar_fluxes_in.carbon_dioxide) ≈ (ecto_output.CO2MOL1)u"mol/s" rtol=1e-6
@test u"mol/s"(molar_fluxes_out.carbon_dioxide) ≈ (ecto_output.CO2MOL2)u"mol/s" rtol=1e-6
@test u"mol/s"(molar_fluxes_in.water) ≈ (ecto_output.WMOL1)u"mol/s" rtol=1e-6
@test u"mol/s"(molar_fluxes_out.water) ≈ (ecto_output.WMOL2)u"mol/s" rtol=1e-6
@test u"mol/s"(molar_fluxes_in.oxygen) ≈ (ecto_output.O2MOL1)u"mol/s" rtol=1e-6
@test u"mol/s"(molar_fluxes_out.oxygen) ≈ (ecto_output.O2MOL2)u"mol/s" rtol=1e-6
@test u"mol/s"(molar_fluxes_in.air) ≈ (ecto_output.AIRML1)u"mol/s" rtol=1e-6
@test u"mol/s"(molar_fluxes_out.air) ≈ (ecto_output.AIRML2)u"mol/s" rtol=1e-6
@test resp_out.respiration_heat_flow ≈ (ecto_output.QRESP2)u"W" rtol=1e-3
@test resp_out.respiration_mass ≈ (ecto_output.GEVAP)u"g/s" rtol=1e-3

respiration_heat_flow = resp_out.respiration_heat_flow
net_generated_heat_flow = metabolic_heat_flow - respiration_heat_flow
specific_generated_heat_flow = net_generated_heat_flow / geometry.geometry.volume

# compute skin and lung temperature
(; surface_temperature, lung_temperature) = surface_and_lung_temperature(; body=geometry, flesh_conductivity, specific_metabolic_heat_production=specific_generated_heat_flow, core_temperature)
skin_temperature = surface_temperature  # alias for later use

# test lung temperature
@test lung_temperature ≈ u"K"((ecto_output.TLUNG)u"°C") rtol=1e-7

# solar radiation
diffuse_fraction = ecto_input.PDIF
absorptivities = Absorptivities(;
    body=DorsalVentral(body_absorptivity_dorsal, body_absorptivity_ventral),
    ground=ground_albedo,
)
view_factors = ViewFactors(sky_view_factor, ground_view_factor, 0.0, 0.0)
solar_conditions = SolarConditions(;
    zenith_angle,
    global_radiation,
    diffuse_fraction,
    shade,
)
solar_out = solar(geometry, absorptivities, view_factors, solar_conditions, A_silhouette, A_conduction)
solar_flow = solar_out.solar_flow

# longwave radiation
emissivities = Emissivities(;
    body=DorsalVentral(body_emissivity_dorsal, body_emissivity_ventral),
    ground=ground_emissivity,
    sky=sky_emissivity,
)
env_temps = EnvironmentTemperatures(air_temperature, sky_temperature, ground_temperature, air_temperature, air_temperature, substrate_temperature)
ir_gain = radiation_in(geometry, view_factors, emissivities, env_temps; conduction_fraction)
longwave_flow_in = ir_gain.longwave_flow_in
ir_loss = radiation_out(geometry, view_factors, emissivities, conduction_fraction, skin_temperature, skin_temperature)
longwave_flow_out = ir_loss.longwave_flow_out

# conduction
conduction_flow = conduction(; conduction_area=A_conduction, L=Le, surface_temperature=skin_temperature, substrate_temperature, substrate_conductivity)

# convection
conv_out = convection(;
    body=geometry,
    area=A_convection,
    air_temperature,
    surface_temperature=skin_temperature,
    wind_speed,
    atmospheric_pressure,
    fluid,
    gas_fractions,
    convection_enhancement,
)
convection_heat_flow = conv_out.heat_flow

# evaporation
evap_pars_test = AnimalEvaporationParameters(; skin_wetness, eye_fraction, bare_skin_fraction=1.0)
atmos_test = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
evap_out = evaporation(evap_pars_test, conv_out.mass, atmos_test, A_convection, skin_temperature, air_temperature; water_potential=ψ_org, gas_fractions)
evaporation_heat_flow = evap_out.evaporation_heat_flow

# energy balance test
@test solar_flow ≈ (ecto_output.QSOL)u"W" rtol=1e-9
@test conduction_flow ≈ (ecto_output.QCOND)u"W" rtol=1e-5
@test convection_heat_flow ≈ (ecto_output.QCONV)u"W" rtol=1e-4
@test evaporation_heat_flow ≈ (ecto_output.QEVAP)u"W" rtol=1e-4

heat_flow_in = solar_flow + longwave_flow_in + metabolic_heat_flow
heat_flow_out = longwave_flow_out + convection_heat_flow + evaporation_heat_flow + respiration_heat_flow + conduction_flow
heat_flow_in - heat_flow_out

# This is broken
# (For it to work we either need an anonymous function 
# or a functor struct - you can make a struct into a function!
# then we have all the parameters attached to it.
#@test_broken T_core_s = find_zero(ectotherm, (T_air - 40K, T_air + 100K), Bisection())
#@test_broken T_core_C = (Unitful.ustrip(T_core_s) - 273.15)°C

# using structs to pass parameters

conduction_pars_external = ExternalConductionParameters(; conduction_fraction)

conduction_pars_internal = InternalConductionParameters(; flesh_conductivity)

radiation_pars = RadiationParameters(;
    body_absorptivity_dorsal,
    body_absorptivity_ventral,
    body_emissivity_dorsal,
    body_emissivity_ventral,
    silhouette_area=A_silhouette,
    total_area=A_total,
    conduction_area=A_conduction,
    sky_view_factor,
    ground_view_factor,
    ventral_fraction,
    solar_orientation,
)

convection_pars = ConvectionParameters(; convection_area=A_convection)

evaporation_pars = AnimalEvaporationParameters(; skin_wetness, eye_fraction)

hydraulic_pars = HydraulicParameters(; water_potential=ψ_org)

respiration_pars = RespirationParameters(; oxygen_extraction_efficiency, pant, respiratory_quotient)

metabolism_pars = MetabolismParameters(; metabolic_heat_flow, model=AndrewsPough2())

traits = HeatExchangeTraits(
    shape_pars,
    InsulationParameters(),
    conduction_pars_external,
    conduction_pars_internal,
    radiation_pars,
    convection_pars,
    evaporation_pars,
    hydraulic_pars,
    respiration_pars,
    metabolism_pars,
    SolveMetabolicRateOptions(),
)
# construct the Model which holds the parameters of the organism in the Organism concrete struct, of type AbstractOrganism
#lizard = Model(Organism(geometry, traits))
lizard = Organism(geometry, traits)

# get the environmental parameters
environment_pars = EnvironmentalPars(; ground_albedo, ground_emissivity, sky_emissivity, elevation, fluid, gas_fractions)

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
    diffuse_fraction,
    shade,
)

environment = (; environment_pars, environment_vars)
# define the method 'heat_balance' for passing to find_zero, which dispatches off 'lizard'
air_temperature = environment_vars.air_temperature
heat_balance(air_temperature, lizard, environment)

core_temperature_s = find_zero(
    t -> heat_balance(t, lizard, environment), (air_temperature - 40u"K", air_temperature + 100u"K"), Bisection()
)
core_temperature_C = (Unitful.ustrip(core_temperature_s) - 273.15)u"°C"

# test core temperature calculation
@test core_temperature_C ≈ (ecto_output.TC)u"°C" rtol=1e-4

heat_balance_out = heat_balance(core_temperature_s, lizard, environment)

@test heat_balance_out.core_temperature ≈ (ecto_output.TC + 273.15)u"K" rtol=1e-5
@test heat_balance_out.surface_temperature ≈ (ecto_output.TSKIN + 273.15)u"K" rtol=1e-5
@test heat_balance_out.lung_temperature ≈ (ecto_output.TLUNG + 273.15)u"K" rtol=1e-5

@test heat_balance_out.energy_balance.solar_flow ≈ (ecto_output.QSOL)u"W" rtol=1e-9
@test heat_balance_out.energy_balance.longwave_flow_in ≈ (ecto_output.QIRIN)u"W" rtol=1e-3 # TODO make better?
@test heat_balance_out.energy_balance.longwave_flow_out ≈ (ecto_output.QIROUT)u"W" rtol=1e-3 # TODO make better?
@test heat_balance_out.energy_balance.conduction_flow ≈ (ecto_output.QCOND)u"W" rtol=1e-3 # TODO make better?
@test heat_balance_out.energy_balance.convection_heat_flow ≈ (ecto_output.QCONV)u"W" rtol=1e-4
@test heat_balance_out.energy_balance.evaporation_heat_flow ≈ (ecto_output.QEVAP)u"W" rtol=1e-4
@test heat_balance_out.energy_balance.metabolic_heat_flow ≈ (ecto_output.QMET)u"W" rtol=1e-4 # TODO make better?
@test heat_balance_out.energy_balance.respiration_heat_flow ≈ (ecto_output.QRESP)u"W" rtol=1e-3 # TODO make better?

@test heat_balance_out.mass_balance.oxygen_flow ≈ (ecto_output.O2_ml)u"ml/hr" rtol=1e-4 # TODO make better?
@test heat_balance_out.mass_balance.cutaneous_mass ≈ (ecto_output.H2OCut_g / 3600)u"g/s" rtol=1e-4 # TODO make better?
@test heat_balance_out.mass_balance.eye_mass ≈ (ecto_output.H2OEyes_g / 3600)u"g/s" rtol=1e-4 # TODO make better?

molar_fluxes_in = heat_balance_out.respiration_out.molar_fluxes_in
molar_fluxes_out = heat_balance_out.respiration_out.molar_fluxes_out
@test u"mol/s"(molar_fluxes_in.carbon_dioxide) ≈ (ecto_output.CO2MOL1)u"mol/s" rtol=1e-4
@test u"mol/s"(molar_fluxes_out.carbon_dioxide) ≈ (ecto_output.CO2MOL2)u"mol/s" rtol=1e-4
@test u"mol/s"(molar_fluxes_in.water) ≈ (ecto_output.WMOL1)u"mol/s" rtol=1e-4
@test u"mol/s"(molar_fluxes_out.water) ≈ (ecto_output.WMOL2)u"mol/s" rtol=1e-3
@test u"mol/s"(molar_fluxes_in.oxygen) ≈ (ecto_output.O2MOL1)u"mol/s" rtol=1e-4
@test u"mol/s"(molar_fluxes_out.oxygen) ≈ (ecto_output.O2MOL2)u"mol/s" rtol=1e-4
@test u"mol/s"(molar_fluxes_in.air) ≈ (ecto_output.AIRML1)u"mol/s" rtol=1e-4
@test u"mol/s"(molar_fluxes_out.air) ≈ (ecto_output.AIRML2)u"mol/s" rtol=1e-4
@test heat_balance_out.respiration_out.respiration_heat_flow ≈ (ecto_output.QRESP2)u"W" rtol=1e-3
@test heat_balance_out.respiration_out.respiration_mass ≈ (ecto_output.GEVAP)u"g/s" rtol=1e-3

@test heat_balance_out.evaporation_out.m_cut ≈ (ecto_output.WCUT)u"g/s" rtol=1e-4 # TODO check if this can be better
@test heat_balance_out.evaporation_out.m_eyes ≈ (ecto_output.WEYES)u"g/s" rtol=1e-4 # TODO check if this can be better
@test heat_balance_out.evaporation_out.evaporation_heat_flow ≈ (ecto_output.QEVAP)u"W" rtol=1e-4 # TODO check if this can be better
