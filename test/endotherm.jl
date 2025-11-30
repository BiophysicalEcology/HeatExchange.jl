using HeatExchange
using ModelParameters
using Unitful, UnitfulMoles
using FluidProperties
using Test
using DataFrames, CSV
#using Plots

# ellipsoid model
testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

ellipsoid_input_vec = DataFrame(CSV.File("$testdir/data/ellipsoid_input.csv"))[:, 2]
ellipsoid_output_NicheMapR = DataFrame(CSV.File("$testdir/data/ellipsoid_output.csv"))[:, 2:end]

names = [
    :posture, :mass, :density, :coreT, :furdepth, :furcond, :O2eff, :stress, :windspd, :rh, :Q10, 
    :basal, :basmult
]

ellipsoid_input = (; zip(names, ellipsoid_input_vec[1:13])...)

# TODO make ellipsoid_endotherm broadcastable
air_temperatures = (ellipsoid_input_vec[14:end])u"°C"
ellipsoid_output = DataFrame(map(air_temperatures) do Ta
    ellipsoid_endotherm(;
        mrate_equation = Kleiber(),
        posture = ellipsoid_input.posture, 
        mass = (ellipsoid_input.mass)u"kg", 
        density = (ellipsoid_input.density)u"kg/m^3", 
        core_temperature = (ellipsoid_input.coreT)u"°C",
        insulation_depth = (ellipsoid_input.furdepth)u"mm", 
        k_insulation = (ellipsoid_input.furcond)u"W/m/K", 
        oxygen_extraction_efficiency = ellipsoid_input.O2eff, 
        stress_factor = ellipsoid_input.stress,
        air_temperature = Ta,      # ← element
        wind_speed = (ellipsoid_input.windspd)u"m/s", 
        relative_humidity = ellipsoid_input.rh / 100, 
        P_atmos = 101325.0u"Pa", 
        q10 = ellipsoid_input.Q10,
        metabolic_multiplier = ellipsoid_input.basmult, 
        lethal_desiccation = 0.15, 
        f_O2 = 0.2094
    )
end)

@testset "ellipsoid comparisons" begin
    @test ellipsoid_output_NicheMapR.Tskin ≈ ustrip.(u"°C", ellipsoid_output.skin_temperature) rtol = 1e-6
    @test ellipsoid_output_NicheMapR.LCT ≈ ustrip.(u"°C", ellipsoid_output.lower_critical_air_temperature) rtol = 1e-8
    @test ellipsoid_output_NicheMapR.UCT ≈ ustrip.(u"°C", ellipsoid_output.upper_critical_air_temperature) rtol = 1e-8
    @test ellipsoid_output_NicheMapR.QgenFinal ≈ ustrip.(u"W", ellipsoid_output.Q_gen_final) rtol = 1e-7
    @test ellipsoid_output_NicheMapR.Qresp_W ≈ ustrip.(u"W", ellipsoid_output.Q_respiration) rtol = 1e-4
    @test ellipsoid_output_NicheMapR.H2Oloss_W ≈ ustrip.(u"W", ellipsoid_output.Q_evap) rtol = 1e-7
    @test ellipsoid_output_NicheMapR.mlO2ph ≈ ustrip.(u"ml/hr", ellipsoid_output.O2_consumption_rate) rtol = 1e-7
    @test ellipsoid_output_NicheMapR.H2O_gph ≈ ustrip.(u"g/hr", ellipsoid_output.total_water_loss_rate) rtol = 1e-3
    @test ellipsoid_output_NicheMapR.PctBasal / 100 ≈ ellipsoid_output.basal_metabolic_rate_fraction rtol = 1e-7
    @test ellipsoid_output_NicheMapR.massph_percent / 100 ≈ ustrip.(u"hr^-1", ellipsoid_output.fractional_mass_loss) rtol = 1e-3
end

# solvendo

T_air = u"K"(20.0u"°C") # air temperature at local height
T_air_reference = T_air # air temperature at reference height
T_substrate = T_air # ground temperature
T_sky = T_air # sky temperature
T_ground = T_air # ground temperature
T_vegetation = T_air_reference
wind_speed = 0.1u"m/s" # wind speed
rh = 0.05 # relative humidity (fractional)
solar_radiation = 0.0u"W/m^2" # solar radiation, horizontal plane
fraction_diffuse = 0.15 # fraction of solar radiation that is diffuse
diffuse_radiation = solar_radiation * fraction_diffuse
direct_radiation = solar_radiation - diffuse_radiation
zenith_angle = 20u"°" # zenith angle of sun (degrees from overhead)
elevation = 0.0u"m" # elevation
α_ground = 0.8 # solar absorptivity of substrate, fractional
ϵ_ground = 1.0 # substrate emissivity, fractional
ϵ_sky = 1.0 # sky emissivity, fractional
shade = 0.0 # shade, fractional

# other environmental variables
fluid = 0 # fluid type: 0 = air; 1 = fresh water; 2 = salt water
T_conduction = T_substrate # surface temperature for conduction
k_substrate = 2.79u"W/m/K" # substrate thermal conductivity
T_bush = T_air # bush temperature
P_atmos = atmospheric_pressure(elevation) # Pa, negative means elevation is used
fO2 = 0.2095 # oxygen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (fractional)
fN2 = 0.7902 # nitrogen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (fractional)
fCO2 = 0.000412 # carbon dioxide concentration of air, to account for non-atmospheric concentrations e.g. in burrows (fractional)
p_diffuse = 0.15 # proportion of solar radiation that is diffuse 

# BEHAVIOUR

shade = 0 # shade level (%)
shape_b_step = 0.1 # allows the animal to uncurl to shape_b_max, the value being the increment shape_b is increased per iteration
T_core_step = 0.1u"K" # turns on core temperature elevation, the value being the increment by which T_core is increased per iteration
skin_wetness_step = 0.001 # turns on sweating, the value being the increment by which PCTWET is increased per iteration
skin_wetness_max = 1.0 # maximum surface area that can be wet (%)
k_flesh_step = 0.1u"W/m/K" # turns on thermal conductivity increase (W/mK), the value being the increment by which k_flesh is increased per iteration
k_flesh_max = 2.8u"W/m/K" # maximum flesh conductivity (W/mK)
pant = 1 # multiplier on breathing rate to simulate panting (-)
pant_step = 0.1 # increment for multiplier on breathing rate to simulate panting (-)
pant_multiplier = 1.05 # multiplier on basal metabolic rate at maximum panting level (-)
insulation_step = 1.0

# MORPHOLOGY

# geometry
mass = 65.0u"kg" # kg
ρ_flesh = 1000.0u"kg/m^3" # kg/m3
ρ_fat = 901.0u"kg/m^3" # kg/m3
fat_fraction = 0.0 # proportion body fat (subcutaneous)
shape_b = 1.1 # current ratio between long and short axis, must be > 1 (-)
shape_b_max = 5.0 # max possible ratio between long and short axis, must be > 1 (-)
shape_c = shape_b # current ratio of length:height (plate)
ventral_fraction = 0.5 # fraction of surface area that is ventral fur (fractional, 0-1)
conduction_fraction = 0.0 # fraction of surface area that is touching the substrate (fractional, 0-1)
#SAMODE = 0 # if 0, uses surface area for SHAPE parameter geometry, if 1, uses bird skin surface area allometry from Walsberg & King. 1978. JEB 76:185–189, if 2 uses mammal surface area from Stahl 1967.J. App. Physiol. 22, 453–460.
#ORIENT = 0 # if 1 = normal to sun's rays (heat maximising), if 2 = parallel to sun's rays (heat minimising), 3 = vertical and changing with solar altitude, or 0 = average

# fur properties
insulation_conductivity_dorsal = nothing # user-specified fur thermal conductivity (W/mK), not used if 0
insulation_conductivity_ventral = nothing # user-specified fur thermal conductivity (W/mK), not used if 0
fibre_diameter_dorsal = 30.0u"μm" # hair diameter, dorsal (m)
fibre_diameter_ventral = 30.0u"μm" # hair diameter, ventral (m)
fibre_length_dorsal = 23.9u"mm" # hair length, dorsal (m)
fibre_length_ventral = 23.9u"mm" # hair length, ventral (m)
insulation_depth_dorsal = 2.0u"mm" # fur depth, dorsal (m)
insulation_depth_ventral = 2.0u"mm" # fur depth, ventral (m)
max_insulation_depth_dorsal = fibre_length_dorsal # max fur depth, dorsal (m)
max_insulation_depth_ventral = fibre_length_ventral # max fur depth, ventral (m)fibre_density_dorsal = 3000E+04 # hair density, dorsal (1/m2)
fibre_density_dorsal = 3000u"cm^-2" # hair density, dorsal (1/m2)
fibre_density_ventral = 3000u"cm^-2" # hair density, ventral (1/m2)
insulation_reflectance_dorsal = 0.301  # fur reflectivity dorsal (fractional, 0-1)
insulation_reflectance_ventral = 0.301  # fur reflectivity ventral (fractional, 0-1)
insulation_depth_compressed = insulation_depth_ventral # depth of compressed fur (for conduction) (m)
fibre_conductivity = 0.209u"W/m/K" # hair thermal conductivity (W/m°C)
longwave_depth_fraction = 1 # fractional depth of fur at which longwave radiation is exchanged (0-1)

# radiation exchange
ϵ_body_dorsal = 0.99 # animal emissivity (-)
ϵ_body_ventral = 0.99 # animal emissivity (-)
F_bush = 0.0 # this is for veg below/around animal (at TALOC)
F_ground = 0.5 # reference configuration factor to ground
F_sky = 0.5 # configuration factor to sky
F_vegetation = 0 # configuration factor to sky

# PHYSIOLOGY

# thermal
T_core_target = u"K"(37.0u"°C") # core temperature (°C)
T_core_max = u"K"(39.0u"°C") # maximum core temperature (°C)
k_flesh = 0.9u"W/m/K" # initial thermal conductivity of flesh (0.412 - 2.8 W/m°C)
k_fat = 0.230u"W/m/K" # conductivity of fat (W/mK)

# evaporation
skin_wetness = 0.005 # part of the skin surface that is wet (fractional)
insulation_wetness = 0.0 # part of the fur/feathers that is wet after rain (fractional)
bare_skin_fraction = 0.0 # surface area for evaporation that is skin, e.g. licking paws (fractional)
eye_fraction = 0.0 # surface area made up by the eye (fractional) - make zero if sleeping
Δ_breath = u"K"(0.0u"°C") # offset between air temperature and breath (°C)
rh_exit = 1.0 # relative humidity of exhaled air, fractional
ψ_org = 0.0u"J/kg"

# metabolism/respiration
Q_minimum = (70 * ustrip(u"kg", mass)^0.75) * (4.185 / (24 * 3.6))u"W" # basal heat generation (W) from Kleiber (1947)
rq = 0.80 # respiratory quotient (fractional, 0-1)
fO2_extract = 0.2 # O2 extraction efficiency (fractional)
pant_max = 5 # maximum breathing rate multiplier to simulate panting (-)
fur_step = 1 # # incremental fractional reduction in insulation_depth from piloerect state (-) (a value greater than zero triggers piloerection response)
q10 = 2 # q10 factor for adjusting BMR for T_core
T_core_min = 19 # minimum core temperature during torpor (TORPOR = 1)
tolerance_torpor = 0.05 # allowable tolerance of heat balance as a fraction of torpid metabolic rate

# initial conditions
T_skin = T_core_target - 3u"K" # skin temperature (°C)
T_insulation = T_air # fur/air interface temperature (°C)

# other model settings
convection_enhancement = 1 # convective enhancement factor for turbulent conditions, typically 1.4
tolerance = 0.001u"K" # tolerance for SIMULSOL
thermoregulate = true # invoke thermoregulatory response
respire = true # compute respiration and associated heat loss
thermoregulation_mode = 1 # 1 = raise core then pant then sweat, 2 = raise core and pant simultaneously, then sweat
torpor = false # go into torpor if possible (drop T_core down to TC_MIN)

bodyshape = Ellipsoid(mass, ρ_flesh, shape_b, shape_b) # define shape

environmental_vars = EnvironmentalVars(
    T_air,
    T_air_reference,
    T_sky,
    T_ground,
    T_substrate,
    T_bush,
    T_vegetation,
    rh,
    wind_speed,
    P_atmos,
    zenith_angle,
    k_substrate,
    solar_radiation,
    direct_radiation,
    diffuse_radiation,
)

environmental_pars = EnvironmentalPars(;
    α_ground,
    ϵ_ground,
    ϵ_sky,
    elevation,
    fluid,
    fN2,
    fO2,
    fCO2,
    convection_enhancement,
    )

body_pars = BodyPars(;
    mass,
    ρ_flesh,
    ρ_fat,
    fat_fraction,
    shape_b,
    shape_c,
)

integument_pars = IntegumentPars(;
    α_body_dorsal = 1 - insulation_reflectance_dorsal,
    α_body_ventral = 1 - insulation_reflectance_ventral,
    ϵ_body_dorsal,
    ϵ_body_ventral,
    F_sky,
    F_ground,
    F_vegetation,
    F_bush,
    eye_fraction,
    skin_wetness,
    bare_skin_fraction,
    conduction_fraction,
    ventral_fraction,
)

physio_pars = PhysioPars(;
    Q_minimum,
    q10,
    k_flesh,
    k_fat,
    fO2_extract,
    rq,
    Δ_breath,
    rh_exit,
)

thermoreg_pars = ThermoregulationPars(;
    insulation_step,
    shape_b_step,
    shape_b_max,
    T_core_target,
    T_core_max,
    T_core_min,
    T_core_step,
    k_flesh_step,
    k_flesh_max,
    pant,
    pant_step,
    pant_multiplier,
    pant_max,
    skin_wetness_step,
    skin_wetness_max,
    )

insulation_pars = InsulationPars(;
    insulation_conductivity_dorsal,
    insulation_conductivity_ventral,
    fibre_diameter_dorsal,
    fibre_diameter_ventral,
    fibre_length_dorsal,
    fibre_length_ventral,
    insulation_depth_dorsal,
    insulation_depth_ventral,
    fibre_density_dorsal,
    fibre_density_ventral,
    insulation_reflectance_dorsal,
    insulation_reflectance_ventral,
    insulation_depth_compressed,
    fibre_conductivity,
    longwave_depth_fraction,
    insulation_wetness,
)

organism_vars = OrganismalVars(;
    T_core = T_core_target,
    T_skin,
    T_insulation,
    T_lung = T_core_target,
    ψ_org,
)

model_pars = EndoModelPars(
    thermoregulation_mode,
    thermoregulate,
    respire,
    torpor,
    tolerance
    )

endotherm_out = endotherm(; model_pars, bodyshape, body_pars, integument_pars, insulation_pars, 
    physio_pars, thermoreg_pars, environmental_pars, organism_vars, environmental_vars)

endotherm_out.thermoregulation







T_air = u"K"(20.0u"°C") # air temperature at local height
T_air_reference = T_air # air temperature at reference height
T_substrate = T_air # ground temperature
T_sky = T_air # sky temperature
T_ground = T_air # ground temperature
T_vegetation = T_air_reference
wind_speed = 0.1u"m/s" # wind speed
rh = 0.05 # relative humidity (fractional)
solar_radiation = 0.0u"W/m^2" # solar radiation, horizontal plane
fraction_diffuse = 0.15 # fraction of solar radiation that is diffuse
diffuse_radiation = solar_radiation * fraction_diffuse
direct_radiation = solar_radiation - diffuse_radiation
zenith_angle = 20u"°" # zenith angle of sun (degrees from overhead)
elevation = 0.0u"m" # elevation
α_ground = 0.8 # solar absorptivity of substrate, fractional
ϵ_ground = 1.0 # substrate emissivity, fractional
ϵ_sky = 1.0 # sky emissivity, fractional
shade = 0.0 # shade, fractional

# other environmental variables
fluid = 0 # fluid type: 0 = air; 1 = fresh water; 2 = salt water
T_conduction = T_substrate # surface temperature for conduction
k_substrate = 2.79u"W/m/K" # substrate thermal conductivity
T_bush = T_air # bush temperature
P_atmos = atmospheric_pressure(elevation) # Pa, negative means elevation is used
fO2 = 0.2095 # oxygen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (fractional)
fN2 = 0.7902 # nitrogen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (fractional)
fCO2 = 0.000412 # carbon dioxide concentration of air, to account for non-atmospheric concentrations e.g. in burrows (fractional)
p_diffuse = 0.15 # proportion of solar radiation that is diffuse 

# BEHAVIOUR

shade = 0 # shade level (%)
shape_b_step = 0.1 # allows the animal to uncurl to shape_b_max, the value being the increment shape_b is increased per iteration
T_core_step = 0.1u"K" # turns on core temperature elevation, the value being the increment by which T_core is increased per iteration
skin_wetness_step = 0.001 # turns on sweating, the value being the increment by which PCTWET is increased per iteration
skin_wetness_max = 1.0 # maximum surface area that can be wet (%)
k_flesh_step = 0.1u"W/m/K" # turns on thermal conductivity increase (W/mK), the value being the increment by which k_flesh is increased per iteration
k_flesh_max = 2.8u"W/m/K" # maximum flesh conductivity (W/mK)
pant = 1 # multiplier on breathing rate to simulate panting (-)
pant_step = 0.1 # increment for multiplier on breathing rate to simulate panting (-)
pant_multiplier = 1.05 # multiplier on basal metabolic rate at maximum panting level (-)
insulation_step = 1.0

# MORPHOLOGY

# geometry
mass = 65.0u"kg" # kg
ρ_flesh = 1000.0u"kg/m^3" # kg/m3
ρ_fat = 901.0u"kg/m^3" # kg/m3
fat_fraction = 0.0 # proportion body fat (subcutaneous)
shape_b = 1.1 # current ratio between long and short axis, must be > 1 (-)
shape_b_max = 5.0 # max possible ratio between long and short axis, must be > 1 (-)
shape_c = shape_b # current ratio of length:height (plate)
ventral_fraction = 0.5 # fraction of surface area that is ventral fur (fractional, 0-1)
conduction_fraction = 0.0 # fraction of surface area that is touching the substrate (fractional, 0-1)
#SAMODE = 0 # if 0, uses surface area for SHAPE parameter geometry, if 1, uses bird skin surface area allometry from Walsberg & King. 1978. JEB 76:185–189, if 2 uses mammal surface area from Stahl 1967.J. App. Physiol. 22, 453–460.
#ORIENT = 0 # if 1 = normal to sun's rays (heat maximising), if 2 = parallel to sun's rays (heat minimising), 3 = vertical and changing with solar altitude, or 0 = average

# fur properties
insulation_conductivity_dorsal = nothing # user-specified fur thermal conductivity (W/mK), not used if 0
insulation_conductivity_ventral = nothing # user-specified fur thermal conductivity (W/mK), not used if 0
fibre_diameter_dorsal = 30.0u"μm" # hair diameter, dorsal (m)
fibre_diameter_ventral = 30.0u"μm" # hair diameter, ventral (m)
fibre_length_dorsal = 23.9u"mm" # hair length, dorsal (m)
fibre_length_ventral = 23.9u"mm" # hair length, ventral (m)
insulation_depth_dorsal = 2.0u"mm" # fur depth, dorsal (m)
insulation_depth_ventral = 2.0u"mm" # fur depth, ventral (m)
max_insulation_depth_dorsal = fibre_length_dorsal # max fur depth, dorsal (m)
max_insulation_depth_ventral = fibre_length_ventral # max fur depth, ventral (m)fibre_density_dorsal = 3000E+04 # hair density, dorsal (1/m2)
fibre_density_dorsal = 3000u"cm^-2" # hair density, dorsal (1/m2)
fibre_density_ventral = 3000u"cm^-2" # hair density, ventral (1/m2)
insulation_reflectance_dorsal = 0.301  # fur reflectivity dorsal (fractional, 0-1)
insulation_reflectance_ventral = 0.301  # fur reflectivity ventral (fractional, 0-1)
insulation_depth_compressed = insulation_depth_ventral # depth of compressed fur (for conduction) (m)
fibre_conductivity = 0.209u"W/m/K" # hair thermal conductivity (W/m°C)
longwave_depth_fraction = 1 # fractional depth of fur at which longwave radiation is exchanged (0-1)

# radiation exchange
ϵ_body_dorsal = 0.99 # animal emissivity (-)
ϵ_body_ventral = 0.99 # animal emissivity (-)
F_bush = 0.0 # this is for veg below/around animal (at TALOC)
F_ground = 0.5 # reference configuration factor to ground
F_sky = 0.5 # configuration factor to sky
F_vegetation = 0 # configuration factor to sky

# PHYSIOLOGY

# thermal
T_core_target = u"K"(37.0u"°C") # core temperature (°C)
T_core_max = u"K"(39.0u"°C") # maximum core temperature (°C)
k_flesh = 0.9u"W/m/K" # initial thermal conductivity of flesh (0.412 - 2.8 W/m°C)
k_fat = 0.230u"W/m/K" # conductivity of fat (W/mK)

# evaporation
skin_wetness = 0.005 # part of the skin surface that is wet (fractional)
insulation_wetness = 0.0 # part of the fur/feathers that is wet after rain (fractional)
bare_skin_fraction = 0.0 # surface area for evaporation that is skin, e.g. licking paws (fractional)
eye_fraction = 0.0 # surface area made up by the eye (fractional) - make zero if sleeping
Δ_breath = u"K"(0.0u"°C") # offset between air temperature and breath (°C)
rh_exit = 1.0 # relative humidity of exhaled air, fractional
ψ_org = 0.0u"J/kg"

# metabolism/respiration
Q_minimum = (70 * ustrip(u"kg", mass)^0.75) * (4.185 / (24 * 3.6))u"W" # basal heat generation (W) from Kleiber (1947)
rq = 0.80 # respiratory quotient (fractional, 0-1)
fO2_extract = 0.2 # O2 extraction efficiency (fractional)
pant_max = 5 # maximum breathing rate multiplier to simulate panting (-)
fur_step = 1 # # incremental fractional reduction in insulation_depth from piloerect state (-) (a value greater than zero triggers piloerection response)
q10 = 2 # q10 factor for adjusting BMR for T_core
T_core_min = 19 # minimum core temperature during torpor (TORPOR = 1)
tolerance_torpor = 0.05 # allowable tolerance of heat balance as a fraction of torpid metabolic rate

# initial conditions
T_skin = T_core_target - 3u"K" # skin temperature (°C)
T_insulation = T_air # fur/air interface temperature (°C)

# other model settings
convection_enhancement = 1 # convective enhancement factor for turbulent conditions, typically 1.4
tolerance = 0.001u"K" # tolerance for SIMULSOL
thermoregulate = true # invoke thermoregulatory response
respire = true # compute respiration and associated heat loss
thermoregulation_mode = 1 # 1 = raise core then pant then sweat, 2 = raise core and pant simultaneously, then sweat
torpor = false # go into torpor if possible (drop T_core down to TC_MIN)
