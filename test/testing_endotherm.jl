using HeatExchange
using ModelParameters
using Unitful, UnitfulMoles
using FluidProperties
#using Plots

ellipsoid_endotherm(; 
    posture=4.5, 
    mass=0.5u"kg", 
    density=1000.0u"kg/m^3", 
    core_temperature=37u"°C",
    fur_depth=5u"mm", 
    fur_conductivity=0.04u"W/m/K", 
    oxygen_extraction_efficiency=0.2, 
    stress_factor=0.6,
    air_temperature=20.0u"°C", 
    wind_speed=1.0u"m/s", 
    relative_humidity=0.5, 
    P_atmos = 101325.0u"Pa", 
    Q10=3,
    minimum_metabolic_rate=missing, 
    metabolic_multiplier=1, 
    lethal_desiccation=0.15, 
    f_O2=0.2094)

insulation = InsulationPars(;
    fibre_diameter_dorsal=30.0u"μm", # hair diameter, dorsal (m)
    fibre_diameter_ventral=30.0u"μm", # hair diameter, ventral (m)
    fibre_length_dorsal=23.9u"mm", # hair length, dorsal (m)
    fibre_length_ventral=23.9u"mm", # hair length, ventral (m)
    insulation_depth_dorsal=9.0u"mm", # fur depth, dorsal (m)
    insulation_depth_ventral=9.0u"mm", # fur depth, ventral (m)
    fibre_density_dorsal=3000u"cm^-2", # hair density, dorsal (1/m2)
    fibre_density_ventral=3000u"cm^-2", # hair density, ventral (1/m2)
    insulation_reflectance_dorsal=0.301,  # fur reflectivity dorsal (fractional, 0-1)
    insulation_reflectance_ventral=0.301,  # fur reflectivity ventral (fractional, 0-1)
    insulation_depth_compressed=9.0u"mm", # depth of compressed fur (for conduction) (m)
    fibre_conductivity=0.209u"W/m/K", # hair thermal conductivity (W/m°C)
    longwave_depth_fraction=1.0,
)

insulation_out = insulation_properties(; insulation, insulation_temperature=(273.15 + 20.0)u"K", ventral_fraction=0.3)

density = 1000.0u"kg/m^3"
trunkmass = 0.04u"kg"
trunkshapeb = 2
trunkshape = Cylinder(trunkmass, density, trunkshapeb) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = (fat, fur)

trunk = Body(trunkshape, Naked()) # construct a Body, which is furred
trunk = Body(trunkshape, fur) # construct a Body, which is furred
trunk = Body(trunkshape, fat) # construct a Body, which has a fat layer
trunk = Body(trunkshape, CompositeInsulation(fur, fat))

trunk_geometry = geometry(trunkshape, fur)

θ = 90u"°"
trunksilhouette = silhouette_area(trunk, θ)

trunkshape = Sphere(trunkmass, density) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = CompositeInsulation(fat, fur)

trunk = Body(trunkshape, fur) # construct a Body, which is furred
trunk = Body(trunkshape, fat) # construct a Body, which has a fat layer
trunk = Body(trunkshape, CompositeInsulation(fur, fat))


trunkshape = Plate(trunkmass, density, trunkshapeb, trunkshapeb) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = (fat, fur)

trunk = Body(trunkshape, Naked()) # construct a Body, which is furred
trunk = Body(trunkshape, fur) # construct a Body, which is furred
trunk = Body(trunkshape, fat) # construct a Body, which has a fat layer
trunk = Body(trunkshape, CompositeInsulation(fur, fat))

trunkshape = Ellipsoid(trunkmass, density, trunkshapeb, trunkshapeb) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = (fat, fur)

trunk = Body(trunkshape, Naked()) # construct a Body, which is furred
trunk = Body(trunkshape, fur) # construct a Body, which is furred
trunk = Body(trunkshape, fat) # construct a Body, which has a fat layer
trunk = Body(trunkshape, CompositeInsulation(fur, fat))

#### solvendo ####

T_air = u"K"(20.0u"°C") # air temperature at local height
T_air_reference = T_air # air temperature at reference height
T_substrate = T_air # ground temperature
T_sky = T_air # sky temperature
wind_speed = 0.1u"m/s" # wind speed
rh = 0.05 # relative humidity (fractional)
solar_radiation = 0.0u"W/m^2" # solar radiation, horizontal plane
fraction_diffuse = 0.15 # fraction of solar radiation that is diffuse
diffuse_radiation = solar_radiation * fraction_diffuse
direct_radiation = solar_radiation - diffuse_radiation
zenith_angle = 20u"°" # zenith angle of sun (degrees from overhead)
elevation = 0.0u"m" # elevation
α_substrate = 0.8 # solar absorptivity of substrate, fractional
ϵ_substrate = 1.0 # substrate emissivity, fractional
ϵ_sky = 1.0 # sky emissivity, fractional
shade = 0.0 # shade, fractional
fluid = 0

# other environmental variables
fluid_type = 0 # fluid type: 0 = air; 1 = fresh water; 2 = salt water
T_conduction = T_substrate # surface temperature for conduction
k_substrate = 2.79u"W/m/K" # substrate thermal conductivity
T_bush = T_air # bush temperature
P_atmos = atmospheric_pressure(elevation) # Pa, negative means elevation is used
fO2 = 0.2095 # oxygen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (fractional)
fN2 = 0.7902 # nitrogen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (fractional)
fCO2 = 0.000412 # carbon dioxide concentration of air, to account for non-atmospheric concentrations e.g. in burrows (fractional)
p_diffuse = 0.15 # proportion of solar radiation that is diffuse 
g = u"gn" # acceleration due to gravity

# BEHAVIOUR

shade = 0 # shade level (%)
uncurl = 0.1 # allows the animal to uncurl to SHAPE_B_MAX, the value being the increment SHAPE_B is increased per iteration
Tc_inc = 0.1 # turns on core temperature elevation, the value being the increment by which T_core is increased per iteration
skin_wetness_step = 0.001 # turns on sweating, the value being the increment by which PCTWET is increased per iteration
skin_wetness_max = 1.0 # maximum surface area that can be wet (%)
flesh_conductivity_step = 0.1u"W/m/K" # turns on thermal conductivity increase (W/mK), the value being the increment by which k_flesh is increased per iteration
flesh_conductivity_max = 2.8u"W/m/K" # maximum flesh conductivity (W/mK)
pant = 1 # multiplier on breathing rate to simulate panting (-)
pant_step = 0.1 # increment for multiplier on breathing rate to simulate panting (-)
pant_multiplier = 1.05 # multiplier on basal metabolic rate at maximum panting level (-)

# MORPHOLOGY

# geometry
mass = 65u"kg" # kg
ρ_body = 1000u"kg/m^3" # kg/m3
ρ_fat = 901u"kg/m^3" # kg/m3
subcutaneous_fat = false # is subcutaneous fat present? (0 is no, 1 is yes)
fat_fraction = 0.20 # proportion body fat
shape_b = 1.1 # current ratio between long and short axis, must be > 1 (-)
shape_b_max = 5 # max possible ratio between long and short axis, must be > 1 (-)
shape_c = shape_b # current ratio of length:height (plate)
ventral_fraction = 0.5 # fraction of surface area that is ventral fur (fractional, 0-1)
conduction_fraction = 0 # fraction of surface area that is touching the substrate (fractional, 0-1)
#SAMODE = 0 # if 0, uses surface area for SHAPE parameter geometry, if 1, uses bird skin surface area allometry from Walsberg & King. 1978. JEB 76:185–189, if 2 uses mammal surface area from Stahl 1967.J. App. Physiol. 22, 453–460.
#ORIENT = 0 # if 1 = normal to sun's rays (heat maximising), if 2 = parallel to sun's rays (heat minimising), 3 = vertical and changing with solar altitude, or 0 = average

# fur properties
#insulation_conductivity = 0 # user-specified fur thermal conductivity (W/mK), not used if 0
fibre_diameter_dorsal = 30.0u"μm" # hair diameter, dorsal (m)
fibre_diameter_ventral = 30.0u"μm" # hair diameter, ventral (m)
fibre_length_dorsal = 23.9u"mm" # hair length, dorsal (m)
fibre_length_ventral = 23.9u"mm" # hair length, ventral (m)
insulation_depth_dorsal = 9.0u"mm" # fur depth, dorsal (m)
insulation_depth_ventral = 9.0u"mm" # fur depth, ventral (m)
max_insulation_depth_dorsal = fibre_length_dorsal # max fur depth, dorsal (m)
max_insulation_depth_ventral = fibre_length_ventral # max fur depth, ventral (m)fibre_density_dorsal = 3000E+04 # hair density, dorsal (1/m2)
fibre_density_dorsal = 3000u"cm^-2" # hair density, dorsal (1/m2)
fibre_density_ventral = 3000u"cm^-2" # hair density, ventral (1/m2)
insulation_reflectance_dorsal = 0.301  # fur reflectivity dorsal (fractional, 0-1)
insulation_reflectance_ventral = 0.301  # fur reflectivity ventral (fractional, 0-1)
insulation_depth_compressed = insulation_depth_ventral # depth of compressed fur (for conduction) (m)
fibre_conductivity = 0.209u"W/m/K" # hair thermal conductivity (W/m°C)
ventral_fraction = 0.3
conduction_fraction = 0.0
longwave_depth_fraction = 1 # fractional depth of fur at which longwave radiation is exchanged (0-1)

# radiation exchange
ϵ_body = 0.99 # animal emissivity (-)
F_bush = 0 # this is for veg below/around animal (at TALOC)
F_ground_reference = 0.5 # reference configuration factor to ground
F_sky_reference = 0.5 # configuration factor to sky

# PHYSIOLOGY

# thermal
T_core = u"K"(37.0u"°C") # core temperature (°C)
T_core_maximum = u"K"(39.0u"°C") # maximum core temperature (°C)
flesh_conductivity = 0.9u"W/m/K" # initial thermal conductivity of flesh (0.412 - 2.8 W/m°C)
fat_conductivity = 0.230u"W/m/K" # conductivity of fat (W/mK)
insulation_conductivity = nothing # conductivity of insulation (W/mK)

# evaporation
skin_wetness = 0.005 # part of the skin surface that is wet (fractional)
fur_wetness = 0.0 # part of the fur/feathers that is wet after rain (fractional)
bare_skin_fraction = 0.0 # surface area for evaporation that is skin, e.g. licking paws (fractional)
eye_fraction = 0.0 # surface area made up by the eye (fractional) - make zero if sleeping
Δ_breath = u"K"(0.0u"°C") # offset between air temperature and breath (°C)
rh_breath = 1.0 # relative humidity of exhaled air, fractional
ψ_org = 0.0u"J/kg"

# metabolism/respiration
Q_minimum = (70 * ustrip(u"kg", mass)^0.75) * (4.185 / (24 * 3.6))u"W" # basal heat generation (W) from Kleiber (1947)
respiratory_quotient = 0.80 # respiratory quotient (fractional, 0-1)
oxygen_extraction_efficiency = 0.2 # O2 extraction efficiency (fractional)
pant_max = 5 # maximum breathing rate multiplier to simulate panting (-)
#AIRVOL_MAX = 1e12 # maximum absolute breathing rate to simulate panting (L/s), can override PANT_MAX
fur_step = 1 # # incremental fractional reduction in insulation_depth from piloerect state (-) (a value greater than zero triggers piloerection response)
Q10 = 2 # Q10 factor for adjusting BMR for T_core
T_core_minimum = 19 # minimum core temperature during torpor (TORPOR = 1)
tolerance_torpor = 0.05 # allowable tolerance of heat balance as a fraction of torpid metabolic rate

# initial conditions
T_skin = T_core - 3u"K" # skin temperature (°C)
T_insulation = T_air # fur/air interface temperature (°C)

# other model settings
convection_enhancement = 1 # convective enhancement factor for turbulent conditions, typically 1.4
tolerance_simusol = 0.001 # tolerance for SIMULSOL
tolerance_zbrent = 1e-5 # tolerance for ZBRENT
thermoregulate = true # invoke thermoregulatory response
respire = true # compute respiration and associated heat loss
thermoregulation_mode = 1 # 1 = raise core then pant then sweat, 2 = raise core and pant simultaneously, then sweat
torpor = false # go into torpor if possible (drop T_core down to TC_MIN)

# insulation properties
insulation = InsulationPars(;
    insulation_temperature = T_insulation * 0.7 + T_skin * 0.3,
    fibre_diameter_dorsal, # hair diameter, dorsal (m)
    fibre_diameter_ventral, # hair diameter, ventral (m)
    fibre_length_dorsal, # hair length, dorsal (m)
    fibre_length_ventral, # hair length, ventral (m)
    insulation_depth_dorsal, # fur depth, dorsal (m)
    insulation_depth_ventral, # fur depth, ventral (m)
    fibre_density_dorsal, # hair density, dorsal (1/m2)
    fibre_density_ventral, # hair density, ventral (1/m2)
    insulation_reflectance_dorsal,  # fur reflectivity dorsal (fractional, 0-1)
    insulation_reflectance_ventral,  # fur reflectivity ventral (fractional, 0-1)
    insulation_depth_compressed, # depth of compressed fur (for conduction) (m)
    ventral_fraction,
    fibre_conductivity, # hair thermal conductivity (W/m°C)
)
insulation_out = insulation_properties(insulation)
effective_conductivities = insulation_out.effective_conductivities # effective thermal conductivity of fur array, mean, dorsal, ventral (W/mK)
absorption_coefficients = insulation_out.absorption_coefficients # term involved in computing optical thickess (1/mK2)
optical_thickness_factors = insulation_out.optical_thickness_factors # optical thickness array, mean, dorsal, ventral (m)
fibre_diameters = insulation_out.fibre_diameters # fur diameter array, mean, dorsal, ventral (m)
fibre_lengths = insulation_out.fibre_lengths # fur length array, mean, dorsal, ventral (m)
fibre_densities = insulation_out.fibre_densities # fur density array, mean, dorsal, ventral (1/m2)
insulation_depths = insulation_out.insulation_depths # fur depth array, mean, dorsal, ventral (m)
insulation_reflectance = insulation_out.insulation_reflectance # fur reflectivity array, mean, dorsal, ventral (fractional, 0-1)
insulation_test = insulation_out.insulation_test # test of presence of fur (length x diamater x density x depth) (-)
insulation_conductivity_compressed = insulation_out.insulation_conductivity_compressed # effective thermal conductivity of compressed ventral fur (W/mK)

# geometry
bodyshape = Ellipsoid(mass, density, shape_b, shape_b) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_depths[1], fibre_diameters[1], fibre_densities[1])
composite_insulation = (fat, fur)

geometry_out = Body(bodyshape, CompositeInsulation(fur, fat))

area_mammal_fur = mammal_fur_area(geometry_out)
area_mammal_skin = mammal_skin_area(bodyshape)
area_bird_plumage = bird_plumage_area(bodyshape)
area_bird_skin = bird_skin_area(bodyshape)

volume = geometry_out.geometry.volume # volume, m3
characteristic_dimension = geometry_out.geometry.characteristic_dimension # characteristic dimension for convection, m
area_total = geometry_out.geometry.area.total # total area, m2
area_silhouette = silhouette_area(bodyshape, geometry_out, 0u"°")
area_skin = geometry_out.geometry.area.skin # area of skin, m2
area_skin_evaporation = geometry_out.geometry.area.convection # area of skin for convection/evaporation (total skin area - hair area), m2
area_convection = area_total * (1 - ventral_fraction) # area of skin for convection/evaporation (total skin area - hair area), m2
area_conduction = area_total * conduction_fraction # area of skin for convection/evaporation (total skin area - hair area), m2
fat_thickness = geometry_out.geometry.lengths.fat
r1 = geometry_out.geometry.lengths[2] / 2 # shape-specific core-skin radius in shortest dimension, m
r2 = geometry_out.geometry.lengths[2] / 2 + insulation_depths[1] # shape-specific core-fur/feather interface radius in shortest dimension, m


# solar radiation normal to sun's rays
normal_radiation = zenith_angle < 90u"°" ? solar_radiation / cos(zenith_angle) : solar_radiation

α_body_dorsal = 1 - insulation_reflectance_dorsal #solar absorptivity of dorsal fur (fractional, 0-1)
α_body_ventral = 1 - insulation_reflectance_ventral # solar absorptivity of ventral fur (fractional, 0-1)

# correct F_sky for vegetaion overhead
F_vegetation = F_sky_reference * shade
F_sky = F_sky_reference - F_vegetation
F_ground = F_ground_reference

Q_solar, Q_direct, Q_solar_sky, Q_solar_substrate = solar_out = solar(α_body_dorsal, α_body_ventral, area_silhouette, area_total, area_conduction, F_ground, F_sky, α_substrate, shade, normal_radiation, direct_radiation, diffuse_radiation)
Q_dorsal = Q_direct + Q_solar_sky
Q_ventral = Q_solar_substrate
Q_diffuse = Q_solar_sky + Q_solar_substrate


(; Q_conv, hc, hd, Sh,
   Q_free, Q_forc,
   hc_free, hc_forc,
   Sh_free, Sh_forc,
   hd_free, hd_forc) = convection(geometry_out, area_convection, T_air, T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)


environmental_vars = EnvironmentalVars(
    T_air,
    T_sky,
    T_substrate,
    rh,
    wind_speed,
    P_atmos,
    zenith_angle,
    k_substrate,
    solar_radiation,
    direct_radiation,
    diffuse_radiation,
)

environmental_pars = EnvironmentalPars(
    α_substrate,
    shade,
    ϵ_substrate,
    ϵ_sky,
    elevation,
    fluid,
    fN2,
    fO2,
    fCO2,
)

# function SIMULSOL(tolerance_simusol::Float64, geometry_out, FURVARS::Vector{Float64},
# GEOMVARS::Vector{Float64}, ENVVARS::Vector{Float64},
# TRAITS::Vector{Float64}, T_insulation::Float64, skin_wetness::Float64,
# T_skin::Float64)

# Constants
σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)


# Unpack FURVARS
fibre_length, insulation_depth, insulation_conductivity, k_eff = FURVARS[1:4]
BETARA = FURVARS[5:7]
insulation_test, ZL, LHAIR, DHAIR, RHO, REFL, KHAIR, S = FURVARS[8:15]
S = Int(S)
s = 1
insulation_depth = insulation_depths[s + 1]
fibre_length = fibre_lengths[s + 1]
absorption_coefficient = absorption_coefficients[s + 1]

# Unpack GEOMVARS
SHAPE, SUBQFAT, SURFAR, volume, D, area_convection, CONVSK, r_insulation, r_flesh, r_skin, longwave_depth_fraction, RRAD, \
ASEMAJ, BSEMIN, CSEMIN, CD, conduction_fraction, r_compressed, BLCMP, k_compressed, CONV_ENHANCE = GEOMVARS

# Unpack ENVVARS
FLTYPE, T_air, T_skin, T_bush, T_vegetation, TLOWER, T_sky, T_conduction, RH, wind_speed, BP, ALT, \
F_sky, F_bush, F_vegetation, F_ground, Q_sol, GRAV = ENVVARS

# Unpack TRAITS
T_core, k_flesh, k_fat, ϵ_body, FATTHK, FLYHR, FURWET, PCTBAREVAP, PCTEYES = TRAITS

# Initialize
SOLPRO = 1.0
SOLCT = 0.0
NTRY = 0.0
success = false
RESULTS = zeros(16)

if insulation_test > 0
    # Furred body
    while true
        NTRY += 1
        for I in 1:20
            # Convective heat transfer
            (; hc, hd, hd_free) = convection(geometry_out, area_convection, T_air, T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)

                        
            # Evaporative heat loss
            SEVAPRES = zeros(7)
            SEVAP_ENDO(BP, T_air, RH, wind_speed, T_core, T_skin, ALT, skin_wetness, FLYHR,
                CONVSK, hd, hd_free, PCTBAREVAP, PCTEYES, insulation_depth, FURWET,
                T_insulation, area_convection, SEVAPRES)
            Q_evap_skin = SEVAPRES[1]
            Q_evap_insulation = SEVAPRES[7]
            evap_out = evaporation(T_skin, ψ_org, skin_wetness, area_convection, hd, eye_fraction, T_air, rh, P_atmos, fO2, fCO2, fN2)


            # Radiation properties
            (; effective_conductivities, absorption_coefficients) = insulation_out = insulation_properties(;
                insulation_temperature=T_insulation * 0.7 + T_skin * 0.3,
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
                ventral_fraction,
                fibre_conductivity,
            )
            k_eff = effective_conductivities[s + 1]

            # Effective fur conductivity
            if !isnothing(insulation_conductivity)
                k_insulation = insulation_conductivity
            else
                T_rad_approx = T_skin * (1 - longwave_depth_fraction) + T_insulation * longwave_depth_fraction
                k_rad = (16 * σ * T_rad_approx^3) / (3 * absorption_coefficients[1])
                k_insulation = k_eff + k_rad
            end
            T_radiant = radiant_temperature(geometry_out, insulation_out)
            # Geometry-specific calculations
            if Int(IPT) == 1
                # Cylinder geometry
                compression_fraction = (conduction_fraction * 2 * π * k_compressed * fibre_length) / log(r_compressed / r_skin)
                T_ins_compressed = conduction_fraction > 0 ? (compression_fraction * T_skin + CD * T_conduction) / (CD + compression_fraction) : 0.0
                cd1 = (k_compressed / log(r_compressed / r_skin)) * conduction_fraction + (k_insulation / log(r_insulation / r_skin)) * (1 - conduction_fraction)
                cd2 = (k_compressed / log(r_compressed / r_skin)) * conduction_fraction
                cd3 = (k_insulation / log(r_insulation / r_skin)) * (1 - conduction_fraction)
                dv1 = 1 + ((2 * π * fibre_length * r_flesh^2 * cd1) / (4 * k_flesh * volume)) + ((2 * π * fibre_length * r_flesh^2 * cd1) / (2 * k_fat * volume)) * log(r_skin / r_flesh)
                dv2 = Q_evap_skin * ((r_flesh^2 * cd1) / (4 * k_flesh * volume)) + Q_evap_skin * ((r_flesh^2 * cd1) / (2 * k_fat * volume)) * log(r_skin / r_flesh)
                dv3 = ((2 * π * fibre_length) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) * r_flesh^2 / (2 * volume)
                dv4 = longwave_depth_fraction < 1 ? cd2 + (k_insulation / log(r_insulation / RRAD)) * (1 - conduction_fraction) : 1.0
                T_radiant = longwave_depth_fraction < 1 ? dv3 / dv4 + (T_ins_compressed * cd2) / dv4 + (T_insulation * ((k_insulation / log(r_insulation / RRAD)) * (1 - conduction_fraction))) / dv4 : T_insulation
            elseif Int(IPT) == 2
                # Sphere geometry
                compression_fraction = (conduction_fraction * 4 * π * k_compressed * r_compressed * r_skin) / (r_compressed - r_skin)
                T_ins_compressed = conduction_fraction > 0 ? (compression_fraction * T_skin + CD * T_conduction) / (CD + compression_fraction) : 0.0
                cd1 = ((k_compressed * r_compressed) / (r_compressed - r_skin)) * conduction_fraction + ((k_insulation * r_insulation) / (r_insulation - r_skin)) * (1 - conduction_fraction)
                cd2 = ((k_compressed * r_compressed) / (r_compressed - r_skin)) * conduction_fraction
                cd3 = ((k_insulation * r_insulation) / (r_insulation - r_skin)) * (1 - conduction_fraction)
                dv1 = 1 + ((4 * π * r_skin * r_flesh^2 * cd1) / (6 * k_flesh * volume)) + ((4 * π * r_skin * r_flesh^3 * cd1) / (3 * k_fat * volume)) * ((r_skin - r_flesh) / (r_flesh * r_skin))
                dv2 = Q_evap_skin * ((r_flesh^2 * cd1) / (6 * k_flesh * volume)) + Q_evap_skin * ((r_flesh^3 * cd1) / (3 * k_fat * volume)) * ((r_skin - r_flesh) / (r_flesh * r_skin))
                dv3 = ((4 * π * r_skin) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) * r_flesh^3 / (3 * volume * RRAD)
                dv4 = longwave_depth_fraction < 1 ? cd2 + ((k_insulation * r_insulation) / (r_insulation - RRAD)) * (1 - conduction_fraction) : 1.0
                T_radiant = longwave_depth_fraction < 1 ? dv3 / dv4 + (T_ins_compressed * cd2) / dv4 + (T_insulation * ((k_insulation * r_insulation) / (r_insulation - RRAD) * (1 - conduction_fraction))) / dv4 : T_insulation
            else
                # Ellipsoid geometry
                FLSHASEMAJ = ASEMAJ - FATTHK
                FLSHBSEMIN = BSEMIN - FATTHK
                FLSHCSEMIN = CSEMIN - FATTHK
                ASQG = (Int(SUBQFAT) == 1 && FATTHK > 0) ? FLSHASEMAJ^2 : ASEMAJ^2
                BSQG = (Int(SUBQFAT) == 1 && FATTHK > 0) ? FLSHBSEMIN^2 : BSEMIN^2
                CSQG = (Int(SUBQFAT) == 1 && FATTHK > 0) ? FLSHCSEMIN^2 : CSEMIN^2
                SSQG = (ASQG * BSQG * CSQG) / (ASQG * BSQG + ASQG * CSQG + BSQG * CSQG)
                BG = (Int(SUBQFAT) == 1 && FATTHK > 0) ? FLSHBSEMIN : BSEMIN
                BS = BSEMIN
                BL = BSEMIN + ZL
                BR = BS + longwave_depth_fraction * ZL
                compression_fraction = (conduction_fraction * 3 * k_compressed * volume * BLCMP * BS) / ((sqrt(3 * SSQG))^3 * (BLCMP - BS))
                T_ins_compressed = conduction_fraction > 0 ? (compression_fraction * T_skin + CD * T_conduction) / (CD + compression_fraction) : 0.0
                cd1 = ((k_compressed * BLCMP) / (BLCMP - BS)) * conduction_fraction + ((k_insulation * BL) / (BL - BS)) * (1 - conduction_fraction)
                cd2 = ((k_compressed * BLCMP) / (BLCMP - BS)) * conduction_fraction
                cd3 = ((k_insulation * BL) / (BL - BS)) * (1 - conduction_fraction)
                dv1 = 1 + (3 * BS * SSQG * cd1) / (2 * k_flesh * (sqrt(3 * SSQG)^3)) + (BS * cd1) / k_fat * ((BS - BG) / (BS * BG))
                dv2 = Q_evap_skin * ((SSQG * cd1) / (2 * k_flesh * volume)) + Q_evap_skin * ((sqrt(3 * SSQG)^3 * cd1) / (3 * k_fat * volume)) * ((BS - BG) / (BS * BG))
                dv3 = (BS / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) / BR
                dv4 = longwave_depth_fraction < 1 ? cd2 + ((k_insulation * BL) / (BL - BR)) * (1 - conduction_fraction) : 1.0
                T_radiant = longwave_depth_fraction < 1 ? dv3 / dv4 + (T_ins_compressed * cd2) / dv4 + (T_insulation * ((k_insulation * BL) / (BL - BR) * (1 - conduction_fraction))) / dv4 : T_insulation
            end

            # Radiative heat fluxes
            Q_rad1 = area_convection * F_sky * 4 * ϵ_body * σ * ((T_radiant + T_sky) / 2)^3
            Q_rad2 = area_convection * F_bush * 4 * ϵ_body * σ * ((T_radiant + T_bush) / 2)^3
            Q_rad3 = area_convection * F_vegetation * 4 * ϵ_body * σ * ((T_radiant + T_vegetation) / 2)^3
            Q_rad4 = area_convection * F_ground * 4 * ϵ_body * σ * ((T_radiant + TLOWER) / 2)^3

            if conduction_fraction < 1
                Q_cond = CD * (T_ins_compressed - T_conduction)
                Q_conv = hc * area_convection * (T_insulation - T_air)
                Q_rad = Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4
                Q_env = Q_rad + Q_conv + Q_cond + Q_evap_insulation - Q_sol
            else
                # Fully conductive solution
                Q_env = 0.0
                Q_cond = CD * (T_ins_compressed - T_conduction)
                Q_conv = 0.0
                Q_rad = 0.0
                Q_evap_insulation = 0.0
            end

            # Skin temperature calculations
            T_skin1 = T_core - Q_env * r_flesh^2 / (4 * k_flesh * volume) - Q_env * r_flesh^2 / (2 * k_fat * volume) * log(r_skin / r_flesh)
            T_skin2 = T_insulation # Approximation for now
            T_skin_ave = (T_skin1 + T_skin2) / 2

            ΔT_insulation = abs(T_insulation - T_ins_compressed)
            ΔT_skin = abs(T_skin - T_skin_ave)

            if ΔT_insulation < tolerance_simusol && ΔT_skin < tolerance_simusol
                break
            end

            # Update guesses
            T_insulation = T_ins_compressed
            T_skin = T_skin_ave
            SOLCT += 1
        end
        break
    end
else
    # Bare body case
    # (Translation omitted for brevity, similar structure as above)
    QGENNET = 0.0
end

# Prepare results
RESULTS = [T_insulation, T_skin_ave, Q_conv, Q_cond, QGENNET, Q_evap_skin, Q_rad, Q_sol,
           Q_rad1, Q_rad2, Q_rad3, Q_rad4, Q_evap_insulation, NTRY, success, k_insulation]

return RESULTS

#end