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
mass = 0.04u"kg"
shapeb = 2
shape = Cylinder(mass, density, shapeb) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = (fat, fur)

body = Body(shape, Naked()) # construct a Body, which is furred
body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))

body_geometry = geometry(shape, fur)

θ = 90u"°"
silhouette = silhouette_area(body, θ)

shape = Sphere(mass, density) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = (fat, fur)

body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))


shape = Plate(mass, density, shapeb, shapeb) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = (fat, fur)

body = Body(shape, Naked()) # construct a Body, which is furred
body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))

shape = Ellipsoid(mass, density, shapeb, shapeb) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = (fat, fur)

body = Body(shape, Naked()) # construct a Body, which is furred
body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))

#### solvendo ####

T_air = u"K"(20.0u"°C") # air temperature at local height
T_air_reference = T_air # air temperature at reference height
T_substrate = T_air # ground temperature
T_sky = T_air # sky temperature
T_vegetation = T_air_reference
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
flesh_conductivity_step = 0.1u"W/m/K" # turns on thermal conductivity increase (W/mK), the value being the increment by which flesh_conductivity is increased per iteration
flesh_conductivity_max = 2.8u"W/m/K" # maximum flesh conductivity (W/mK)
pant = 1 # multiplier on breathing rate to simulate panting (-)
pant_step = 0.1 # increment for multiplier on breathing rate to simulate panting (-)
pant_multiplier = 1.05 # multiplier on basal metabolic rate at maximum panting level (-)

# MORPHOLOGY

# geometry
mass = 65.0u"kg" # kg
ρ_body = 1000.0u"kg/m^3" # kg/m3
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
insulation_conductivity = nothing # user-specified fur thermal conductivity (W/mK), not used if 0
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
ventral_fraction = 0.3
conduction_fraction = 0.0
longwave_depth_fraction = 1 # fractional depth of fur at which longwave radiation is exchanged (0-1)

# radiation exchange
ϵ_body = 0.99 # animal emissivity (-)
F_bush_reference = 0.0 # this is for veg below/around animal (at TALOC)
F_ground_reference = 0.5 # reference configuration factor to ground
F_sky_reference = 0.5 # configuration factor to sky
F_vegetation_reference = 0 # configuration factor to sky

# PHYSIOLOGY

# thermal
T_core = u"K"(37.0u"°C") # core temperature (°C)
T_core_maximum = u"K"(39.0u"°C") # maximum core temperature (°C)
flesh_conductivity = 0.9u"W/m/K" # initial thermal conductivity of flesh (0.412 - 2.8 W/m°C)
fat_conductivity = 0.230u"W/m/K" # conductivity of fat (W/mK)
insulation_conductivity = nothing # conductivity of insulation (W/mK)

# evaporation
skin_wetness = 0.005 # part of the skin surface that is wet (fractional)
insulation_wetness = 0.0 # part of the fur/feathers that is wet after rain (fractional)
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
tolerance_simusol = 0.001u"K" # tolerance for SIMULSOL
tolerance_zbrent = 1e-5u"W" # tolerance for ZBRENT
thermoregulate = true # invoke thermoregulatory response
respire = true # compute respiration and associated heat loss
thermoregulation_mode = 1 # 1 = raise core then pant then sweat, 2 = raise core and pant simultaneously, then sweat
torpor = false # go into torpor if possible (drop T_core down to TC_MIN)

bodyshape = Ellipsoid(mass, ρ_body, shape_b, shape_b) # define shape as a Cylinder struct of type 'Shape' and give it required values


# # insulation properties
# insulation_pars = InsulationPars(;
#     fibre_diameter_dorsal, # hair diameter, dorsal (m)
#     fibre_diameter_ventral, # hair diameter, ventral (m)
#     fibre_length_dorsal, # hair length, dorsal (m)
#     fibre_length_ventral, # hair length, ventral (m)
#     insulation_depth_dorsal, # fur depth, dorsal (m)
#     insulation_depth_ventral, # fur depth, ventral (m)
#     fibre_density_dorsal, # hair density, dorsal (1/m2)
#     fibre_density_ventral, # hair density, ventral (1/m2)
#     insulation_reflectance_dorsal,  # fur reflectivity dorsal (fractional, 0-1)
#     insulation_reflectance_ventral,  # fur reflectivity ventral (fractional, 0-1)
#     insulation_depth_compressed, # depth of compressed fur (for conduction) (m)
#     fibre_conductivity, # hair thermal conductivity (W/m°C)
#     longwave_depth_fraction,
# )
# insulation_temperature = T_insulation * 0.7 + T_skin * 0.3
# insulation_out = insulation_properties(; insulation = insulation_pars, insulation_temperature, ventral_fraction)
# effective_conductivities = insulation_out.effective_conductivities # effective thermal conductivity of fur array, mean, dorsal, ventral (W/mK)
# absorption_coefficients = insulation_out.absorption_coefficients # term involved in computing optical thickess (1/mK2)
# optical_thickness_factors = insulation_out.optical_thickness_factors # optical thickness array, mean, dorsal, ventral (m)
# fibre_diameters = insulation_out.fibre_diameters # fur diameter array, mean, dorsal, ventral (m)
# fibre_lengths = insulation_out.fibre_lengths # fur length array, mean, dorsal, ventral (m)
# fibre_densities = insulation_out.fibre_densities # fur density array, mean, dorsal, ventral (1/m2)
# insulation_depths = insulation_out.insulation_depths # fur depth array, mean, dorsal, ventral (m)
# insulation_reflectance = insulation_out.insulation_reflectance # fur reflectivity array, mean, dorsal, ventral (fractional, 0-1)
# insulation_test = insulation_out.insulation_test # test of presence of fur (length x diamater x density x depth) (-)
# insulation_conductivity_compressed = insulation_out.insulation_conductivity_compressed # effective thermal conductivity of compressed ventral fur (W/mK)

# # geometry
# bodyshape = Ellipsoid(mass, density, shape_b, shape_b) # define shape as a Cylinder struct of type 'Shape' and give it required values
# fat = Fat(fat_fraction, fat_density)
# fur = Fur(insulation_depths[1], fibre_diameters[1], fibre_densities[1])
# composite_insulation = (fat, fur)

# geometry_out = Body(bodyshape, CompositeInsulation(fur, fat))

# area_mammal_fur = mammal_fur_area(geometry_out)
# area_mammal_skin = mammal_skin_area(bodyshape)
# area_bird_plumage = bird_plumage_area(bodyshape)
# area_bird_skin = bird_skin_area(bodyshape)

# volume = geometry_out.geometry.volume # volume, m3
# characteristic_dimension = geometry_out.geometry.characteristic_dimension # characteristic dimension for convection, m
# area_total = get_total_area(geometry_out) # total area, m2
# area_silhouette = silhouette_area(bodyshape, geometry_out, 0u"°")
# area_skin = get_skin_area(geometry_out) # area of skin, m2
# area_skin_evaporation = get_evaporation_area(geometry_out) # area of skin for convection/evaporation (total skin area - hair area), m2
# area_evaporation = area_skin_evaporation # area of skin for convection/evaporation (total skin area - hair area), m2
# area_conduction = area_total * conduction_fraction # area of skin for convection/evaporation (total skin area - hair area), m2
# r1 = get_r_skin(geometry_out) # shape-specific core-skin radius in shortest dimension, m
# r2 = get_r_insulation(geometry_out)# shape-specific core-fur/feather interface radius in shortest dimension, m


# # solar radiation normal to sun's rays
# normal_radiation = zenith_angle < 90u"°" ? solar_radiation / cos(zenith_angle) : solar_radiation

# α_body_dorsal = 1 - insulation_reflectance_dorsal #solar absorptivity of dorsal fur (fractional, 0-1)
# α_body_ventral = 1 - insulation_reflectance_ventral # solar absorptivity of ventral fur (fractional, 0-1)

# # correct F_sky for vegetaion overhead
# F_vegetation = F_sky_reference * shade
# F_sky = F_sky_reference - F_vegetation
# F_ground = F_ground_reference

# Q_solar, Q_direct, Q_solar_sky, Q_solar_substrate = solar_out = solar(α_body_dorsal, α_body_ventral, area_silhouette, area_total, area_conduction, F_ground, F_sky, α_substrate, shade, normal_radiation, direct_radiation, diffuse_radiation)
# Q_dorsal = Q_direct + Q_solar_sky
# Q_ventral = Q_solar_substrate
# Q_diffuse = Q_solar_sky + Q_solar_substrate


# (; Q_conv, hc, hd, Sh,
#     Q_free, Q_forc,
#     hc_free, hc_forc,
#     Sh_free, Sh_forc,
#     hd_free, hd_forc) = convection(geometry_out, area_evaporation, T_air, T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)


#####################################
#             solvendo              #
#####################################

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

# infrared radiation properties of insulation
insulation_pars = InsulationPars(;
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
    fibre_conductivity, # hair thermal conductivity (W/m°C)
    longwave_depth_fraction,
)
insulation_temperature = T_insulation * 0.7 + T_skin * 0.3
insulation_out = insulation_properties(; insulation = insulation_pars, insulation_temperature, ventral_fraction)
fibre_diameters = insulation_out.fibre_diameters # fur diameter array, mean, dorsal, ventral (m)
fibre_densities = insulation_out.fibre_densities # fur density array, mean, dorsal, ventral (1/m2)
insulation_depths = insulation_out.insulation_depths # fur depth array, mean, dorsal, ventral (m)

fat = Fat(fat_fraction, ρ_fat)
fur = Fur(insulation_depths[1], fibre_diameters[1], fibre_densities[1])
composite_insulation = (fat, fur)

# geometric properties
geometry_out = Body(bodyshape, CompositeInsulation(fur, fat))

# solar radiation
normal_radiation = zenith_angle < 90u"°" ? solar_radiation / cos(zenith_angle) : solar_radiation

α_body_dorsal = 1 - insulation_reflectance_dorsal #solar absorptivity of dorsal fur (fractional, 0-1)
α_body_ventral = 1 - insulation_reflectance_ventral # solar absorptivity of ventral fur (fractional, 0-1)

# correct F_sky for vegetaion overhead
F_vegetation = F_sky_reference * shade
F_sky = F_sky_reference - F_vegetation
F_ground = F_ground_reference

area_silhouette = silhouette_area(geometry_out, 0u"°")
area_total = get_total_area(geometry_out) # total area, m2
area_conduction = area_total * conduction_fraction # area of skin for convection/evaporation (total skin area - hair area), m2
area_evaporation = get_evaporation_area(geometry_out)
area_convection = area_total * (1 - conduction_fraction)

Q_solar, Q_direct, Q_solar_sky, Q_solar_substrate = solar_out = solar(α_body_dorsal, α_body_ventral, area_silhouette, area_total, area_conduction, F_ground, F_sky, α_substrate, shade, normal_radiation, direct_radiation, diffuse_radiation)
Q_dorsal = Q_direct + Q_solar_sky
Q_ventral = Q_solar_substrate
Q_diffuse = Q_solar_sky + Q_solar_substrate

# convection
(; Q_conv, hc, hd, Sh,
    Q_free, Q_forc,
    hc_free, hc_forc,
    Sh_free, Sh_forc,
    hd_free, hd_forc) = convection(geometry_out, area_convection, T_air, T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)

# set infrared environment
T_veg = T_air_reference # assume vegetation casting shade is at reference (e.g. 1.2m or 2m) air temperature (deg C)
T_lower = T_substrate

side = 1
for side in 1:2

# Calculating solar intensity entering fur. This will depend on whether we are calculating the fur
# temperature for the dorsal side or the ventral side. The dorsal side will have solar inputs from
# the direct beam hitting the silhouette area as well as diffuse solar scattered from the sky.
# The ventral side will have diffuse solar scattered off the substrate.
# Resetting config factors and solar depending on whether the dorsal side (side=1) or ventral side (side=2) is being estimated.
if Q_solar > 0.0u"W" then
 if side == 1
  F_sky = F_sky_reference * 2.0 # proportion of upward view that is sky
  F_vegetation = F_vegetation_reference * 2.0 # proportion of upward view that is vegetation (shade)
  F_ground = 0.0
  F_bush = 0.0
  Q_solar = 2.0 * Q_direct + ((Q_solar_sky / F_sky_reference) * F_sky) # direct x 2 because assuming sun in both directions, and unadjusting Q_solar_sky for config factor imposed in SOLAR_ENDO and back to new larger one in both directions
 else
  F_sky = 0.0
  F_vegetation = 0.0
  F_ground = F_ground_reference * 2.0
  F_bush = F_bush_reference * 2.0
  Q_solar = (Q_ventral / (1.0 - F_sky_reference - F_vegetation_reference)) * (1.0 - (2.0 * conduction_fraction)) # unadjust by config factor imposed in SOLAR_ENDO to have it coming in both directions, but also cutting down according to fractional area conducting to ground (in both directions)
 end
else
 Q_solar = 0.0u"W"
 if side == 1
  F_sky = F_sky_reference * 2.0
  F_vegetation = F_vegetation_reference * 2.0
  F_ground = 0.0
  F_bush = 0.0
 else
  F_sky = 0.0
  F_vegetation = 0.0
  F_ground = F_ground_reference * 2.0
  F_bush = F_bush_reference
 end
end

# set fur depth and conductivity
# index for effective_conductivities, is the average (1), front/dorsal (2), back/ventral(3) of the body part
if Q_solar > 0.0u"W" || (insulation_depths[2] != insulation_depths[3])
 if side == 1
    insulation = Fur(insulation_depths[2], fibre_diameters[2], fibre_densities[2])
 else
    insulation = Fur(insulation_depths[3], fibre_diameters[3], fibre_densities[3])
 end
else
 insulation = Fur(insulation_depths[1], fibre_diameters[1], fibre_densities[1])
end

composite_insulation = (fat, fur)
geometry_out = Body(bodyshape, CompositeInsulation(insulation, fat))
# volume, accounting for any subcutaneous fat
volume = geometry_out.geometry.volume - fat_fraction * mass / ρ_fat

A_total = get_total_area(geometry_out)
r_skin = get_r_skin(geometry_out) # body radius (including fat), m
r_flesh = get_r_flesh(geometry_out) # body radius flesh only (no fat), m
r_insulation = r_skin + insulation.thickness # body radius including fur, m
diameter = 2 * r_insulation # diameter, m
r_radiation = r_skin + (longwave_depth_fraction * insulation.thickness) # effective radiation radius, m
if geometry_out.shape isa Cylinder || geometry_out.shape isa Sphere
  r_compressed = r_skin + insulation_depth_compressed
else
  r_compressed = r_insulation # Note that this value is never used if conduction not being modeled, but need to have a value for the calculations
end

#LEN = ALENTH # length, m

# Calculating the "Cd" variable: Qcond = Cd(Tskin-Tsub), where Cd = Conduction area*ksub/subdepth
if side == 2 # doing ventral side, add conduction
 A_conduction = A_total * (conduction_fraction * 2)
 cd = (A_conduction * k_substrate) / 0.025u"m" # assume conduction happens from 2.5 cm depth
else  # doing dorsal side, no conduction. No need to adjust areas used for convection.
 A_conduction = 0.0u"m^2"
 cd = 0.0u"W/K"
end

# package up inputs
#fur_vars = (; side)#(/LEN,ZFUR,FURTHRMK,KEFF,BETARA,FURTST,insulation_depth,LHAR(S+1),DHAR(S+1),RHOAR(S+1),REFLFR(S+1),KHAIR,REAL(S,8)/)
geom_vars = (; side, cd, conduction_fraction, longwave_depth_fraction, convection_enhancement)#(/SHAPE,SUBQFAT,area_evaporation,volume,diameter,area_evaporation,CONVSK,r_insulation,r_flesh,r_skin,longwave_depth_fraction,r_radiation,ASEMAJ,BSEMIN,CSEMIN,cd,conduction_fraction,r_compressed,r_insulation_compressed,k_compressed,CONV_ENHANCE/)
env_vars = (; fluid_type, T_air, T_substrate, T_bush, T_vegetation, T_lower, T_sky, T_conduction, rh, wind_speed, P_atmos, F_sky, F_ground, F_bush, F_vegetation, Q_solar, fO2, fCO2, fN2)#(/FLTYPE,T_air,TS,T_bush,T_vegetation,T_lower,T_sky,T_conduction,RH,VEL,BP,ELEV,F_sky,F_bush,FAVEG,F_ground,Q_solar,GRAV/)
traits = (; T_core, k_flesh = flesh_conductivity, k_fat = fat_conductivity, ϵ_body, skin_wetness, insulation_wetness, bare_skin_fraction, eye_fraction, insulation_conductivity)#(/T_core,k_flesh,k_fat,EMISAN,FATTHK,FLYHR,FURWET,PCTBAREVAP,PCTEYES/)
    
# function SIMULSOL(tolerance_simusol::Float64, geometry_out, FURVARS::Vector{Float64},
# GEOMVARS::Vector{Float64}, ENVVARS::Vector{Float64},
# TRAITS::Vector{Float64}, T_insulation::Float64, skin_wetness::Float64,
# T_skin::Float64)

function simulsol(tolerance_simusol, geometry_out, insulation_pars, insulation_out, geom_vars, env_vars, traits)

    #fibre_length = insulation_out.fibre_lengths[side + 1]
    #insulation_depth = insulation_out.insulation_depths[side + 1]
    #k_eff = insulation_out.effective_conductivities[side + 1]
    #absorption_coefficient = insulation_out.absorption_coefficients[side + 1]
    insulation_test = u"m"(insulation_out.insulation_test)

    # # Initialize
    # solpro = 1
    # solution_count = 0.0
    # ntry = 0
    # success = true
    # qgennet = 0.0u"W"
    success = true

    T_skin = T_core - 3u"K" # skin temperature (°C)
T_insulation = T_air # fur/air interface temperature (°C)
    if insulation_test > 0.0u"m"
        (; T_insulation, T_skin, Q_convection, Q_conduction, Q_gen_net, Q_evap_skin, 
        Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, 
        Q_rad_ground, Q_evap_insulation, success, ntry) = solve_with_insulation!(T_skin, T_insulation, 
        geometry_out, insulation_pars, insulation_out, geom_vars, env_vars, traits, 
        tolerance_simusol)
        return
    else
        (; T_insulation, T_skin, Q_convection, Q_conduction, Q_gen_net, Q_evap_skin, 
        Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, 
        Q_rad_ground, Q_evap_insulation, success, ntry) = solve_without_insulation!(T_skin, T_insulation, 
            geometry_out, insulation_pars, insulation_out, geom_vars, env_vars, traits, 
            tolerance_simusol)
        return
    end





    # if insulation_test > 0.0u"m"
    #     # Insulated body
    #     while true
    #         ntry += 1
    #         for i in 1:20
    #             # Convective heat transfer
    #             (; hc, hd, hd_free) = convection(geometry_out, area_evaporation, T_air, T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)

    #             # Evaporative heat loss
    #             Q_evap_skin = evaporation(; T_surface = T_skin, wetness = skin_wetness, area = area_evaporation, hd, eye_fraction, T_air, rh, P_atmos, fO2, fCO2, fN2).Q_evap
    #             Q_evap_insulation = evaporation(; T_surface = T_insulation, wetness = insulation_wetness, area = area_evaporation, hd, eye_fraction, T_air, rh, P_atmos, fO2, fCO2, fN2).Q_evap

    #             # Radiation properties
    #             (; effective_conductivities, absorption_coefficients) = insulation_out = insulation_properties(;
    #                 insulation = insulation_pars,
    #                 insulation_temperature = T_insulation * 0.7 + T_skin * 0.3,
    #                 ventral_fraction,
    #                 )

    #             k_eff = effective_conductivities[side + 1]

    #             # Effective fur conductivity
    #             if !isnothing(insulation_conductivity)
    #                 k_insulation = insulation_conductivity
    #             else
    #                 T_rad_approx = T_skin * (1 - longwave_depth_fraction) + T_insulation * longwave_depth_fraction
    #                 k_rad = (16 * σ * T_rad_approx^3) / (3 * absorption_coefficient)
    #                 k_insulation = k_eff + k_rad
    #             end
    #             (; T_radiant, T_ins_compressed, cd1, cd2, cd3, dv1, dv2, dv3, dv4) = radiant_temperature(; body = geometry_out, insulation = insulation_out, insulation_pars, Q_evap = Q_evap_skin, T_core, T_skin, T_conduction, T_insulation, k_flesh, k_fat, cd, longwave_depth_fraction, conduction_fraction)

    #             # Radiative heat fluxes
    #             Q_rad1 = area_evaporation * F_sky * 4 * ϵ_body * σ * ((T_radiant + T_sky) / 2)^3
    #             Q_rad2 = area_evaporation * F_bush * 4 * ϵ_body * σ * ((T_radiant + T_bush) / 2)^3
    #             Q_rad3 = area_evaporation * F_vegetation * 4 * ϵ_body * σ * ((T_radiant + T_vegetation) / 2)^3
    #             Q_rad4 = area_evaporation * F_ground * 4 * ϵ_body * σ * ((T_radiant + T_lower) / 2)^3

    #             if conduction_fraction < 1
    #                 # These calculations are for when there is less than 100% conduction.
    #                 # The term Q_evap_insulation is included for heat lost due to evaporation from
    #                 # the insulation surface
    #                 (; T_insulation_calc, T_radiant2) = insulation_radiant_temperature(; body=geometry_out, insulation=insulation_out, insulation_pars, T_core, T_ins_compressed, T_sky, T_lower, T_vegetation, T_bush, T_conduction, area_evaporation, hc, cd, k_insulation, Q_solar, Q_evap_insulation, Q_rad1, Q_rad2, Q_rad3, Q_rad4, cd1, cd2, cd3, dv1, dv2, dv3, dv4, longwave_depth_fraction, conduction_fraction)
    #                 Q_rad_sky = Q_rad1 * (T_radiant2 - T_sky)
    #                 Q_rad_bush = Q_rad2 * (T_radiant2 - T_bush)
    #                 Q_rad_vegetation = Q_rad3 * (T_radiant2 - T_vegetation)
    #                 Q_rad_ground = Q_rad4 * (T_radiant2 - T_lower)
    #                 Q_radiation = Q_rad_sky + Q_rad_bush + Q_rad_vegetation + Q_rad_ground
    #                 Q_convection = hc * area_evaporation * (T_insulation_calc - T_air)
    #                 Q_conduction = cd * (T_ins_compressed - T_conduction)
    #             else

    #                 (; cf1, T_ins_compressed) = compressed_radiant_temperature(; body=geometry_out, insulation=insulation_out, insulation_pars, k_flesh, k_fat, T_core, T_conduction, cd)
    #                 Q_radiation=0.0u"W"
    #                 Q_convection=0.0u"W"
    #                 Q_evap_insulation = 0.0u"W"
    #                 Q_solar = 0.0u"W"
    #                 Q_conduction=cd*(T_ins_compressed - T_conduction)
    #             end

    #             Q_env = Q_radiation+Q_convection+Q_conduction+Q_evap_insulation-Q_solar
    #             T_skin_mean = mean_skin_temperature(; body=geometry_out, insulation=insulation_out, insulation_pars, Q_env, Q_evap_skin, k_flesh, k_fat, T_core, T_insulation_calc, T_ins_compressed, cd1, cd2, cd3, conduction_fraction)

    #             ΔT_insulation = abs(T_insulation - T_insulation_calc)
    #             ΔT_skin = abs(T_skin - T_skin_mean)
    #                 # CHECK TO SEE IF THE T_insulation GUESS AND THE CALCULATED GUESS ARE SIMILAR
    #             if ΔT_insulation < tolerance_simusol
    #                 return T_insulation   # equivalent to "GO TO 16" → exit with no further changes
    #             else
    #                 update_T_insulation!(T_insulation, T_insulation_calc, ΔT_insulation, tolerance_simusol, solpro)
    #                 T_skin = T_skin_mean
    #                 solution_count += solution_count
    #                 solution_count >= 100
    #                 if solpro != 3
    #                     solution_count = 0
    #                     solpro += solpro
    #                 else
    #                 # even the second way of solving for balance doesn't work, increase tolerance
    #                     if tolerance_simusol <= 0.001
    #                         tolerance_simusol=0.01
    #                         solution_count = 0
    #                         solpro = 1
    #                     else
    #                         success = false
    #                         #QGENNET=0
    #                         #GOTO 150
    #                     end
    #                 end
    #             end
    #         end
    #     end
    # else
    #     # Bare body case
    #     ntry = ntry + 1

    #     i = 1
    #     while i <= 20
    #     @label label125

    #     (; hc, hd, hd_free) = convection(geometry_out, area_evaporation, T_air, T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)

    #     # Evaporative heat loss
    #     Q_evap_skin = evaporation(; T_surface = T_skin, wetness = skin_wetness, area = area_evaporation, hd, eye_fraction, T_air, rh, P_atmos, fO2, fCO2, fN2).Q_evap

    #     # Q_rad variables for radiant exchange
    #     Q_rad1 = area_evaporation * (F_sky  * 4.0 * ϵ_body * σ * ((T_skin + T_sky)/2)^3)
    #     Q_rad2 = area_evaporation * (F_bush * 4.0 * ϵ_body * σ * ((T_skin + T_bush)/2)^3)
    #     Q_rad3 = area_evaporation * (F_vegetation  * 4.0 * ϵ_body * σ * ((T_skin + T_vegetation)/2)^3)
    #     Q_rad4 = area_evaporation * (F_ground  * 4.0 * ϵ_body * σ * ((T_skin + T_lower)/2)^3)

    #     T_skin1 = ((4.0 * k_flesh * volume) / (r_skin^2) * T_core) - Q_evap_skin + hc * area_evaporation * T_air + Q_solar
    #     T_skin2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_lower
    #     T_skin3 = ((4.0 * k_flesh * volume) / (r_skin^2)) + hc * area_evaporation + Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4

    #     T_skin_calc = (T_skin1 + T_skin2) / T_skin3

    #     Q_rad_sky = Q_rad1 * (T_skin_calc - T_sky)
    #     Q_rad_bush = Q_rad2 * (T_skin_calc - T_bush)
    #     Q_rad_vegetation = Q_rad3 * (T_skin_calc - T_vegetation)
    #     Q_rad_ground = Q_rad4 * (T_skin_calc - T_lower)

    #     Q_radiation = Q_rad_sky + Q_rad_bush + Q_rad_vegetation + Q_rad_ground
    #     Q_conduction = hc * area_evaporation * (T_skin_calc - T_air)
    #     Q_env = qrad + qconv - Q_solar

    #     tskdiff = abs(T_skin - T_skin_calc)

    #     # check if tsk guess and calculated tsk are similar
    #     if tskdiff < diftol
    #         # T_insulation and tsk guesses are similar to calculated values
    #         qgennet = ((4.0 * k_flesh * volume) / r_skin^2) * (T_core - T_skin_calc)
    #         @goto label150
    #     else
    #         # try another initial T_skin guess, restart loop
    #         T_skin = T_skin_calc
    #         T_skin_mean = T_skin_calc
    #         T_insulation = T_skin_calc

    #         ntry = ntry + 1

    #         if ntry == 101
    #          if diftol <= 0.001u"K"
    #             tolerance_simusol = 0.01u"K"
    #             ntry = 0
    #         else
    #             # can't find solution → quit
    #             success = false
    #             qgennet = 0.0u"W"
    #             @goto label150
    #         end
    #     end

    #     @goto label125
    # end

    # # Prepare results
    # #RESULTS = [T_insulation, T_skin_ave, Q_conv, Q_cond, QGENNET, Q_evap_skin, Q_rad, Q_sol,
    # #           Q_rad1, Q_rad2, Q_rad3, Q_rad4, Q_evap_insulation, ntry, success, k_insulation]

    return
end