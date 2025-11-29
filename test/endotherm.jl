using HeatExchange
using ModelParameters
using Unitful, UnitfulMoles
using FluidProperties
#using Plots

# ellipsoid model

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

# solvendo

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
ϵ_body = 0.99 # animal emissivity (-)
F_bush_reference = 0.0 # this is for veg below/around animal (at TALOC)
F_ground_reference = 0.5 # reference configuration factor to ground
F_sky_reference = 0.5 # configuration factor to sky
F_vegetation_reference = 0 # configuration factor to sky

# PHYSIOLOGY

# thermal
T_core = u"K"(37.0u"°C") # core temperature (°C)
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
respiratory_quotient = 0.80 # respiratory quotient (fractional, 0-1)
fO2_extract = 0.2 # O2 extraction efficiency (fractional)
pant_max = 5 # maximum breathing rate multiplier to simulate panting (-)
fur_step = 1 # # incremental fractional reduction in insulation_depth from piloerect state (-) (a value greater than zero triggers piloerection response)
q10 = 2 # q10 factor for adjusting BMR for T_core
T_core_min = 19 # minimum core temperature during torpor (TORPOR = 1)
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

# initial conditions
T_skin = T_core - 3u"K" # skin temperature (°C)
T_insulation = T_air # fur/air interface temperature (°C)

# other model settings
convection_enhancement = 1 # convective enhancement factor for turbulent conditions, typically 1.4
tolerance_simusol = 0.001u"K" # tolerance for SIMULSOL
tolerance_torpor = 0.05 # allowable tolerance of heat balance as a fraction of torpid metabolic rate
thermoregulate = true # invoke thermoregulatory response
respire = true # compute respiration and associated heat loss
thermoregulation_mode = 1 # 1 = raise core then pant then sweat, 2 = raise core and pant simultaneously, then sweat
torpor = false # go into torpor if possible (drop T_core down to TC_MIN)

bodyshape = Ellipsoid(mass, ρ_body, shape_b, shape_b) # define shape as a Cylinder struct of type 'Shape' and give it required values


#function endotherm()
    Q_minimum_ref = Q_minimum
    T_core_ref = T_core
    insulation_depth_dorsal_ref = insulation_depth_dorsal
    insulation_depth_ventral_ref = insulation_depth_ventral
    Q_gen = 0.0u"W"
    pant_cost = 0.0u"W"

    while Q_gen < Q_minimum * 0.995
        insulation_temperature = T_insulation * 0.7 + T_skin * 0.3
        insulation_out = insulation_properties(; insulation=insulation_pars, insulation_temperature, ventral_fraction)
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

        # # convection
        # (; Q_conv, hc, hd, Sh,
        #     Q_free, Q_forc,
        #     hc_free, hc_forc,
        #     Sh_free, Sh_forc,
        #     hd_free, hd_forc) = convection(geometry_out, area_convection, T_air, T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)

        # set infrared environment
        T_veg = T_air_reference # assume vegetation casting shade is at reference (e.g. 1.2m or 2m) air temperature (deg C)
        T_lower = T_substrate

        simulsol_out = Vector{NamedTuple}(undef, 2)
        for side in 1:2

            # Calculating solar intensity entering fur. This will depend on whether we are calculating the fur
            # temperature for the dorsal side or the ventral side. The dorsal side will have solar inputs from
            # the direct beam hitting the silhouette area as well as diffuse solar scattered from the sky.
            # The ventral side will have diffuse solar scattered off the substrate.
            # Resetting config factors and solar depending on whether the dorsal side (side=1) or ventral side (side=2) is being estimated.
            if Q_solar > 0.0u"W"
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

            if side == 1
                insulation_conductivity = insulation_conductivity_dorsal
            else
                insulation_conductivity = insulation_conductivity_ventral                
            end
            composite_insulation = (fat, fur)
            geometry_out = Body(bodyshape, CompositeInsulation(insulation, fat))
            # volume, accounting for any subcutaneous fat
            #volume = geometry_out.geometry.volume - fat_fraction * mass / ρ_fat

            A_total = get_total_area(geometry_out)
            r_skin = get_r_skin(geometry_out) # body radius (including fat), m
            #r_flesh = get_r_flesh(geometry_out) # body radius flesh only (no fat), m
            r_insulation = r_skin + insulation.thickness # body radius including fur, m
            #diameter = 2 * r_insulation # diameter, m
            #r_radiation = r_skin + (longwave_depth_fraction * insulation.thickness) # effective radiation radius, m
            if geometry_out.shape isa Cylinder || geometry_out.shape isa Sphere
                r_compressed = r_skin + insulation_depth_compressed
            else
                r_compressed = r_insulation # Note that this value is never used if conduction not being modeled, but need to have a value for the calculations
            end

            #LEN = ALENTH # length, m

            # Calculating the "cd" variable: Qcond = cd(Tskin-Tsub), where cd = Conduction area*ksub/subdepth
            if side == 2 # doing ventral side, add conduction
                A_conduction = A_total * (conduction_fraction * 2)
                cd = (A_conduction * k_substrate) / 0.025u"m" # assume conduction happens from 2.5 cm depth
            else  # doing dorsal side, no conduction. No need to adjust areas used for convection.
                A_conduction = 0.0u"m^2"
                cd = 0.0u"W/K"
            end

            # package up inputs
            geom_vars = (; side, cd, conduction_fraction, longwave_depth_fraction, convection_enhancement)
            env_vars = (; fluid_type, T_air, T_substrate, T_bush, T_vegetation, T_lower, T_sky, T_conduction,
                rh, wind_speed, P_atmos, F_sky, F_ground, F_bush, F_vegetation, Q_solar, fO2, fCO2, fN2)
            traits = (; T_core, k_flesh, k_fat, ϵ_body, skin_wetness,
                insulation_wetness, bare_skin_fraction, eye_fraction, insulation_conductivity)

            simulsol_out[side,] = simulsol(; T_skin, T_insulation, tolerance_simusol, geometry_out, insulation_pars, insulation_out, geom_vars,
                env_vars, traits)
        end

        T_skin_max = max(simulsol_out[1].T_skin, simulsol_out[2].T_skin)

        # zbrent and respfun

        # Now compute a weighted mean heat generation for all the parts/components = (dorsal value *(F_sky+F_vegetation))+(ventral value*F_ground)
        gen_d = simulsol_out[1].Q_gen_net
        gen_v = simulsol_out[2].Q_gen_net
        dmult = F_sky_reference + F_vegetation_reference
        vmult = 1 - dmult # assume that reflectivity of veg below equals reflectivity of soil so vmult left as 1 - dmult
        x = gen_d * dmult + gen_v * vmult # weighted estimate of metabolic heat generation
        Q_sum = x

        # reset configuration factors
        F_bush = F_bush_reference # nearby bush
        F_sky = F_sky_reference # sky
        F_ground = F_ground_reference # ground
        F_vegetation = F_vegetation_reference # vegetation

        # lung temperature and temperature of exhaled air
        T_skin = (simulsol_out[1].T_skin + simulsol_out[2].T_skin) * 0.5
        T_insulation = (simulsol_out[1].T_insulation + simulsol_out[2].T_insulation) * 0.5
        T_lung = (T_core + T_skin) * 0.5 # average of skin and core
        T_air_exit = min(T_air + Δ_breath, T_lung) # temperature of exhaled air, deg C

        if respire
            # now guess for metabolic rate that balances the heat budget while allowing metabolic rate
            # to remain at or above Q_basal, via root-finder ZBRENT
            Q_min = Q_minimum
            Q_m1 = Q_minimum * (-2.)
            Q_m2 = Q_minimum * 10.
            if T_skin_max >= T_core
                Q_m2 = Q_minimum * 1.01
            end
            #tolerance = Q_minimum * tolerance_zbrent

            Q_gen = find_zero(x -> respiration_endotherm(x; T_air_reference=T_air, fO2, fN2, fCO2, P_atmos, Q_min, respiratory_quotient,
                    T_lung, mass, fO2_extract, rh, rh_exit, T_air_exit, pant, Q_sum).balance, (Q_m1, Q_m2), Bisection())

            zbrent_out = respiration_endotherm(Q_gen; T_air_reference=T_air, fO2, fN2, fCO2, P_atmos, Q_min, respiratory_quotient,
                T_lung, mass, fO2_extract, rh, rh_exit, T_air_exit, pant, Q_sum)

            Q_gen = zbrent_out.Q_gen # Q_gen_net
        else
            Q_gen = Q_sum
        end
        # store last-step values
        shape_b_last = shape_b
        k_flesh_last = k_flesh
        T_core_last = T_core
        pant_last = pant
        skin_wetness_last = skin_wetness

        if thermoregulate

            if (insulation_depth_dorsal > insulation_depth_dorsal_ref) &&
               (insulation_depth_ventral > insulation_depth_ventral_ref)
                insulation_depth_dorsal = max(insulation_depth_dorsal_ref, insulation_depth_dorsal -
                                                                           insulation_step * fibre_length_dorsal)
                insulation_depth_ventral = max(insulation_depth_ventral_ref, insulation_depth_ventral -
                                                                             insulation_step * fibre_length_ventral)
            else
                insulation_depth_dorsal = insulation_depth_dorsal_ref
                insulation_depth_ventral = insulation_depth_ventral_ref
                if shape_b < shape_b_max
                    shape_b += shape_b_step
                else
                    shape_b = shape_b_max

                    if k_flesh < k_flesh_max
                        k_flesh += k_flesh_step
                    else
                        k_flesh = k_flesh_max

                        if T_core < T_core_max
                            T_core += T_core_step
                            q10mult = q10^((ustrip(u"K", (T_core - T_core_ref))) / 10)

                            if (thermoregulation_mode >= 2) && (pant < pant_max)
                                pant += pant_step
                                pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                            (pant_multiplier - 1) * Q_minimum_ref

                                if thermoregulation_mode == 3
                                    skin_wetness += skin_wetness_step
                                    if skin_wetness > skin_wetness_max
                                        skin_wetness = skin_wetness_max
                                    end
                                end
                            end

                            Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                        else
                            T_core = T_core_max
                            q10mult = q10^((ustrip(u"K", (T_core - T_core_ref))) / 10)

                            if pant < pant_max
                                pant += pant_step
                                pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                            (pant_multiplier - 1) * Q_minimum_ref
                                Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                                if thermoregulation_mode == 3
                                    skin_wetness += skin_wetness_step
                                    if skin_wetness > skin_wetness_max
                                        skin_wetness = skin_wetness_max
                                    end
                                end
                            else
                                pant = pant_max
                                pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                            (pant_multiplier - 1) * Q_minimum_ref

                                Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                                skin_wetness += skin_wetness_step
                                if (skin_wetness > skin_wetness_max) || (skin_wetness_step <= 0)
                                    skin_wetness = skin_wetness_max
                                    return
                                end
                            end
                        end
                    end
                end
            end
        end
    end
#end