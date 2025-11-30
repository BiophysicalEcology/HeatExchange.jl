library(NicheMapR)

endoR_input = list(
    TA = 40.0, # air temperature at local height (°C)
    TAREF = 40.0, # air temperature at reference height (°C)
    TGRD = 60.0, # ground temperature (°C)
    TSKY = 5.0, # sky temperature (°C)
    VEL = 1.0, # wind speed (m/s)
    RH = 5.0, # relative humidity (%)
    QSOLR = 0.0, # solar radiation, horizontal plane (W/m2)
    Z = 20.0, # zenith angle of sun (degrees from overhead)
    ELEV = 0.0, # elevation (m)
    ABSSB = 0.8, # solar absorptivity of substrate (fractional, 0-1)
    
    # other environmental variables
    FLTYPE = 0, # fluid type: 0 = air; 1 = fresh water; 2 = salt water
    TCONDSB = 35.0, # surface temperature for conduction (°C)
    KSUB = 2.79, # substrate thermal conductivity (W/m°C)
    TBUSH = 40.0, # bush temperature (°C)
    BP = 101325.0, # Pa, negative means elevation is used
    O2GAS = 20.95, # oxygen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)
    N2GAS = 79.02, # nitrogen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)
    CO2GAS = 0.0412, # carbon dioxide concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)
    R_PCO2 = 0.0412 / 100, # reference atmospheric dioxide concentration of air (proportion), to allow for anthropogenic change (\%)
    PDIF = 0.15, # proportion of solar radiation that is diffuse (fractional, 0-1)
    GRAV = 9.80665, # acceleration due to gravity, m/s^2
    
    # BEHAVIOUR
    
    SHADE = 0.0, # shade level (%)
    FLYHR = 0, # is flight occurring this hour? (imposes forced evaporative loss)
    UNCURL = 0.1, # allows the animal to uncurl to SHAPE_B_MAX, the value being the increment SHAPE_B is increased per iteration
    TC_INC = 0.1, # turns on core temperature elevation, the value being the increment by which TC is increased per iteration
    PCTWET_INC = 0.1, # turns on sweating, the value being the increment by which PCTWET is increased per iteration
    PCTWET_MAX = 100.0, # maximum surface area that can be wet (%)
    AK1_INC = 0.1, # turns on thermal conductivity increase (W/mK), the value being the increment by which AK1 is increased per iteration
    AK1_MAX = 2.8, # maximum flesh conductivity (W/mK)
    PANT = 1.0, # multiplier on breathing rate to simulate panting (-)
    PANT_INC = 0.1, # increment for multiplier on breathing rate to simulate panting (-)
    PANT_MULT = 1.05, # multiplier on basal metabolic rate at maximum panting level (-)
    
    # MORPHOLOGY
    
    # geometry
    AMASS = 65.0, # kg
    ANDENS = 1000.0, # kg/m3
    FATDEN = 901.0, # kg/m3
    SUBQFAT = 0, # is subcutaneous fat present? (0 is no, 1 is yes)
    FATPCT = 0.0, # % body fat
    SHAPE = 4, # shape, 1 is cylinder, 2 is sphere, 3 is plate, 4 is ellipsoid
    SHAPE_B = 1.1, # current ratio between long and short axis, must be > 1 (-)
    SHAPE_B_MAX = 5.0, # max possible ratio between long and short axis, must be > 1 (-)
    SHAPE_C = 1.1, # current ratio of length:height (plate)
    PVEN = 0.5, # fraction of surface area that is ventral fur (fractional, 0-1)
    PCOND = 0.0, # fraction of surface area that is touching the substrate (fractional, 0-1)
    SAMODE = 0.0, # if 0, uses surface area for SHAPE parameter geometry, if 1, uses bird skin surface area allometry from Walsberg & King. 1978. JEB 76:185–189, if 2 uses mammal surface area from Stahl 1967.J. App. Physiol. 22, 453–460.
    ORIENT = 0.0, # if 1 = normal to sun's rays (heat maximising), if 2 = parallel to sun's rays (heat minimising), 3 = vertical and changing with solar altitude, or 0 = average
    
    # fur properties
    FURTHRMK = 0.0, # user-specified fur thermal conductivity (W/mK), not used if 0
    DHAIRD = 30E-06, # hair diameter, dorsal (m)
    DHAIRV = 30E-06, # hair diameter, ventral (m)
    LHAIRD = 23.9E-03, # hair length, dorsal (m)
    LHAIRV = 23.9E-03, # hair length, ventral (m)
    ZFURD_MAX = 23.9E-03, # max fur depth, dorsal (m)
    ZFURV_MAX = 23.9E-03, # max fur depth, ventral (m)
    ZFURD = 2E-03, # fur depth, dorsal (m)
    ZFURV = 2E-03, # fur depth, ventral (m)
    RHOD = 3000E+04, # hair density, dorsal (1/m2)
    RHOV = 3000E+04, # hair density, ventral (1/m2)
    REFLD = 0.2,  # fur reflectivity dorsal (fractional, 0-1)
    REFLV = 0.2,  # fur reflectivity ventral (fractional, 0-1)
    ZFURCOMP = 2E-03, # depth of compressed fur (for conduction) (m)
    KHAIR = 0.209, # hair thermal conductivity (W/m°C)
    XR = 1.0, # fractional depth of fur at which longwave radiation is exchanged (0-1)
    
    # radiation exchange
    EMISAN = 0.99, # animal emissivity (-)
    FABUSH = 0.0, # this is for veg below/around animal (at TALOC)
    FGDREF = 0.5, # reference configuration factor to ground
    FSKREF = 0.5, # configuration factor to sky
    
    # PHYSIOLOGY
    
    # thermal
    TC = 37.0, # core temperature (°C)
    TC_MAX = 39.0, # maximum core temperature (°C)
    AK1 = 0.9, # initial thermal conductivity of flesh (0.412 - 2.8 W/m°C)
    AK2 = 0.23, # conductivity of fat (W/mK)
    
    # evaporation
    PCTWET = 0.5, # part of the skin surface that is wet (%)
    FURWET = 5.0, # part of the fur/feathers that is wet after rain (%)
    PCTBAREVAP = 5.0, # surface area for evaporation that is skin, e.g. licking paws (%)
    PCTEYES = 0.03, # surface area made up by the eye (%) - make zero if sleeping
    DELTAR = 0.0, # offset between air temperature and breath (°C)
    RELXIT = 100.0, # relative humidity of exhaled air, %
    
    # metabolism/respiration
    QBASAL = (70 * 65 ^ 0.75) * (4.185 / (24 * 3.6)), # basal heat generation (W) from Kleiber (1947)
    RQ = 0.80, # respiratory quotient (fractional, 0-1)
    EXTREF = 20.0, # O2 extraction efficiency (%)
    PANT_MAX = 5.0, # maximum breathing rate multiplier to simulate panting (-)
    AIRVOL_MAX = 1e12, # maximum absolute breathing rate to simulate panting (L/s), can override PANT_MAX
    PZFUR = 1.0, # # incremental fractional reduction in ZFUR from piloerect state (-) (a value greater than zero triggers piloerection response)
    Q10 = 2.0, # Q10 factor for adjusting BMR for TC
    TC_MIN = 19.0, # minimum core temperature during torpor (TORPOR = 1)
    TORPTOL = 0.05, # allowable tolerance of heat balance as a fraction of torpid metabolic rate
    
    # initial conditions
    TS = 37.0 - 3.0, # skin temperature (°C)
    TFA = 40.0, # fur/air interface temperature (°C)
    
    # other model settings
    CONV_ENHANCE = 1.0, # convective enhancement factor for turbulent conditions, typically 1.4
    DIFTOL = 0.001, # tolerance for SIMULSOL
    BRENTOL = 1e-5, # tolerance for ZBRENT
    THERMOREG = 1, # invoke thermoregulatory response
    RESPIRE = 1, # compute respiration and associated heat loss
    TREGMODE = 1, # 1 = raise core then pant then sweat, 2 = raise core and pant simultaneously, then sweat
    TORPOR = 0, # go into torpor if possible (drop TC down to TC_MIN)
    WRITE_INPUT = 0
)

TAs <- c(endoR_input$TA)
endoR_out = lapply(1:length(TAs), function(x){
endoR(
    TA = TAs[x], # air temperature at local height (°C)
    TAREF = endoR_input$TAREF, # air temperature at reference height (°C)
    TGRD = endoR_input$TGRD, # ground temperature (°C)
    TSKY = endoR_input$TSKY, # sky temperature (°C)
    VEL = endoR_input$VEL, # wind speed (m/s)
    RH = endoR_input$RH, # relative humidity (%)
    QSOLR = endoR_input$QSOLR, # solar radiation, horizontal plane (W/m2)
    Z = endoR_input$Z, # zenith angle of sun (degrees from overhead)
    ELEV = endoR_input$ELEV, # elevation (m)
    ABSSB = endoR_input$ABSSB, # solar absorptivity of substrate (fractional, 0-1)
    
    # other environmental variables
    FLTYPE = endoR_input$FLTYPE, # fluid type: 0 = air; 1 = fresh water; 2 = salt water
    TCONDSB = endoR_input$TCONDSB, # surface temperature for conduction (°C)
    KSUB = endoR_input$KSUB, # substrate thermal conductivity (W/m°C)
    TBUSH = endoR_input$TBUSH, # bush temperature (°C)
    BP = endoR_input$BP, # Pa, negative means elevation is used
    O2GAS = endoR_input$O2GAS, # oxygen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)
    N2GAS = endoR_input$N2GAS, # nitrogen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)
    CO2GAS = endoR_input$CO2GAS, # carbon dioxide concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)
    R_PCO2 = endoR_input$R_PCO2, # reference atmospheric dioxide concentration of air (proportion), to allow for anthropogenic change (\%)
    PDIF = endoR_input$PDIF, # proportion of solar radiation that is diffuse (fractional, 0-1)
    GRAV = endoR_input$GRAV, # acceleration due to gravity, m/s^2
    
    # BEHAVIOUR
    
    SHADE = endoR_input$SHADE, # shade level (%)
    FLYHR = endoR_input$FLYHR, # is flight occurring this hour? (imposes forced evaporative loss)
    UNCURL = endoR_input$UNCURL, # allows the animal to uncurl to SHAPE_B_MAX, the value being the increment SHAPE_B is increased per iteration
    TC_INC = endoR_input$TC_INC, # turns on core temperature elevation, the value being the increment by which TC is increased per iteration
    PCTWET_INC = endoR_input$PCTWET_INC, # turns on sweating, the value being the increment by which PCTWET is increased per iteration
    PCTWET_MAX = endoR_input$PCTWET_MAX, # maximum surface area that can be wet (%)
    AK1_INC = endoR_input$AK1_INC, # turns on thermal conductivity increase (W/mK), the value being the increment by which AK1 is increased per iteration
    AK1_MAX = endoR_input$AK1_MAX, # maximum flesh conductivity (W/mK)
    PANT = endoR_input$PANT, # multiplier on breathing rate to simulate panting (-)
    PANT_INC = endoR_input$PANT_INC, # increment for multiplier on breathing rate to simulate panting (-)
    PANT_MULT = endoR_input$PANT_MULT, # multiplier on basal metabolic rate at maximum panting level (-)
    
    # MORPHOLOGY
    
    # geometry
    AMASS = endoR_input$AMASS, # kg
    ANDENS = endoR_input$ANDENS, # kg/m3
    FATDEN = endoR_input$FATDEN, # kg/m3
    SUBQFAT = endoR_input$SUBQFAT, # is subcutaneous fat present? (0 is no, 1 is yes)
    FATPCT = endoR_input$FATPCT, # % body fat
    SHAPE = endoR_input$SHAPE, # shape, 1 is cylinder, 2 is sphere, 3 is plate, 4 is ellipsoid
    SHAPE_B = endoR_input$SHAPE_B, # current ratio between long and short axis, must be > 1 (-)
    SHAPE_B_MAX = endoR_input$SHAPE_B_MAX, # max possible ratio between long and short axis, must be > 1 (-)
    SHAPE_C = endoR_input$SHAPE_C, # current ratio of length:height (plate)
    PVEN = endoR_input$PVEN, # fraction of surface area that is ventral fur (fractional, 0-1)
    PCOND = endoR_input$PCOND, # fraction of surface area that is touching the substrate (fractional, 0-1)
    SAMODE = endoR_input$SAMODE, # if 0, uses surface area for SHAPE parameter geometry, if 1, uses bird skin surface area allometry from Walsberg & King. 1978. JEB 76:185–189, if 2 uses mammal surface area from Stahl 1967.J. App. Physiol. 22, 453–460.
    ORIENT = endoR_input$ORIENT, # if 1 = normal to sun's rays (heat maximising), if 2 = parallel to sun's rays (heat minimising), 3 = vertical and changing with solar altitude, or 0 = average
    
    # fur properties
    FURTHRMK = endoR_input$FURTHRMK, # user-specified fur thermal conductivity (W/mK), not used if 0
    DHAIRD = endoR_input$DHAIRD, # hair diameter, dorsal (m)
    DHAIRV = endoR_input$DHAIRV, # hair diameter, ventral (m)
    LHAIRD = endoR_input$LHAIRD, # hair length, dorsal (m)
    LHAIRV = endoR_input$LHAIRV, # hair length, ventral (m)
    ZFURD_MAX = endoR_input$ZFURD_MAX, # max fur depth, dorsal (m)
    ZFURV_MAX = endoR_input$ZFURV_MAX, # max fur depth, ventral (m)
    ZFURD = endoR_input$ZFURD, # fur depth, dorsal (m)
    ZFURV = endoR_input$ZFURV, # fur depth, ventral (m)
    RHOD = endoR_input$RHOD, # hair density, dorsal (1/m2)
    RHOV = endoR_input$RHOV, # hair density, ventral (1/m2)
    REFLD = endoR_input$REFLD,  # fur reflectivity dorsal (fractional, 0-1)
    REFLV = endoR_input$REFLV,  # fur reflectivity ventral (fractional, 0-1)
    ZFURCOMP = endoR_input$ZFURCOMP, # depth of compressed fur (for conduction) (m)
    KHAIR = endoR_input$KHAIR, # hair thermal conductivity (W/m°C)
    XR = endoR_input$XR, # fractional depth of fur at which longwave radiation is exchanged (0-1)
    
    # radiation exchange
    EMISAN = endoR_input$EMISAN, # animal emissivity (-)
    FABUSH = endoR_input$FABUSH, # this is for veg below/around animal (at TALOC)
    FGDREF = endoR_input$FGDREF, # reference configuration factor to ground
    FSKREF = endoR_input$FSKREF, # configuration factor to sky
    
    # PHYSIOLOGY
    
    # thermal
    TC = endoR_input$TC, # core temperature (°C)
    TC_MAX = endoR_input$TC_MAX, # maximum core temperature (°C)
    AK1 = endoR_input$AK1, # initial thermal conductivity of flesh (0.412 - 2.8 W/m°C)
    AK2 = endoR_input$AK2, # conductivity of fat (W/mK)
    
    # evaporation
    PCTWET = endoR_input$PCTWET, # part of the skin surface that is wet (%)
    FURWET = endoR_input$FURWET, # part of the fur/feathers that is wet after rain (%)
    PCTBAREVAP = endoR_input$PCTBAREVAP, # surface area for evaporation that is skin, e.g. licking paws (%)
    PCTEYES = endoR_input$PCTEYES, # surface area made up by the eye (%) - make zero if sleeping
    DELTAR = endoR_input$DELTAR, # offset between air temperature and breath (°C)
    RELXIT = endoR_input$RELXIT, # relative humidity of exhaled air, %
    
    # metabolism/respiration
    QBASAL = endoR_input$QBASAL, # basal heat generation (W) from Kleiber (1947)
    RQ = endoR_input$RQ, # respiratory quotient (fractional, 0-1)
    EXTREF = endoR_input$EXTREF, # O2 extraction efficiency (%)
    PANT_MAX = endoR_input$PANT_MAX, # maximum breathing rate multiplier to simulate panting (-)
    AIRVOL_MAX = endoR_input$AIRVOL_MAX, # maximum absolute breathing rate to simulate panting (L/s), can override PANT_MAX
    PZFUR = endoR_input$PZFUR, # # incremental fractional reduction in ZFUR from piloerect state (-) (a value greater than zero triggers piloerection response)
    Q10 = endoR_input$Q10, # Q10 factor for adjusting BMR for TC
    TC_MIN = endoR_input$TC_MIN, # minimum core temperature during torpor (TORPOR = 1)
    TORPTOL = endoR_input$TORPTOL, # allowable tolerance of heat balance as a fraction of torpid metabolic rate
    
    # initial conditions
    TS = endoR_input$TS, # skin temperature (°C)
    TFA = endoR_input$TFA, # fur/air interface temperature (°C)
    
    # other model settings
    CONV_ENHANCE = endoR_input$CONV_ENHANCE, # convective enhancement factor for turbulent conditions, typically 1.4
    DIFTOL = endoR_input$DIFTOL, # tolerance for SIMULSOL
    BRENTOL = endoR_input$BRENTOL, # tolerance for ZBRENT
    THERMOREG = endoR_input$THERMOREG, # invoke thermoregulatory response
    RESPIRE = endoR_input$RESPIRE, # compute respiration and associated heat loss
    TREGMODE = endoR_input$TREGMODE, # 1 = raise core then pant then sweat, 2 = raise core and pant simultaneously, then sweat
    TORPOR = endoR_input$TORPOR, # go into torpor if possible (drop TC down to TC_MIN)
    WRITE_INPUT = endoR_input$WRITE_INPUT
)
})

endoR_out1 <- do.call("rbind", lapply(endoR_out, data.frame)) # turn results into data frame
treg <- endoR_out1[, grep(pattern = "treg", colnames(endoR_out1))]
colnames(treg) <- gsub(colnames(treg), pattern = "treg.", replacement = "")
morph <- endoR_out1[, grep(pattern = "morph", colnames(endoR_out1))]
colnames(morph) <- gsub(colnames(morph), pattern = "morph.", replacement = "")
enbal <- endoR_out1[, grep(pattern = "enbal", colnames(endoR_out1))]
colnames(enbal) <- gsub(colnames(enbal), pattern = "enbal.", replacement = "")
masbal <- endoR_out1[, grep(pattern = "masbal", colnames(endoR_out1))]
colnames(masbal) <- gsub(colnames(masbal), pattern = "masbal.", replacement = "")

write.csv(unlist(endoR_input), file = '../data/endoR_input.csv')
write.csv(names(endoR_input), file = '../data/endoR_input_names.csv')
write.csv(treg, file = '../data/endoR_treg.csv')
write.csv(morph, file = '../data/endoR_morph.csv')
write.csv(enbal, file = '../data/endoR_enbal.csv')
write.csv(masbal, file = '../data/endoR_masbal.csv')
