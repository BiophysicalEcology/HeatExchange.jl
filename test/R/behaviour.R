library(NicheMapR)

loc <- c(140, -35)
x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))
ALREF <- abs(trunc(x[1]))
HEMIS <- ifelse(x[2]<0,2.,1.) # 1 is northern hemisphere
# break decimal degree lat/lon into deg and min
ALAT <- abs(trunc(x[2]))
AMINUT <- (abs(x[2])-ALAT)*60
ALONG <- abs(trunc(x[1]))
ALMINT <- (abs(x[1])-ALONG)*60

writecsv <- 0 # make Fortran program write output as csv files
microdaily <- 0 # run microclimate model as normal, where each day is iterated 3 times starting with the initial condition of uniform soil temp at mean monthly temperature
runshade <- 1 # run the model twice, once for each shade level (1) or just for the first shade level (0)?
runmoist <- 0 # run soil moisture model (0=no, 1=yes)?
snowmodel <- 0 # run the snow model (0=no, 1=yes)? - note that this runs slower
hourly <- 0 # run from hourly weather input (0=no, 1=yes)?
rainhourly <- 0 # run from hourly weather input (0=no, 1=yes)?
IR <- 0 # compute clear-sky longwave radiation using Campbell and Norman (1998) eq. 10.10 (includes humidity)
message <- 0 # do not allow the Fortran integrator to output warnings
fail <- 24 * 365 # how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)

doynum <- 2 # number of time intervals to generate predictions for over a year (must be 12 <= x <=365)
doy<-c(15, 46) # middle day of each month
idayst <- 1 # start month
ida <- 2 # end month
EC <- 0.0167238 # Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)

RUF <- 0.004 # Roughness height (m), , e.g. sand is 0.05, grass may be 2.0, current allowed range: 0.001 (snow) - 2.0 cm.
Refhyt <- 2 # Reference height (m), reference height at which air temperature, wind speed and relative humidity input data are measured
Usrhyt <- 0.01# local height (m) at which air temperature, relative humidity and wind speed calculations will be made
ZH <- 0 # heat transfer roughness height (m) for Campbell and Norman air temperature/wind speed profile (invoked if greater than 1, 0.02 * canopy height in m if unknown)
D0 <- 0 # zero plane displacement correction factor (m) for Campbell and Norman air temperature/wind speed profile (0.6 * canopy height in m if unknown)
# Next four parameters are segmented velocity profiles due to bushes, rocks etc. on the surface
#IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO! (then roughness height is based on the parameter RUF)
Z01 <- 0 # Top (1st) segment roughness height(m)
Z02 <- 0 # 2nd segment roughness height(m)
ZH1 <- 0 # Top of (1st) segment, height above surface(m)
ZH2 <- 0 # 2nd segment, height above surface(m)

SLE <- 0.96 # substrate longwave IR emissivity (decimal %), typically close to 1
REFL <- 0.20 # substrate solar reflectivity (decimal %)
CMH2O <- 1 # precipitable cm H2O in air column, 0.1 = VERY DRY; 1.0 = MOIST AIR CONDITIONS; 2.0 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)
# Aerosol extinction coefficient profile
# the values extracted from GADS for Madison
TAI<-c(0.269904738,0.266147825, 0.262442906, 0.258789404, 0.255186744, 0.251634356, 0.248131676, 0.2412732, 0.234606887, 0.228128378, 0.221833385, 0.215717692, 0.20977715, 0.204007681, 0.198405272, 0.187685927, 0.177588357, 0.168082846, 0.159140695, 0.150734206, 0.142836655, 0.135422274, 0.128466227, 0.12194459, 0.115834329, 0.110113284, 0.104760141, 0.099754417, 0.09507644, 0.090707328, 0.086628967, 0.082823998, 0.07927579, 0.075968428, 0.072886691, 0.070016034, 0.067342571, 0.064853053, 0.062534858, 0.060375964, 0.058364941, 0.056490925, 0.054743609, 0.053113222, 0.051590514, 0.050166738, 0.046408775, 0.045302803, 0.044259051, 0.043271471, 0.042334415, 0.041442618, 0.040591184, 0.039775572, 0.038991583, 0.038235345, 0.037503301, 0.036792197, 0.036099067, 0.034101935, 0.033456388, 0.032817888, 0.032184949, 0.031556287, 0.030930816, 0.030307633, 0.029065372, 0.027825562, 0.027205981, 0.026586556, 0.025967391, 0.025348692, 0.024114005, 0.023498886, 0.021669152, 0.021066668, 0.019292088, 0.018144698, 0.016762709, 0.015451481, 0.014949794, 0.014224263, 0.013093462, 0.012670686, 0.012070223, 0.011164062, 0.010241734, 0.009731103, 0.009507687, 0.009212683, 0.008965785, 0.008827751, 0.008710756, 0.008574128, 0.008462605, 0.008446967, 0.008539475, 0.009015237, 0.009748444, 0.010586023, 0.011359647, 0.011901268, 0.012062153, 0.011735443, 0.010882215, 0.009561062, 0.007961182, 0.006438984, 0.005558204, 0.006133532, 0.009277754)

ALTT <- 10 # altitude (m)
slope <- 0 # slope (degrees, range 0-90)
azmuth <- 180 # aspect (degrees, 0 = North, range 0-360)
hori <- rep(0, 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
VIEWF <- 1 - sum(sin(hori * pi / 180)) / length(hori) # convert horizon angles to radians and calc view factor(s)
solonly <- 0 # Only run SOLRAD to get solar radiation? 1=yes, 0=no
lamb <- 0 # Return wavelength-specific solar radiation output?
IUV <- 0 # Use gamma function for scattered solar radiation? (computationally intensive)
ndmax <- 3 # iterations of first day to get a steady periodic
minshade <- 0 # minimum available shade (%)
maxshade <- 90 # maximum available shade (%)
PCTWET <- 0 # percentage of surface area acting as a free water surface (%)

DEP <- c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
ERR <- 1.5 # Integrator error for soil temperature calculations

TIMINS <- c(0, 0, 1, 1)   # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
TIMAXS <- c(1, 1, 0, 0)   # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
TMINN <- c(10.0, 8.0) # minimum air temperatures (°C)
TMAXX <- c(30.0, 25) # maximum air temperatures (°C)
RHMINN <- c(20.0, 30.0) # min relative humidity (%)
RHMAXX <- c(80.0, 90.0) # max relative humidity (%)
WNMINN <- c(0.1, 0.2) # min wind speed (m/s)
WNMAXX <- c(1.0, 1.4) # max wind speed (m/s)
CCMINN <- c(20.0, 23.0) # min cloud cover (%)
CCMAXX <- c(90.0, 100.0) # max cloud cover (%)
RAINFALL <- c(28, 28.2) # monthly mean rainfall (mm)
TAIRhr <- rep(0, 24*doynum) # hourly air temperatures (°C), not used unless 'hourly=1'
RHhr <- rep(0, 24*doynum) # hourly relative humidity (%), not used unless 'hourly=1'
WNhr <- rep(0, 24*doynum) # hourly wind speed (m/s), not used unless 'hourly=1'
CLDhr <- rep(0, 24*doynum) # hourly cloud cover (%), not used unless 'hourly=1'
SOLRhr <- rep(0, 24*doynum) # hourly solar radiation (W/m2, not used unless 'hourly=1'
RAINhr <- rep(0, 24*doynum) # hourly rainfall (mm), not used unless 'hourly=1'
ZENhr <- rep(-1, 24*doynum) # hourly zenith angle (degrees), not used unless 'hourly=1'
IRDhr <- rep(-1, 24*doynum) # hourly downwelling longwave radiation (W/m2), not used if '-1'
tannul <- mean(c(TMAXX, TMINN)) # annual mean temperature for getting monthly deep soil temperature (°C)
tannulrun <- rep(tannul, doynum) # monthly deep soil temperature (2m) (°C)
SoilMoist <- c(0.2, 0.2) # soil moisture (decimal %, 1 means saturated)
# creating the arrays of environmental variables that are assumed not to change with month for this simulation
MAXSHADES <- rep(maxshade, doynum) # daily max shade (%)
MINSHADES <- rep(minshade, doynum) # daily min shade (%)
SLES <- rep(SLE, doynum) # set up vector of ground emissivities for each day
REFLS <- rep(REFL, doynum) # set up vector of soil reflectances for each day
PCTWET <- rep(PCTWET, doynum) # set up vector of soil wetness for each day

# set up a profile of soil properites with depth for each day to be run
Numtyps <- 1 # number of soil types
Nodes <- matrix(data = 0, nrow = 10, ncol = doynum) # array of all possible soil nodes
Nodes[1, 1:doynum] <- 10 # deepest node for first substrate type

# soil thermal parameters 
Thcond <- 1.25 # soil minerals thermal conductivity (W/mC)
Density <- 2.560 # soil minerals density (Mg/m3)
SpecHeat <- 870 # soil minerals specific heat (J/kg-K)
BulkDensity <- 2.56 # soil bulk density (Mg/m3)
SatWater <- 0.26 # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)

# now make the depth-specific soil properties matrix
# columns are:
#1) bulk density (Mg/m3)
#2) volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
#3) thermal conductivity (W/mK)
#4) specific heat capacity (J/kg-K)
#5) mineral density (Mg/m3)
soilprops<-matrix(data = 0, nrow = 10, ncol = 5) # create an empty soil properties matrix
soilprops[1, 1]<-BulkDensity # insert soil bulk density to profile 1
soilprops[1, 2]<-SatWater # insert saturated water content to profile 1
soilprops[1, 3]<-Thcond # insert thermal conductivity to profile 1
soilprops[1, 4]<-SpecHeat # insert specific heat to profile 1
soilprops[1, 5]<-Density # insert mineral density to profile 1
soilinit<-rep(tannul, 20) # make initial soil temps equal to mean annual

# note that these are set for sand (Table 9.1 in Campbell and Norman, 1995)
PE <- rep(0.7, 19) #air entry potential J/kg
KS <- rep(0.0058, 19) #saturated conductivity, kg s/m3
BB <- rep(1.7, 19) #soil 'b' parameter
BD <- rep(1.3, 19) # soil bulk density, Mg/m3
DD <- rep(Density, 19)
L <- c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000 # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131)
R1 <- 0.001 #root radius, m}\cr\cr
RW <- 2.5e+10 #resistance per unit length of root, m3 kg-1 s-1
RL <- 2e+6 #resistance per unit length of leaf, m3 kg-1 s-1
PC <- -1500 #critical leaf water potential for stomatal closure, J kg-1
SP <- 10 #stability parameter for stomatal closure equation, -
IM <- 1e-06 #maximum allowable mass balance error, kg
MAXCOUNT <- 500 #maximum iterations for mass balance, -
LAI <- rep(0.1, doynum) # leaf area index, used to partition traspiration/evaporation from PET
rainmult <- 1 # rainfall multiplier to impose catchment
maxpool <- 10 # max depth for water pooling on the surface, mm (to account for runoff)
evenrain <- 1 # spread daily rainfall evenly across 24hrs (1) or one event at midnight (2)
SoilMoist_Init <- rep(0.2, 10) # initial soil water content for each node, m3/m3
moists <- matrix(nrow = 10, ncol = doynum, data = 0) # set up an empty vector for soil moisture values through time
moists[1:10, ] <- SoilMoist_Init # insert initial soil moisture
spinup <- 1 # repeat first day ndmax times for steady periodic
dewrain <- 0 # don't feed dew back into soil as rain
moiststep <- 360 # how many steps within the hour is soil moisture solved over
maxsurf <- 95 # what is the maximum allowable soil surface temp (for stability purposes), deg C

snowtemp <- 1.5 # temperature at which precipitation falls as snow (used for snow model)
snowdens <- 0.375 # snow density (mg/m3)
densfun <- c(0.5979, 0.2178, 0.001, 0.0038) # slope and intercept of linear model of snow density as a function of day of year - if it is c(0,0) then fixed density used
snowmelt <- 1 # proportion of calculated snowmelt that doesn't refreeze
undercatch <- 1 # undercatch multipier for converting rainfall to snow
rainmelt <- 0.0125 # parameter in equation from Anderson's SNOW-17 model that melts snow with rainfall as a function of air temp
snowcond <- 0 # effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)
intercept <- 0 # snow interception fraction for when there's shade (0-1)
grasshade <- 0 # if 1, means shade is removed when snow is present, because shade is cast by grass/low veg

# intertidal simulation input vector (col 1 = tide in(1)/out(0), col 2 = sea water temperature in °C, col 3 = % wet from wave splash)
tides <- matrix(data = 0, nrow = 24 * doynum, ncol = 3) # matrix for tides

# input parameter vector
microinput<-c(doynum, RUF, ERR, Usrhyt, Refhyt, Numtyps, Z01, Z02, ZH1, ZH2, idayst, ida, HEMIS, ALAT, AMINUT, ALONG, ALMINT, ALREF, slope, azmuth, ALTT, CMH2O, microdaily, tannul, EC, VIEWF, snowtemp, snowdens, snowmelt, undercatch, rainmult, runshade, runmoist, maxpool, evenrain, snowmodel, rainmelt, writecsv, densfun, hourly, rainhourly, lamb, IUV, RW, PC, RL, SP, R1, IM, MAXCOUNT, IR, message, fail, snowcond, intercept, grasshade, solonly, ZH, D0, TIMAXS, TIMINS, spinup, dewrain, moiststep, maxsurf, ndmax)

# Final input list - all these variables are expected by the input argument of the Fortran microclimate subroutine
micro_in<-list(microinput = microinput, tides = tides, doy = doy, SLES = SLES, DEP = DEP, Nodes = Nodes, MAXSHADES = MAXSHADES, MINSHADES = MINSHADES, TMAXX = TMAXX, TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN, CCMAXX = CCMAXX, CCMINN = CCMINN, WNMAXX = WNMAXX, WNMINN = WNMINN, TAIRhr = TAIRhr, RHhr = RHhr, WNhr = WNhr, CLDhr = CLDhr, SOLRhr = SOLRhr, RAINhr = RAINhr, ZENhr = ZENhr, IRDhr = IRDhr, REFLS = REFLS, PCTWET = PCTWET, soilinit = soilinit, hori = hori, TAI = TAI, soilprops = soilprops, moists = moists, RAINFALL = RAINFALL, tannulrun = tannulrun, PE = PE, KS = KS, BB = BB, BD = BD, DD = DD, L = L, LAI = LAI)

### Executing the microclimate model
microut<-microclimate(micro_in)

metout <- microut$metout # retrieve above ground microclimatic conditions, min shade
shadmet <- microut$shadmet # retrieve above ground microclimatic conditions, max shade
soil <- microut$soil # retrieve soil temperatures, minimum shade
shadsoil <- microut$shadsoil # retrieve soil temperatures, maximum shade
tcond <- microut$tcond
shadtcond <- microut$shadtcond
specheat <- microut$specheat
shadspecheat <- microut$shadspecheat
densit <- microut$densit
shaddensit <- microut$shaddensit

if(runmoist == 1){
  soilmoist <- microut$soilmoist # retrieve soil moisture, minimum shade
  shadmoist <- microut$shadmoist # retrieve soil moisture, maximum shade
  humid <- microut$humid # retrieve soil humidity, minimum shade
  shadhumid <- microut$shadhumid # retrieve soil humidity, maximum shade
  soilpot <- microut$soilpot # retrieve soil water potential, minimum shade
  shadpot <- microut$shadpot # retrieve soil water potential, maximum shade
  plant <- microut$plant # retrieve plant output, minimum shade
  shadplant <- microut$shadplant # retrieve plant output, maximum shade
}else{
  soilpot <- soil
  soilmoist <- soil
  shadpot <- soil
  shadmoist <- soil
  humid <- soil
  shadhumid <- soil
  plant <- cbind(soil,soil[,3:4])
  shadplant <- cbind(soil,soil[,3:4])
  soilpot[,3:12] <- 0
  soilmoist[,3:12] <- 0.5
  shadpot[,3:12] <- 0
  shadmoist[,3:12] <- 0.5
  humid[,3:12] <- 0.99
  shadhumid[,3:12] <- 0.99
  plant[,3:14] <- 0
  shadplant[,3:14] <- 0
}
if(snowmodel == 1){
  sunsnow <- microut$sunsnow
  shdsnow <- microut$shdsnow
}
ndays <- length(TMAXX)
nyears <- 1
timeinterval <- ndays
days <- rep(seq(1,timeinterval), 24)
days <- days[order(days)]
dates <- days + metout[, 2] / 60 / 24 - 1 # dates for hourly output
dates2 <- seq(1, timeinterval) # dates for daily output
diffuse_frac <- rep(0.15, ndays * 24)
if(lamb == 1){
  drlam <- as.data.frame(microut$drlam) # retrieve direct solar irradiance
  drrlam <- as.data.frame(microut$drrlam) # retrieve direct Rayleigh component solar irradiance
  srlam <- as.data.frame(microut$srlam) # retrieve scattered solar irradiance
  if(snowmodel == 1){
    micro <- list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,tcond=tcond,shadtcond=shadtcond,specheat=specheat,shadspecheat=shadspecheat,densit=densit,shaddensit=shaddensit,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=MINSHADES,maxshade=MAXSHADES,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam,dates=dates,dates2=dates2,PE=PE,BD=BD,DD=DD,BB=BB,KS=KS, diffuse_frac = diffuse_frac)
  }else{
    micro <- list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,tcond=tcond,shadtcond=shadtcond,specheat=specheat,shadspecheat=shadspecheat,densit=densit,shaddensit=shaddensit,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=MINSHADES,maxshade=MAXSHADES,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam,dates=dates,dates2=dates2,PE=PE,BD=BD,DD=DD,BB=BB,KS=KS, diffuse_frac = diffuse_frac)
  }
}else{
  if(snowmodel == 1){
    micro <- list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,tcond=tcond,shadtcond=shadtcond,specheat=specheat,shadspecheat=shadspecheat,densit=densit,shaddensit=shaddensit,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=MINSHADES,maxshade=MAXSHADES,DEP=DEP,dates=dates,dates2=dates2,PE=PE,BD=BD,DD=DD,BB=BB,KS=KS, diffuse_frac = diffuse_frac)
  }else{
    micro <- list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,tcond=tcond,shadtcond=shadtcond,specheat=specheat,shadspecheat=shadspecheat,densit=densit,shaddensit=shaddensit,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=MINSHADES,maxshade=MAXSHADES,DEP=DEP,dates=dates,dates2=dates2,PE=PE,BD=BD,DD=DD,BB=BB,KS=KS, diffuse_frac = diffuse_frac)
  }
}


plot(soil[, 3], type = 'l')

Ww_g = 40
burrow = 0
shade_seek = 0
diurn = 0
nocturn = 0
crepus = 0
delta_shade = 1
mindepth = 2
T_F_min = 24.5
T_F_max = 34.5
T_pref = 34.5
CT_min = 0.1
CT_max = T_F_max
pct_cond = 0
live = 0
shape = 2
shape_b = 3
shape_c = 2/3

ecto <- ectotherm(Ww_g = Ww_g,
                  live = live,
                  shape = shape,
                  shape_b = shape_b,
                  shape_c = shape_c,
                  burrow = burrow,
                  shade_seek = shade_seek,
                  diurn = diurn,
                  nocturn = nocturn,
                  crepus = crepus,
                  delta_shade = delta_shade,
                  mindepth = mindepth,
                  T_F_min = T_F_min,
                  T_F_max = T_F_max,
                  T_pref = T_pref,
                  T_B_min = T_F_min,
                  T_RB_min = T_F_min,
                  CT_min = CT_min,
                  CT_max = CT_max,
                  pct_cond = pct_cond,
                  nyears = 1,
                  minshades = rep(minshade, 24 * length(TMAXX)),
                  maxshades = rep(maxshade, 24 * length(TMAXX)),
                  alpha_sub = rep(1 - REFL)
                  )

environ <- as.data.frame(ecto$environ)

plot(environ$TC, type = 'l', ylim = c(0, 55))
