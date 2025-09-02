using HeatExchange
using Microclimate
using ModelParameters
using Unitful, UnitfulMoles
using Unitful: °, rad, °C, K, Pa, J, kJ, W, ml, L, g, kg, cm, m, s, hr, d, mol, R
using Roots
using Test
using Plots

# define the geometry
mass = 0.04kg
ρ_body = 1000.0kg/m^3
shapeb = 3
shapec = 2 / 3
shape_body = Ellipsoid(mass, ρ_body, shapeb, shapec) # define body shape as a Cylinder struct of type 'Shape' and give it required values
geometric_traits = Body(shape_body, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct

# construct the Model which holds the parameters of the organism in the Organism concrete struct, of type AbstractOrganism
lizard = Model(Organism(geometric_traits, MorphoPars(), PhysioPars()))
# get the environmental parameters
environmental_params = EnvironmentalPars()

# specify place and time
lat = -30.0°
lon = 140.0°
elev = 10.0m
days = [15, 45]*1.0
hours = collect(0.:1:24.)
heights = [1.0,]u"cm"

# define daily weather and soil moisture
TIMINS = [0, 0, 1, 1] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
TIMAXS = [1, 1, 0, 0] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
TMINN = [10.0, 8.0]u"°C" # minimum air temperatures (°C)
TMAXX = [30.0, 25]u"°C" # maximum air temperatures (°C)
RHMINN = [20.0, 30.0] # min relative humidity (%)
RHMAXX = [80.0, 90.0] # max relative humidity (%)
WNMINN = [0.1, 0.2]u"m/s" # min wind speed (m/s)
WNMAXX = [1.0, 1.4]u"m/s" # max wind speed (m/s)
CCMINN = [20.0, 23.0] # min cloud cover (%)
CCMAXX = [90.0, 100.0] # max cloud cover (%)
SoilMoist = [0.2, 0.2] # fractional
minshade = 0.0
maxshade = 90.0

# run microclimate model in minshade environment
micro_minshade = runmicro(;
    lat,
    elev,
    heights,
    days,
    hours,
    TMINN,
    TMAXX,
    RHMINN,
    RHMAXX,
    WNMINN,
    WNMAXX,
    CCMINN,
    CCMAXX,
    SoilMoist,
    SHADES = fill(minshade, length(TAIRs)),
)

# run microclimate model in maxshade environment
micro_maxshade = runmicro(;
    lat,
    elev,
    heights,
    days,
    hours,
    TMINN,
    TMAXX,
    RHMINN,
    RHMAXX,
    WNMINN,
    WNMAXX,
    CCMINN,
    CCMAXX,
    SoilMoist,
    SHADES = fill(maxshade, length(TAIRs)),
)

env_minshade = EnvironmentalVarsVec(
    T_air = K.(micro_minshade.airtemperature[:, 2]), # second column is first node above surface
    T_sky = micro_minshade.skytemperature,
    T_sub = micro_minshade.soiltemperature[:, 1], # surface temperature
    rh = micro_minshade.relativehumidity[:, 2], # second column is first node above surface
    vel = micro_minshade.windspeed[:, 2], # second column is first node above surface
    Q_sol = micro_minshade.globalsolar .* (1.0-minshade/100.0),
    Q_dir = micro_minshade.directsolar .* (1.0-minshade/100.0),
    Q_dif = micro_minshade.diffusesolar .* (1.0-minshade/100.0),
    zen = micro_minshade.zenithangle
)

env_maxshade = EnvironmentalVarsVec(
    T_air = K.(micro_maxshade.airtemperature[:, 2]),
    T_sky =micro_maxshade.skytemperature,
    T_sub = micro_maxshade.soiltemperature[:, 1],
    rh = micro_maxshade.relativehumidity[:, 2],
    vel = micro_maxshade.windspeed[:, 2],
    Q_sol = micro_maxshade.globalsolar .* (1.0-minshade/100.0),
    Q_dir = micro_maxshade.directsolar .* (1.0-minshade/100.0),
    Q_dif = micro_maxshade.diffusesolar .* (1.0-minshade/100.0),
    zen = micro_maxshade.zenithangle
)

# set shade
environment = env_minshade

# compute body temperature
balances = map(1:n) do i
    env_i = EnvironmentalVars(
        T_air   = environment.T_air[i],
        T_sky   = environment.T_sky[i],
        T_sub   = environment.T_sub[i],
        rh      = environment.rh[i],
        vel     = environment.vel[i],
        P_atmos = environment.P_atmos[i],
        zen     = environment.zen[i],
        k_sub   = environment.k_sub[i],
        Q_sol   = environment.Q_sol[i],
        Q_dir   = environment.Q_dir[i],
        Q_dif   = environment.Q_dif[i],
    )
    variables_i = (organism = OrganismalVars(), environment = env_i)
    get_Tb(lizard, environmental_params, variables_i)
end
balance_out = flip2vectors(balances); # pull out each output as a vector
resp_out = flip2vectors(balance_out.resp_out); # pull out each output as a vector
evap_out = flip2vectors(balance_out.evap_out); # pull out each output as a vector
conv_out = flip2vectors(balance_out.conv_out); # pull out each output as a vector

plot(1:1:n, °C.(balance_out.T_core), ylims=[0.0, 55.0])
plot!(1:1:n, TAIRs)

plot(1:1:n, u"mg/hr".(evap_out.m_cut))
plot!(1:1:n, u"mg/hr".(evap_out.m_resp))
plot!(1:1:n, u"mg/hr".(evap_out.m_eyes))
plot!(1:1:n, u"mg/hr".(evap_out.m_evap))
