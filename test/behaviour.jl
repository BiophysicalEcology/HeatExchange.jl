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
# define body shape as an ellipsoid struct of type 'Shape' and give it required values
shape_body = Ellipsoid(mass, ρ_body, shapeb, shapec) 
# construct a Body, which is naked - this constructor will apply the 'geometry' function 
# to the inputs and return a struct that has the struct for the 'Shape' type, as well 
# as the insulation and the geometry struct
geometric_traits = Body(shape_body, Naked()) 
# construct the Model which holds the parameters of the organism in the Organism concrete struct, 
# of type AbstractOrganism
lizard = Model(Organism(geometric_traits, MorphoPars(), PhysioPars()))

# specify place and time
lat = -30.0°
lon = 140.0°
elev = 10.0m
days = [15, 45]*1.0
hours = collect(0.:1:24.)
heights = [100.0,]u"cm"
α_sub = 0.8

# set the environmental parameters
environmental_params = EnvironmentalPars(
    elev = elev,
    α_sub = Param(α_sub, bounds=(0.0, 1.0)),
)

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
shade = 0.0

# run microclimate model in minshade environment
micro = runmicro(;
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
    SHADES = fill(shade, length(TAIRs)),
)

environment = EnvironmentalVarsVec(
    T_air = K.(micro.airtemperature[:, 2]), # second column is first node above surface
    T_sky = micro.skytemperature,
    T_sub = micro.soiltemperature[:, 1], # surface temperature
    rh = micro.relativehumidity[:, 2], # second column is first node above surface
    vel = micro.windspeed[:, 2], # second column is first node above surface
    Q_sol = micro.globalsolar .* (1.0-shade/100.0),
    Q_dir = micro.directsolar .* (1.0-shade/100.0),
    Q_dif = micro.diffusesolar .* (1.0-shade/100.0),
    zen = micro.zenithangle
)

# compute body temperature
n = length(TAIRs)
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

plot(1:1:n, °C.(balance_out.T_core), ylims=[-10.0, 55.0])
plot!(1:1:n, micro.airtemperature[:, 2])
plot!(1:1:n, TAIRs)
plot!(1:1:n, micro.soiltemperature[:, 1])
plot!(1:1:n, °C.(micro.skytemperature[:, 1]))