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

# specify place and time
latitude = -30.0°
longitude = 140.0°
elevation = 10.0m
days = [15, 45]*1.0
hours = collect(0.:1:24.)
heights = [1.0,]u"cm"
α_sub = 0.8

# set the environmental parameters
environmental_params = EnvironmentalPars(
    elev = elev,
    α_sub = Param(α_sub, bounds=(0.0, 1.0)),
)

# define daily weather and soil moisture
minima_times = [0, 0, 1, 1] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
maxima_times = [1, 1, 0, 0] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
air_temperature_min = [10.0, 8.0]u"°C" # minimum air temperatures (°C)
air_temperature_max = [30.0, 25]u"°C" # maximum air temperatures (°C)
humidity_min = [20.0, 30.0] # min relative humidity (%)
humidity_max = [80.0, 90.0] # max relative humidity (%)
wind_min = [0.1, 0.2]u"m/s" # min wind speed (m/s)
wind_max = [1.0, 1.4]u"m/s" # max wind speed (m/s)
cloud_min = [20.0, 23.0] # min cloud cover (%)
cloud_max = [90.0, 100.0] # max cloud cover (%)
initial_soil_moisture = [0.2, 0.2] # fractional
min_shade = 0.0
max_shade = 90.0

# run microclimate model in minshade environment
micro_minshade = runmicro(;
    latitude,
    elevation,
    heights,
    days,
    hours,
    minima_times,
    maxima_times,
    air_temperature_min,
    air_temperature_max,
    humidity_min,
    humidity_max,
    wind_min,
    wind_max,
    cloud_min,
    cloud_max,
    initial_soil_moisture,
    shades = fill(min_shade, length(days)),
)

# run microclimate model in maxshade environment
micro_maxshade = runmicro(;
    latitude,
    elevation,
    heights,
    days,
    hours,
    minima_times,
    maxima_times,
    air_temperature_min,
    air_temperature_max,
    humidity_min,
    humidity_max,
    wind_min,
    wind_max,
    cloud_min,
    cloud_max,
    initial_soil_moisture,
    shades = fill(max_shade, length(days)),
)

env_minshade = EnvironmentalVarsVec(
    T_air = K.(micro_minshade.air_temperature[:, 2]), # second column is first node above surface
    T_sky = micro_minshade.sky_temperature,
    T_sub = micro_minshade.soil_temperature[:, 1], # surface temperature
    rh = micro_minshade.relative_humidity[:, 2], # second column is first node above surface
    vel = micro_minshade.wind_speed[:, 2], # second column is first node above surface
    Q_sol = micro_minshade.global_solar .* (1.0 - min_shade / 100.0),
    Q_dir = micro_minshade.direct_solar .* (1.0 - min_shade / 100.0),
    Q_dif = micro_minshade.diffuse_solar .* (1.0 - min_shade / 100.0),
    zen = micro_minshade.zenith_angle
)

env_maxshade = EnvironmentalVarsVec(
    T_air = K.(micro_minshade.air_temperature[:, 2]), # second column is first node above surface
    T_sky = micro_maxshade.sky_temperature,
    T_sub = micro_maxshade.soil_temperature[:, 1], # surface temperature
    rh = micro_maxshade.relative_humidity[:, 2], # second column is first node above surface
    vel = micro_maxshade.wind_speed[:, 2], # second column is first node above surface
    Q_sol = micro_maxshade.global_solar .* (1.0 - max_shade / 100.0),
    Q_dir = micro_maxshade.direct_solar .* (1.0 - max_shade / 100.0),
    Q_dif = micro_maxshade.diffuse_solar .* (1.0 - max_shade / 100.0),
    zen = micro_maxshade.zenith_angle
)

# set shade
environment = env_minshade

# compute body temperature
n = length(days) * (length(hours) - 1)
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
plot!(1:1:n, environment.T_air)

plot(1:1:n, u"mg/hr".(evap_out.m_cut))
plot!(1:1:n, u"mg/hr".(evap_out.m_resp))
plot!(1:1:n, u"mg/hr".(evap_out.m_eyes))
plot!(1:1:n, u"mg/hr".(evap_out.m_evap))
