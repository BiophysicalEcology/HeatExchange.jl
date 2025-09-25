using HeatExchange
using Microclimate
using ModelParameters
using Unitful, UnitfulMoles
using Unitful: °, rad, °C, K, Pa, J, kJ, W, ml, L, g, kg, cm, m, s, hr, d, mol, R
using Roots
using Test
using Plots
using CSV, DataFrames

testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

Tb_NMR = (DataFrame(CSV.File("$testdir/data/TC.csv")))[:, 2] .* u"°C"
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
latitude = -30.0°
longitude = 140.0°
elevation = 10.0m
days = [15, 45]*1.0
hours = collect(0.:1:24.)
heights = [0.01]u"m"
α_sub = 0.8
albedos = [0.2, 0.2]

# set the environmental parameters
environmental_params = EnvironmentalPars(
    elev = elevation,
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

n_hours = length(air_temperature_min) * 24

minimum_shade = 0.0
maximum_shade = 0.9

# run microclimate model in minshade environment
micro_min_shade = runmicro(;
    latitude,
    elevation,
    heights,
    days,
    hours,
    air_temperature_min,
    air_temperature_max,
    humidity_min,
    humidity_max,
    wind_min,
    wind_max,
    cloud_min,
    cloud_max,
    minima_times,
    maxima_times,
    initial_soil_moisture,
    albedos,
    shades = fill(minimum_shade, n_hours),
);
plot(1:1:n, u"°C".(micro_min_shade.soil_temperature[:, 1]))

# run microclimate model in minshade environment
micro_max_shade = runmicro(;
    latitude,
    elevation,
    heights,
    days,
    hours,
    air_temperature_min,
    air_temperature_max,
    humidity_min,
    humidity_max,
    wind_min,
    wind_max,
    cloud_min,
    cloud_max,
    minima_times,
    maxima_times,
    initial_soil_moisture,
    albedos,
    shades = fill(maximum_shade, n_hours),
);

min_shade_habitat = EnvironmentalVarsVec(
    T_air = K.(micro_min_shade.air_temperature), # second column is first node above surface
    T_sky = micro_min_shade.sky_temperature,
    T_sub = micro_min_shade.soil_temperature, # surface temperature
    rh = Matrix(micro_min_shade.relative_humidity), # second column is first node above surface
    vel = Matrix(micro_min_shade.wind_speed), # second column is first node above surface
    P_atmos = fill(101325.0Pa, n_hours),
    k_sub = micro_min_shade.soil_thermal_conductivity,
    Q_sol = micro_min_shade.global_solar .* (1.0 - minimum_shade),
    Q_dir = micro_min_shade.direct_solar .* (1.0 - minimum_shade),
    Q_dif = micro_min_shade.diffuse_solar .* (1.0 - minimum_shade),
    zen = micro_min_shade.zenith_angle
);

max_shade_habitat = EnvironmentalVarsVec(
    T_air = K.(micro_max_shade.air_temperature), # second column is first node above surface
    T_sky = micro_max_shade.sky_temperature,
    T_sub = micro_max_shade.soil_temperature, # surface temperature
    rh = Matrix(micro_max_shade.relative_humidity), # second column is first node above surface
    vel = Matrix(micro_max_shade.wind_speed), # second column is first node above surface
    P_atmos = fill(101325.0Pa, n_hours),
    k_sub = micro_max_shade.soil_thermal_conductivity,
    Q_sol = micro_max_shade.global_solar .* (1.0 - maximum_shade),
    Q_dir = micro_max_shade.direct_solar .* (1.0 - maximum_shade),
    Q_dif = micro_max_shade.diffuse_solar .* (1.0 - maximum_shade),
    zen = micro_max_shade.zenith_angle
);

# compute body temperature
shade = 0.0
depth = 1
height = 2
n = n_hours
deep_rh = 99.
deep_vel = 0.01u"m/s"

balances = map(1:n) do i
    P_atmos = min_shade_habitat.P_atmos[i]
    zen     = min_shade_habitat.zen[i]
    if depth > 1
        T_air   = min_shade_habitat.T_sub[i, depth] * (1 - shade) + max_shade_habitat.T_sub[i, depth] * shade
        T_sky   = min_shade_habitat.T_sub[i, depth] * (1 - shade) + max_shade_habitat.T_sub[i, depth] * shade
        T_sub   = min_shade_habitat.T_sub[i, depth] * (1 - shade) + max_shade_habitat.T_sub[i, depth] * shade
        rh      = deep_rh
        vel     = deep_vel
        k_sub   = min_shade_habitat.k_sub[i, depth] * (1 - shade) + max_shade_habitat.k_sub[i, depth] * shade
        Q_sol   = 0.0u"W/m^2"
        Q_dir   = 0.0u"W/m^2"
        Q_dif   = 0.0u"W/m^2"
    else
        T_air   = min_shade_habitat.T_air[i, height + 1] * (1 - shade) + max_shade_habitat.T_air[i, height + 1] * shade
        T_sky   = min_shade_habitat.T_sky[i] * (1 - shade) + max_shade_habitat.T_sky[i] * shade
        T_sub   = min_shade_habitat.T_sub[i, depth] * (1 - shade) + max_shade_habitat.T_sub[i, depth] * shade
        rh      = min_shade_habitat.rh[i, height + 1] * (1 - shade) + max_shade_habitat.rh[i, height + 1] * shade
        vel     = min_shade_habitat.vel[i, height + 1] * (1 - shade) + max_shade_habitat.vel[i, height + 1] * shade
        k_sub   = min_shade_habitat.k_sub[i, depth] * (1 - shade) + max_shade_habitat.k_sub[i, depth] * shade
        Q_sol   = min_shade_habitat.Q_sol[i] * (1 - shade) + max_shade_habitat.Q_sol[i] * shade
        Q_dir   = min_shade_habitat.Q_dir[i] * (1 - shade) + max_shade_habitat.Q_dir[i] * shade
        Q_dif   = min_shade_habitat.Q_dif[i] * (1 - shade) + max_shade_habitat.Q_dif[i] * shade
    end
    env_i = EnvironmentalVars(;
        T_air,
        T_sky,
        T_sub,
        rh,
        vel,
        P_atmos,
        zen,
        k_sub,
        Q_sol,
        Q_dir,
        Q_dif,
    )
    variables_i = (organism = OrganismalVars(), environment = env_i)
    get_Tb(lizard, environmental_params, variables_i);
end
balance_out = flip2vectors(balances); # pull out each output as a vector
resp_out = flip2vectors(balance_out.resp_out); # pull out each output as a vector
evap_out = flip2vectors(balance_out.evap_out); # pull out each output as a vector
conv_out = flip2vectors(balance_out.conv_out); # pull out each output as a vector

plot(1:1:n, °C.(balance_out.T_core), ylims=[-10.0, 55.0])
plot!(1:1:n, micro_min_shade.air_temperature[:, 2])
plot!(1:1:n, micro_min_shade.soil_temperature[:, 1])
plot!(1:1:n, °C.(micro_min_shade.sky_temperature[:, 1]))

function set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
    P_atmos = min_shade_habitat.P_atmos[i]
    zen     = min_shade_habitat.zen[i]
    if depth > 1
        T_air = min_shade_habitat.T_sub[i, depth] * (1 - shade) + max_shade_habitat.T_sub[i, depth] * shade
        T_sky = min_shade_habitat.T_sub[i, depth] * (1 - shade) + max_shade_habitat.T_sub[i, depth] * shade
        T_sub = min_shade_habitat.T_sub[i, depth] * (1 - shade) + max_shade_habitat.T_sub[i, depth] * shade
        rh = deep_rh
        vel = deep_vel
        k_sub = min_shade_habitat.k_sub[i, depth] * (1 - shade) + max_shade_habitat.k_sub[i, depth] * shade
        Q_sol = 0.0u"W/m^2"
        Q_dir = 0.0u"W/m^2"
        Q_dif = 0.0u"W/m^2"
    else
        T_air = min_shade_habitat.T_air[i, height+1] * (1 - shade) + max_shade_habitat.T_air[i, height+1] * shade
        T_sky = min_shade_habitat.T_sky[i] * (1 - shade) + max_shade_habitat.T_sky[i] * shade
        T_sub = min_shade_habitat.T_sub[i, depth] * (1 - shade) + max_shade_habitat.T_sub[i, depth] * shade
        rh = min_shade_habitat.rh[i, height+1] * (1 - shade) + max_shade_habitat.rh[i, height+1] * shade
        vel = min_shade_habitat.vel[i, height+1] * (1 - shade) + max_shade_habitat.vel[i, height+1] * shade
        k_sub = min_shade_habitat.k_sub[i, depth] * (1 - shade) + max_shade_habitat.k_sub[i, depth] * shade
        Q_sol = min_shade_habitat.Q_sol[i] * (1 - shade) + max_shade_habitat.Q_sol[i] * shade
        Q_dir = min_shade_habitat.Q_dir[i] * (1 - shade) + max_shade_habitat.Q_dir[i] * shade
        Q_dif = min_shade_habitat.Q_dif[i] * (1 - shade) + max_shade_habitat.Q_dif[i] * shade
    end
    env_i = EnvironmentalVars(;
        T_air,
        T_sky,
        T_sub,
        rh,
        vel,
        P_atmos,
        zen,
        k_sub,
        Q_sol,
        Q_dir,
        Q_dif,
    )
    return env_i
end

T_forage_min = u"K"(24.5u"°C")
T_forage_max = u"K"(34.5u"°C")
T_lethal_min = u"K"(0.1u"°C")
shade = 0.0
Δshade = 0.01
min_depth = 2
max_height = 2
depth = 1
height = 1
climber = false
burrower = false
shadeseeker = false

Tb = nothing
balances = map(1:n) do i
    #i=12
    height = 1 # reset height
    shade = minimum_shade # reset shade
    depth = 1 # reset depth to surface
    env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
    variables_i = (organism=OrganismalVars(), environment=env_i)
    Tb = get_Tb(lizard, environmental_params, variables_i).T_core;
    # check if too hot and, if so, seek shade
    if shadeseeker
    while Tb > T_forage_max && shade <= maximum_shade && depth == 1
        shade += Δshade
        env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
        variables_i = (organism=OrganismalVars(), environment=env_i)
        Tb = get_Tb(lizard, environmental_params, variables_i).T_core;
    end
    shade = clamp(shade, minimum_shade, maximum_shade)
    end
    if climber
    while Tb > T_forage_max && (shade >= maximum_shade || !shadeseeker) && height < max_height
        shade = minimum_shade # reset to minimum shade burrow
        height = height + 1
        env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
        variables_i = (organism=OrganismalVars(), environment=env_i)
        Tb = get_Tb(lizard, environmental_params, variables_i).T_core;
        while (Tb > T_forage_max || Tb < T_lethal_min) && height < max_height
            height = height + 1
            env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
            variables_i = (organism=OrganismalVars(), environment=env_i)
            Tb = get_Tb(lizard, environmental_params, variables_i).T_core;
        end
    end
    end
    # check if too cold or too hot and shade maxed out, and, if so, go underground
    if burrower
    while (Tb < T_forage_min || (Tb > T_forage_max && shade >= maximum_shade) || (Tb > T_forage_max && height >= max_height)) && depth == 1
        shade = minimum_shade # reset to minimum shade burrow
        depth = max(min_depth, depth + 1)
        env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
        variables_i = (organism=OrganismalVars(), environment=env_i)
        Tb = get_Tb(lizard, environmental_params, variables_i).T_core;
        while (Tb > T_forage_max || Tb < T_lethal_min) && depth < 10
            depth = max(min_depth, depth + 1)
            env_i = set_environment(; i, shade, depth, height, min_shade_habitat, max_shade_habitat)
            variables_i = (organism=OrganismalVars(), environment=env_i)
            Tb = get_Tb(lizard, environmental_params, variables_i).T_core;
        end
    end
    end
    return get_Tb(lizard, environmental_params, variables_i);
end
balance_out = flip2vectors(balances); # pull out each output as a vector

plot(1:1:n, °C.(balance_out.T_core), ylims=[0, 55.0])
plot!(1:1:n, Tb_NMR;
        xlabel="time", ylabel="body temperature", lw=2,
        linestyle=:dash, linecolor="grey"
    )