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
# get the variables for both the organism and environment
variables = (organism=OrganismalVars(), environment=EnvironmentalVars())

# define the method 'heat_balance' for passing to find_zero, which dispatches off 'lizard' 
T_air = EnvironmentalVars().T_air
T_core_s = find_zero(t -> heat_balance(t, lizard, environmental_params, variables), (T_air - 40K, T_air + 100K), Bisection())
T_core_C = (Unitful.ustrip(T_core_s) - 273.15)°C
heat_balance_out = heat_balance(T_core_s, lizard, environmental_params, variables)

# specify place and time
lat = -30.0°
lon = 140.0°
elev = 10.0m
days = [15, 45]*1.0
hours = collect(0.:1:24.)

# compute solar radiation
solrad_out = solrad(;
    days,       # days of year
    hours,      # hours of day
    lat,        # latitude (degrees)
    elev,
    )
solrad_out.Zenith[solrad_out.Zenith.>90u"°"] .= 90u"°"
Q_sol = solrad_out.Global
Zenith = solrad_out.Zenith
Q_dir = solrad_out.Direct
Q_dif = solrad_out.Scattered

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

# interpolate air temperature to hourly
TAIRs, WNs, RHs, CLDs = hourly_vars(
    TMINN,
    TMAXX,
    WNMINN,
    WNMAXX,
    RHMINN,
    RHMAXX,
    CCMINN,
    CCMAXX,
    solrad_out,
    TIMINS,
    TIMAXS
)
RHs[RHs.>100] .= 100
CLDs[CLDs.>100] .= 100

# compute sky temperature for downwelling longwave
Tskys = zeros(length(TAIRs))u"K"  # or whatever you have from your loop
for i in 1:length(TAIRs)
    Tskys[i] = Microclimate.get_longwave(
        elev=elev,
        rh=RHs[i],
        tair=TAIRs[i],
        tsurf=TAIRs[i],
        slep=0.95,
        sle=0.95,
        cloud=CLDs[i],
        viewf=1,
        shade=0.0
    ).Tsky
end

env_vec = EnvironmentalVarsVec(
    #T_air = (collect(15.0:5.0:35.0) .+ 273.15) .* 1.0K,
    #T_air = fill(293.15K, length(Q_sol)),
    T_air = K.(TAIRs),
    T_sky = Tskys,
    T_sub = K.(TAIRs),
    rh = RHs,
    vel = WNs,
    Q_sol = Q_sol,
    Q_dir = Q_dir,
    Q_dif = Q_dif,
    zen = Zenith
)


n = length(env_vec.T_air)
T_core_s = Vector{typeof(env_vec.T_air[1])}(undef, n)   # in Kelvin
T_core_C = Vector{typeof(0.0°C)}(undef, n)
heat_balance_out = Vector{Float64}(undef, n)  # or whatever type heat_balance returns

for i in 1:n
    # Construct a single-environment object for this iteration
    env_i = EnvironmentalVars(
        T_air   = env_vec.T_air[i],
        T_sky   = env_vec.T_sky[i],
        T_sub   = env_vec.T_sub[i],
        rh      = env_vec.rh[i],
        vel     = env_vec.vel[i],
        P_atmos = env_vec.P_atmos[i],
        zen     = env_vec.zen[i],
        k_sub   = env_vec.k_sub[i],
        Q_sol   = env_vec.Q_sol[i],
        Q_dir   = env_vec.Q_dir[i],
        Q_dif   = env_vec.Q_dif[i]
    )

    # Variables tuple as before
    variables_i = (organism = OrganismalVars(), environment = env_i)

    # Root-finding for this environment
    f(T_core) = heat_balance(T_core, lizard, environmental_params, variables_i)

    T_core_s[i] = find_zero(f, (env_i.T_air - 40K, env_i.T_air + 100K), Bisection())

    # Convert to °C
    T_core_C[i] = (Unitful.ustrip(T_core_s[i]) - 273.15)°C

    # Store heat_balance output if needed
    #heat_balance_out[i] = f(T_core_s[i])
end
T_core_C

plot(1:1:length(T_core_C), T_core_C)
plot!(1:1:length(T_core_C), TAIRs)