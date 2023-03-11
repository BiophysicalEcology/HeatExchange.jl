using Revise
using ModelParameters
using Unitful
using Unitful: °, rad, °C, K, Pa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, R
using Roots

include("../src/geometry.jl")
include("../src/organism.jl")
include("../src/environment.jl")
include("../src/biophysics.jl")
include("../src/heat_balance.jl")

# environment
T_air = K(20°C)
T_sky = K(-5°C)
T_sub = K(30°C)
k_sub = 0.5W/m/K
P_atmos = 101325Pa
rh = 5
elev = 0m
vel = 1m/s
fluid = 0
Z = 20°
Q_sol = 1000W/m^2
Q_dif = Q_sol * 0.1
Q_norm = Q_sol / cos(Z) # use this as Q_dir if want organism to be orienting towards beam
Q_dir = Q_norm - Q_dif
α_sub = 0.2
ϵ_sub = 1
ϵ_sky = 1
fO2 = 0.2095
fCO2 = 0.0003
fN2 = 0.79

# organism morphology
mass_organism = 0.04kg
ρ_organism = 1000kg/m^3
shapeb_organism = 3
shapec_organism = 2 / 3
#shape_organism = Cylinder(mass_organism, ρ_organism, shapeb_organism) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
shape_organism = Ellipsoid(mass_organism, ρ_organism, shapeb_organism, shapec_organism) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
body_organism = Body(shape_organism, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct
p_eyes = 0 / 100
p_cond = 0.1
p_cont = 0
A_tot = body_organism.geometry.area
A_cond = A_tot * p_cond
A_conv = A_tot * (1 - p_cond)
A_sil = calc_silhouette_area(body_organism, Z)
A_up = A_tot / 2
A_down = A_tot / 2
F_sky = 0.4
F_sub = 0.4
α_org_dorsal = 0.85
α_org_ventral = 0.85
ϵ_org_dorsal = 0.95
ϵ_org_ventral = 0.95
Le = 0.025m

# organism physiology
k_body = 0.5W/m/K
ψ_org = -7.07 * 100J/kg
p_wet = 0.1/100
fO2_extract = 0.20
pant = 1.0
rq = 0.8
M1 = 0.013
M2 = 0.8
M3 = 0.038

# organism state
T_core = K(35.35236°C)

# calculate heat fluxes
Q_metab = metabolism(mass_organism, T_core, M1, M2, M3)
resp_out = respiration(T_core, Q_metab, fO2_extract, pant, rq, T_air, rh, elev, P_atmos, fO2, fCO2, fN2)
Q_resp = resp_out.Q_resp
Q_gen_net = Q_metab - Q_resp
Q_gen_spec = Q_gen_net / body_organism.geometry.volume
Tsurf_Tlung_out = get_Tsurf_Tlung(body_organism, k_body, Q_gen_spec, T_core)
T_surf = Tsurf_Tlung_out.T_surf
T_lung = Tsurf_Tlung_out.T_lung
Q_solar = solar(α_org_dorsal, α_org_ventral, A_sil, A_tot, A_cond, F_sub, F_sky, α_sub, Q_sol, Q_dir, Q_dif)
Q_IR_in = radin(A_tot, A_cond, F_sky, F_sub, ϵ_org_dorsal, ϵ_org_ventral, ϵ_sub, ϵ_sky, T_sky, T_sub)
Q_IR_out = radout(T_surf, A_tot, A_cond, F_sky, F_sub, ϵ_org_dorsal, ϵ_org_ventral)
Q_metab = metabolism(body_organism.shape.mass, T_core, M1, M2, M3)
Q_cond = conduction(A_cond, Le, T_surf, T_sub, k_sub)
conv_out = convection(body_organism, A_conv, T_air, T_surf, vel, P_atmos, elev, fluid)
evap_out = evaporation(T_core, T_surf, resp_out.m_resp, ψ_org, p_wet, A_conv, conv_out.Hd, p_eyes, T_air, rh, P_atmos)
Q_conv = conv_out.Q_conv
Q_evap = evap_out.Q_evap
Q_resp = resp_out.Q_resp

Q_in = Q_solar + Q_IR_in + Q_metab
Q_out = Q_IR_out + Q_conv + Q_evap + Q_resp + Q_cond
Q_in - Q_out


T_surf_s = find_zero(heat_balance, (T_air - 40K, T_air + 100K), Bisection())
T_surf_C = (Unitful.ustrip(T_surf_s) - 273.15)°C


# using structs to pass parameters

# define the geometry
mass_organism = 0.04kg
ρ_organism = 1000kg/m^3
#shapeb_organism = 2
shapeb_organism = 3
shapec_organism = 2 / 3
#shape_organism = Cylinder(mass_organism, ρ_organism, shapeb_organism)
#body_organism = Body(shape_organism, Naked())
shape_organism = Ellipsoid(mass_organism, ρ_organism, shapeb_organism, shapec_organism) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
body_organism = Body(shape_organism, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct

# construct the Model which holds the parameters of the organism in the Organism concrete struct, of type AbstractOrganism
lizard = Model(Organism(body_organism, OrganismParams()))
# get the environmental parameters
environmental_params = EnvironmentalParams()
# get the variables for both the organism and environment
variables = (organism=OrganismalVars(), environment=EnvironmentalVars())

# define the method 'heat_balance' for passing to find_zero, which dispatches off 'lizard' 
T_air = EnvironmentalVars().T_air
heat_balance(T_air, lizard, environmental_params, variables)
T_surf_s = find_zero(t -> heat_balance(t, lizard, environmental_params, variables), (T_air - 40K, T_air + 100K), Bisection())
T_surf_C = (Unitful.ustrip(T_surf_s) - 273.15)°C
