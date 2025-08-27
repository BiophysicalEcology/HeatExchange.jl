using HeatExchange
using ModelParameters
using Unitful, UnitfulMoles
using Unitful: °, rad, °C, K, Pa, J, kJ, W, ml, L, g, kg, cm, m, s, hr, d, mol, R
using Roots
using Test

# environment
T_air = K(20.0°C)
T_sky = K(-5.0°C)
T_sub = K(30.0°C)
k_sub = 0.5W/m/K
P_atmos = 101325.0Pa
rh = 5.0
elev = 0.0m
vel = 1.0m/s
fluid = 0
Z = 20.0°
Q_sol = 1000.0W/m^2
Q_dif = Q_sol * 0.1
Q_norm = Q_sol / cos(Z) # use this in calculating Q_dir if want organism to be orienting towards beam
Q_dir = Q_norm - Q_dif
α_sub = 0.2
ϵ_sub = 1.0
ϵ_sky = 1.0
fO2 = 0.2095
fCO2 = 0.0003
fN2 = 0.79

# organism morphology
mass = 0.04kg
ρ_body = 1000kg/m^3
shapeb = 3
shapec = 2 / 3
#shape_body = Cylinder(mass, ρ_body, shapeb) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
shape_body = Ellipsoid(mass, ρ_body, shapeb, shapec) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
geometric_traits = Body(shape_body, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct
p_eyes = 0. / 100.
p_cond = 0.1
p_cont = 0.
A_tot = geometric_traits.geometry.area
A_cond = A_tot * p_cond
A_conv = A_tot * (1 - p_cond)
A_sil = calc_silhouette_area(geometric_traits, Z)
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
T_core = K(35.497186492935384°C)

# calculate heat fluxes
metab_out = metabolism(mass, T_core, M1, M2, M3)
Q_metab = metab_out.Q_metab
resp_out = respiration(T_core, Q_metab, fO2_extract, pant, rq, T_air, rh, elev, P_atmos, fO2, fCO2, fN2)
Q_resp = resp_out.Q_resp
Q_gen_net = Q_metab - Q_resp
Q_gen_spec = Q_gen_net / geometric_traits.geometry.volume
Tsurf_Tlung_out = get_Tsurf_Tlung(geometric_traits, k_body, Q_gen_spec, T_core)
T_surf = Tsurf_Tlung_out.T_surf
T_lung = Tsurf_Tlung_out.T_lung
solar_out = solar(α_org_dorsal, α_org_ventral, A_sil, A_tot, A_cond, F_sub, F_sky, α_sub, Q_sol, Q_dir, Q_dif)
Q_solar = solar_out.Q_solar
ir_gain = radin(A_tot, A_cond, F_sky, F_sub, ϵ_org_dorsal, ϵ_org_ventral, ϵ_sub, ϵ_sky, T_sky, T_sub)
Q_ir_in = ir_gain.Q_ir_in
ir_loss = radout(T_surf, A_tot, A_cond, F_sky, F_sub, ϵ_org_dorsal, ϵ_org_ventral)
Q_ir_out = ir_loss.Q_ir_out
Q_cond = conduction(A_cond, Le, T_surf, T_sub, k_sub)
conv_out = convection(geometric_traits, A_conv, T_air, T_surf, vel, P_atmos, elev, fluid, fO2, fCO2, fN2)
evap_out = evaporation(T_core, T_surf, resp_out.m_resp, ψ_org, p_wet, A_conv, conv_out.Hd, p_eyes, T_air, rh, elev, P_atmos, fO2, fCO2, fN2)
Q_conv = conv_out.Q_conv
Q_evap = evap_out.Q_evap
Q_resp = resp_out.Q_resp

Q_in = Q_solar + Q_ir_in + Q_metab
Q_out = Q_ir_out + Q_conv + Q_evap + Q_resp + Q_cond
Q_in - Q_out

# This is broken
# (For it to work we either need an anonymous function 
# or a functor struct - you can make a struct into a function!
# then we have all the parameters attached to it.
@test_broken T_core_s = find_zero(heat_balance, (T_air - 40K, T_air + 100K), Bisection())
@test_broken T_core_C = (Unitful.ustrip(T_core_s) - 273.15)°C

# using structs to pass parameters

# define the geometry
mass = 0.04kg
ρ_body = 1000kg/m^3
#shapeb = 2
shapeb = 3
shapec = 2 / 3
#shape_body = Cylinder(mass, ρ_body, shapeb)
#geometric_traits = Body(shape_body, Naked())
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
heat_balance(T_air, lizard, environmental_params, variables)
T_core_s = find_zero(t -> heat_balance(t, lizard, environmental_params, variables), (T_air - 40K, T_air + 100K), Bisection())
T_core_C = (Unitful.ustrip(T_core_s) - 273.15)°C
heat_balance_out = heat_balance(T_core_s, lizard, environmental_params, variables)
