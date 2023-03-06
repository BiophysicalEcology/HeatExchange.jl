using Revise
using Unitful
using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R
using Roots

#include("./Ectotherm.jl/src/geometry.jl")
#include("./Ectotherm.jl/src/biophysics.jl")
#include("./Ectotherm.jl/src/heat_balance.jl")

# environment
T_air = (20+273.15)K
T_sky = (10+273.15)K
T_sub = (30+273.15)K
P_atmos = 101325Pa
rh = 50
elev = 0m
vel = 1m/s
fluid = 0
Z = 0°
Q_dir = 700W/m^2
Q_dif = 150W/m^2
Q_norm = Q_dir / cos(Z) # use this as Q_dir if want organism to be orienting towards beam
α_sub = 0.85
ϵ_sub = 1
ϵ_sky = 1
pctO2 = 20.95
pctCO2 = 0.03
pctN2 = 79.0

# organism geometry
mass_organism = 0.04kg
ρ_organism = 1000kg/m^3
shapeb_organism = 2
shape_organism = Cylinder(mass_organism, ρ_organism, shapeb_organism) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
body_organism = Body(shape_organism, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct

A_tot = body_organism.geometry.area

T_core = (20+273.15)K
T_surf = (20+273.15)K # skin
F_sky = 0.4
F_sub = 0.4
α_org_dorsal = 0.8
α_org_ventral = 0.8
ϵ_org = 0.95
ψ_org = -7.07 * 100J/kg
p_wet = 0.1/100
p_eyes = 0.03 / 100
p_cond = 0.1
p_cont = 0
O2_ext_ref = 20
pant = 1
rq = 0.8

A_v = A_tot * p_cond
A_t = A_tot * p_cont
A_sil = sil_area_of_cylinder(body_organism.geometry.lengths[2]/2, body_organism.geometry.lengths[1], Z)
A_up = A_tot / 2
A_down = A_tot / 2

Q_norm = Q_dir / cos(Z)
Q_solar = solar(α_org_dorsal, α_org_ventral, A_sil, A_up, A_down, α_sub, F_sub, F_sky, Q_dir, Q_dif)
Q_IR_in = radin(A_tot, F_sky, F_sub, ϵ_org, ϵ_sub, ϵ_sky, T_sky, T_sub)
Q_IR_out = radout(T_surf, A_tot, F_sky, F_sub, ϵ_org)
Q_metab = 0.01241022W
T_x = T_air

resp_out = resp(T_x, mass_organism, Q_metab, O2_ext_ref, pant, rq, T_air, rh, P_atmos, pctO2, pctCO2, pctN2)
M_resp = resp_out.M_resp
#M_resp = 1.177235e-09kg/s # respiratory water loss to be calculated by function resp

conv_out = convection(body_organism, A_v, T_air, T_surf, vel, P_atmos, elev, fluid)
evap_out = evap(T_core, T_surf, M_resp, ψ_org, p_wet, A_tot, conv_out.Hd, p_eyes, T_air, rh, P_atmos)

Q_conv = conv_out.Q_conv
Q_evap = evap_out.Q_evap

Q_in = Q_solar + Q_IR_in
Q_out = Q_IR_out + Q_conv + Q_evap 
Q_in - Q_out

# heat_balance(trunk, params, env)

function energy_bal(T_x)

# compute areas for exchange
A_v = A_tot * p_cond
A_t = A_tot * p_cont
A_sil = sil_area_of_cylinder(body_organism.geometry.lengths[2]/2, body_organism.geometry.lengths[1], Z)
A_up = A_tot / 2
A_down = A_tot / 2

Q_solar = solar(α_org_dorsal, α_org_ventral, A_sil, A_up, A_down, α_sub, F_sub, F_sky, Q_dir, Q_dif)
Q_IR_in = radin(A_tot, F_sky, F_sub, ϵ_org, ϵ_sub, ϵ_sky, T_sky, T_sub)
Q_IR_out = radout(T_x, A_tot, F_sky, F_sub, ϵ_org)

conv_out = convection(body_organism, A_v, T_air, T_x, vel, P_atmos, elev, fluid)
evap_out = evap(T_x, T_x, J_resp, ψ_org, p_wet, A_tot, conv_out.Hd, p_eyes, T_air, rh, P_atmos)

Q_conv = conv_out.Q_conv # convective heat loss
Q_evap = evap_out.Q_evap # evaporative heat loss

Q_in = Q_solar + Q_IR_in # energy in
Q_out = Q_IR_out + Q_conv + Q_evap # energy out
Q_in - Q_out # this must balance

end;

T_x = T_air

T_surf = find_zero(energy_bal, (T_air - 10K, T_air + 100K), Bisection())
T_surf_C = (Unitful.ustrip(T_surf) - 273.15)°C