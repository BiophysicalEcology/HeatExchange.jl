using Revise
using Unitful
using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

include("./Ectotherm.jl/src/GEOM.jl")
include("./Ectotherm.jl/src/biophysics.jl")
include("./Ectotherm.jl/src/heat_balance.jl")

T_air = (20+273.15)K
T_surf = (25.1+273.15)K
T_core = (20.2+273.15)K
F_sky = 0.4
F_sub = 0.4
F_obj = 0
ϵ_org = 0.95
ϵ_sub = 0.95
ϵ_sky = 0.8
T_sky = (10+273.15)K
T_sub = (30+273.15)K
P_atmos = 101325Pa
ψ_org = -7.07 * 100J/kg
SKINW = 0.1
AEFF = 1.192505e-05m^2
A_tot = 0.01325006m^2
hd = 0.02522706m/s
PEYES = 0.03 / 100
rh = 50
elev = 0m
vel = 1m/s
fluid = 0
p_cond = 0.1
p_cont = 0
A_v = A_tot * p_cond
A_t = A_tot * p_cont

RADIN = radin(A_tot, A_v, A_t, F_sky, F_sub, F_obj, ϵ_org, ϵ_sub, ϵ_sky, T_sky, T_sub)
RADOUT = radout(T_surf, A_tot, A_v, A_t, F_sky, F_sub, ϵ_org)

GEVAP = 1.177235e-09kg/s
T_air = (20+273.15)K
T_surf = (25.1+273.15)K
T_core = (25+273.15)K
EVAP = evap(T_core, T_surf, GEVAP, ψ_org, SKINW, AEFF, A_tot, hd, PEYES, T_air, rh, vel, P_atmos)

density = 1000kg/m^3
trunkmass = 0.04kg
trunkshapeb = 3
trunkshape = Cylinder(trunkmass, density, trunkshapeb) 
trunk = Body(trunkshape, Naked())
T_surf = (25+273.15)K
A_v = 0.01192505
CONV = convection(trunk, A_v, T_air, T_surf, vel, P_atmos, elev, fluid)

# heat_balance(trunk, params, env)