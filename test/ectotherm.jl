using HeatExchange
using FluidProperties
using ModelParameters
using Unitful, UnitfulMoles
using Roots
using Test
using DataFrames, CSV

testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

ecto_input_vec = DataFrame(CSV.File("$testdir/data/ectoR_input.csv"))[:, 2]
ecto_output_vec = DataFrame(CSV.File("$testdir/data/ectoR_output.csv"))[:, 2]

names = [
    :Ww_g, :alpha, :epsilon, :rho_body, :fatosk, :fatosb, :shape, :shape_a, :shape_b, :shape_c, :conv_enhance, 
    :custom_shape1, :custom_shape2, :custom_shape3, :custom_shape4, :custom_shape5, :custom_shape6, 
    :custom_shape7, :custom_shape8, :pct_cond, :pct_touch, :postur, :orient, :k_flesh, 
    :M_1, :M_2, :M_3, :M_4, :Q_act, :pct_wet, :pct_eyes, :pct_mouth, :psi_body, :pantmax, 
    :F_O2, :RQ, :delta_air, :leaf, :g_vs_ab, :g_vs_ad, :elevation, :alpha_sub, 
    :epsilon_sub, :epsilon_sky, :pres, :fluid, :O2gas, :CO2gas, :N2gas, :K_sub, 
    :PDIF, :SHADE, :QSOLR, :Z, :TA, :TGRD, :TSUBST, :TSKY, :VEL, :RH
]

ecto_input = (; zip(names, ecto_input_vec)...)

names = [
:TC, :TSKIN, :TLUNG, :QSOL, :QIRIN, :QMET, :QRESP, :QEVAP, :QIROUT, :QCONV, :QCOND, :ENB, 
:O2_ml, :H2OResp_g, :H2OCut_g, :H2OEyes_g, :AREA, :AV, :AT, :AL, :VOL, :R, :R1, 
:ASEMAJR, :BSEMINR, :CSEMINR, :ASILN, :ASILP, :AEFF, :QSOLAR, :QSDIFF, :QSRSB, :QSSKY, :QSOBJ, :QRESP2, 
:GEVAP, :AIRML1, :AIRML2, :WMOL1, :WMOL2, :O2MOL1, :O2MOL2, :CO2MOL1, :CO2MOL2, :CONVAR, :QCONV2, 
:HC, :HD, :SH, :QFREE, :HCFREE, :HCFORC, :SHFREE, :SHFORC, :HDFREE, :HDFORC, :QSEVAP, :WEVAP, 
:WRESP, :WCUT, :WEYES, :G_V
]

ecto_output = (; zip(names, ecto_output_vec)...)


# environment
T_air = u"K"((ecto_input.TA)u"°C")
T_sky = u"K"((ecto_input.TSKY)u"°C")
T_substrate = u"K"((ecto_input.TSUBST)u"°C")
k_substrate = (ecto_input.K_sub)u"W/m/K"
P_atmos = (ecto_input.pres)u"Pa"
rh = (ecto_input.RH/100)
elevation = (ecto_input.elevation)u"m"
wind_speed = (ecto_input.VEL)u"m/s"
fluid = ecto_input.fluid
zenith_angle = (ecto_input.Z)u"°"
solar_radiation = (ecto_input.QSOLR)u"W/m^2"
α_substrate = ecto_input.alpha_sub
shade = ecto_input.SHADE
ϵ_substrate = ecto_input.epsilon_sub
ϵ_sky = ecto_input.epsilon_sky
fO2 = ecto_input.O2gas / 100
fCO2 = ecto_input.CO2gas / 100
fN2 = ecto_input.N2gas / 100

# organism morphology
mass = u"kg"((ecto_input.Ww_g)u"g")
ρ_flesh = (ecto_input.rho_body)u"kg/m^3"
shape_b = ecto_input.shape_b
shape_c = ecto_input.shape_c
#shape_body = Cylinder(mass, ρ_flesh, shape_b) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
eye_fraction = ecto_input.pct_eyes / 100
conduction_fraction = ecto_input.pct_cond / 100

F_sky = ecto_input.fatosk
F_substrate = ecto_input.fatosb
α_body_dorsal = ecto_input.alpha
α_body_ventral = ecto_input.alpha
ϵ_body_dorsal = ecto_input.epsilon
ϵ_body_ventral = ecto_input.epsilon
Le = 0.025u"m"
ventral_fraction = 0.5

# organism physiology
k_flesh = (ecto_input.k_flesh)u"W/m/K"
ψ_org = (ecto_input.psi_body)u"J/kg"
skin_wetness = ecto_input.pct_wet / 100
fO2_extract = ecto_input.F_O2 / 100
rq = ecto_input.RQ
M1 = ecto_input.M_1
M2 = ecto_input.M_2
M3 = ecto_input.M_3
M4 = ecto_input.M_4

pant = 1.0

#shape_body = Ellipsoid(mass, ρ_flesh, shape_b, shape_c) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
shape_body = DesertIguana(mass, ρ_flesh) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
geometric_pars = Body(shape_body, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct
A_total = get_total_area(geometric_pars) #geometric_pars.geometry.area.total
A_conduction = A_total * conduction_fraction
A_convection = A_total * (1 - conduction_fraction)
A_sil_normal, A_sil_parallel = silhouette_area(shape_body)
A_sil = (A_sil_normal + A_sil_parallel) / 2
A_up = A_total / 2
A_down = A_total / 2
A_eff = A_convection * skin_wetness

# geometry test
@test A_total ≈ (ecto_output.AREA)u"m^2" rtol=1e-9
@test A_conduction ≈ (ecto_output.AV)u"m^2" rtol=1e-9
@test A_sil_normal ≈ (ecto_output.ASILN)u"m^2" rtol=1e-9
@test A_sil_parallel ≈ (ecto_output.ASILP)u"m^2" rtol=1e-9
@test A_eff ≈ (ecto_output.AEFF)u"m^2" rtol=1e-9
@test geometric_pars.geometry.characteristic_dimension ≈ (ecto_output.AL)u"m" rtol=1e-9
@test geometric_pars.geometry.length[1] ≈ (ecto_output.R1)u"m" rtol=1e-9
@test geometric_pars.geometry.volume ≈ (ecto_output.VOL)u"m^3" rtol=1e-9

# organism state
T_core = u"K"((ecto_output.TC)u"°C")

# calculate heat fluxes

# metabolism
metab_out = metabolic_rate(AndrewsPough2(), mass, T_core; M1, M2, M3, M4)

# metabolism test
@test metab_out.Q_metab ≈ (ecto_output.QMET)u"W" rtol=1e-9
@test metab_out.V_O2 ≈ (ecto_output.O2_ml)u"mL/hr" rtol=1e-9

Q_metab = metab_out.Q_metab

# respiration
resp_out = respiration_ectotherm(T_core, Q_metab, fO2_extract, pant, rq, T_air, rh, P_atmos, fO2, fCO2, fN2)
Q_resp = resp_out.Q_resp

# respiration test
@test u"mol/s"(resp_out.J_CO2_in) ≈ (ecto_output.CO2MOL1)u"mol/s" rtol=1e-6
@test u"mol/s"(resp_out.J_CO2_out) ≈ (ecto_output.CO2MOL2)u"mol/s" rtol=1e-6
@test u"mol/s"(resp_out.J_H2O_in) ≈ (ecto_output.WMOL1)u"mol/s" rtol=1e-6
@test u"mol/s"(resp_out.J_H2O_out) ≈ (ecto_output.WMOL2)u"mol/s" rtol=1e-6
@test u"mol/s"(resp_out.J_O2_in) ≈ (ecto_output.O2MOL1)u"mol/s" rtol=1e-6
@test u"mol/s"(resp_out.J_O2_out) ≈ (ecto_output.O2MOL2)u"mol/s" rtol=1e-6
@test u"mol/s"(resp_out.J_air_in) ≈ (ecto_output.AIRML1)u"mol/s" rtol=1e-6
@test u"mol/s"(resp_out.J_air_out) ≈ (ecto_output.AIRML2)u"mol/s" rtol=1e-6
@test resp_out.Q_resp ≈ (ecto_output.QRESP2)u"W" rtol=1e-3
@test resp_out.m_resp ≈ (ecto_output.GEVAP)u"g/s" rtol=1e-3

Q_resp = resp_out.Q_resp
Q_gen_net = Q_metab - Q_resp
Q_gen_spec = Q_gen_net / geometric_pars.geometry.volume

# compute skin and lung temperature
T_surface, T_lung = Tsurf_and_Tlung(geometric_pars, k_flesh, Q_gen_spec, T_core)

# test lung temperature
@test T_lung ≈ u"K"((ecto_output.TLUNG)u"°C") rtol=1e-8

# solar radiation
diffuse_radiation = solar_radiation * ecto_input.PDIF
normal_radiation = solar_radiation * (1 - ecto_input.PDIF) / cos(zenith_angle)  # use this in calculating Q_dir if want organism to be orienting towards beam
direct_radiation = normal_radiation - diffuse_radiation
solar_out = solar(α_body_dorsal, α_body_ventral, A_sil_normal, A_total, A_conduction, F_substrate, F_sky, α_substrate, shade, solar_radiation, normal_radiation, diffuse_radiation)
Q_solar = solar_out.Q_solar

# longwave radiation
ir_gain = radin(A_total, A_conduction, F_sky, F_substrate, ϵ_body_dorsal, ϵ_body_ventral, ϵ_substrate, ϵ_sky, T_sky, T_substrate)
Q_ir_in = ir_gain.Q_ir_in
ir_loss = radout(T_surface, A_total, A_conduction, F_sky, F_substrate, ϵ_body_dorsal, ϵ_body_ventral)
Q_ir_out = ir_loss.Q_ir_out

# conduction
Q_cond = conduction(A_conduction, Le, T_surface, T_substrate, k_substrate)

# convection
conv_out = convection(geometric_pars, A_convection, T_air, T_surface, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)
Q_conv = conv_out.Q_conv

# evaporation
evap_out = evaporation(; T_surface, ψ_org, wetness=skin_wetness, area=A_convection, conv_out.hd, eye_fraction, T_air, rh, P_atmos, fO2, fCO2, fN2)
Q_evap = evap_out.Q_evap


# energy balance test
@test solar_out.Q_solar ≈ (ecto_output.QSOL)u"W" rtol=1e-9
@test Q_cond ≈ (ecto_output.QCOND)u"W" rtol=1e-6
@test conv_out.Q_conv ≈ (ecto_output.QCONV)u"W" rtol=1e-4
@test evap_out.Q_evap ≈ (ecto_output.QEVAP)u"W" rtol=1e-4


Q_in = Q_solar + Q_ir_in + Q_metab
Q_out = Q_ir_out + Q_conv + Q_evap + Q_resp + Q_cond
Q_in - Q_out

# This is broken
# (For it to work we either need an anonymous function 
# or a functor struct - you can make a struct into a function!
# then we have all the parameters attached to it.
#@test_broken T_core_s = find_zero(heat_balance, (T_air - 40K, T_air + 100K), Bisection())
#@test_broken T_core_C = (Unitful.ustrip(T_core_s) - 273.15)°C

# using structs to pass parameters
bodypars = BodyPars(;
    mass,
    ρ_flesh,
    shape_b,
    shape_c,
)

integument_pars = IntegumentPars(;
    α_body_dorsal,
    α_body_ventral,
    ϵ_body_dorsal,
    ϵ_body_ventral,
    F_sky,
    F_substrate,
    eye_fraction,
    skin_wetness,
    conduction_fraction,
    ventral_fraction,
)

physio_pars = PhysioPars(;
    fO2_extract,
    rq,
    k_flesh,
)

thermoreg_pars = ThermoregulationPars(;
    pant,
)

# construct the Model which holds the parameters of the organism in the Organism concrete struct, of type AbstractOrganism
lizard = Model(Organism(geometric_pars, integument_pars, physio_pars, thermoreg_pars))

# get the environmental parameters
environmental_params = EnvironmentalPars(
    α_substrate,
    shade,
    ϵ_substrate,
    ϵ_sky,
    elevation,
    fluid,
    fN2,
    fO2,
    fCO2,
)

# get the variables for both the organism and environment
organism=OrganismalVars(
    T_core,
    T_surface,
    T_lung,
    ψ_org,
)

environment=EnvironmentalVars(
    T_air,
    T_sky,
    T_substrate,
    rh,
    wind_speed,
    P_atmos,
    zenith_angle,
    k_substrate,
    solar_radiation,
    normal_radiation,
    diffuse_radiation,
)
variables = (organism = organism, environment = environment)

# define the method 'heat_balance' for passing to find_zero, which dispatches off 'lizard' 
T_air = environment.T_air
heat_balance(T_air, lizard, environmental_params, variables)
T_core_s = find_zero(t -> heat_balance(t, lizard, environmental_params, variables), (T_air - 40u"K", T_air + 100u"K"), Bisection())
T_core_C = (Unitful.ustrip(T_core_s) - 273.15)u"°C"

# test core temperature calculation
@test T_core_C ≈ (ecto_output.TC)u"°C" rtol=1e-4

heat_balance_out = heat_balance(T_core_s, lizard, environmental_params, variables)

@test heat_balance_out.T_core ≈ (ecto_output.TC + 273.15)u"K" rtol=1e-5
@test heat_balance_out.T_surface ≈ (ecto_output.TSKIN + 273.15)u"K" rtol=1e-5
@test heat_balance_out.T_lung ≈ (ecto_output.TLUNG + 273.15)u"K" rtol=1e-5

@test heat_balance_out.enbal.Q_solar ≈ (ecto_output.QSOL)u"W" rtol=1e-9
@test heat_balance_out.enbal.Q_ir_in ≈ (ecto_output.QIRIN)u"W" rtol=1e-3 # TODO make better?
@test heat_balance_out.enbal.Q_ir_out ≈ (ecto_output.QIROUT)u"W" rtol=1e-3 # TODO make better?
@test heat_balance_out.enbal.Q_cond ≈ (ecto_output.QCOND)u"W" rtol=1e-3 # TODO make better?
@test heat_balance_out.enbal.Q_conv ≈ (ecto_output.QCONV)u"W" rtol=1e-4
@test heat_balance_out.enbal.Q_evap ≈ (ecto_output.QEVAP)u"W" rtol=1e-6
@test heat_balance_out.enbal.Q_metab ≈ (ecto_output.QMET)u"W" rtol=1e-3 # TODO make better?
@test heat_balance_out.enbal.Q_resp ≈ (ecto_output.QRESP)u"W" rtol=1e-3 # TODO make better?

@test heat_balance_out.masbal.V_O2 ≈ (ecto_output.O2_ml)u"ml/hr" rtol=1e-3 # TODO make better?
@test heat_balance_out.masbal.m_cut ≈ (ecto_output.H2OCut_g / 3600)u"g/s" rtol=1e-3 # TODO make better?
@test heat_balance_out.masbal.m_eye ≈ (ecto_output.H2OEyes_g / 3600)u"g/s" rtol=1e-3 # TODO make better?

@test u"mol/s"(heat_balance_out.resp_out.J_CO2_in) ≈ (ecto_output.CO2MOL1)u"mol/s" rtol=1e-3
@test u"mol/s"(heat_balance_out.resp_out.J_CO2_out) ≈ (ecto_output.CO2MOL2)u"mol/s" rtol=1e-3
@test u"mol/s"(heat_balance_out.resp_out.J_H2O_in) ≈ (ecto_output.WMOL1)u"mol/s" rtol=1e-3
@test u"mol/s"(heat_balance_out.resp_out.J_H2O_out) ≈ (ecto_output.WMOL2)u"mol/s" rtol=1e-3
@test u"mol/s"(heat_balance_out.resp_out.J_O2_in) ≈ (ecto_output.O2MOL1)u"mol/s" rtol=1e-3
@test u"mol/s"(heat_balance_out.resp_out.J_O2_out) ≈ (ecto_output.O2MOL2)u"mol/s" rtol=1e-3
@test u"mol/s"(heat_balance_out.resp_out.J_air_in) ≈ (ecto_output.AIRML1)u"mol/s" rtol=1e-3
@test u"mol/s"(heat_balance_out.resp_out.J_air_out) ≈ (ecto_output.AIRML2)u"mol/s" rtol=1e-3
@test heat_balance_out.resp_out.Q_resp ≈ (ecto_output.QRESP2)u"W" rtol=1e-3
@test heat_balance_out.resp_out.m_resp ≈ (ecto_output.GEVAP)u"g/s" rtol=1e-3

@test heat_balance_out.evap_out.m_cut ≈ (ecto_output.WCUT)u"g/s" rtol=1e-4 # TODO check if this can be better
@test heat_balance_out.evap_out.m_eyes ≈ (ecto_output.WEYES)u"g/s" rtol=1e-4 # TODO check if this can be better
@test heat_balance_out.evap_out.Q_evap ≈ (ecto_output.QEVAP)u"W" rtol=1e-6 # TODO check if this can be better