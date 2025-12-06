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
T_ground = T_substrate
k_substrate = (ecto_input.K_sub)u"W/m/K"
P_atmos = (ecto_input.pres)u"Pa"
rh = (ecto_input.RH/100)
elevation = (ecto_input.elevation)u"m"
wind_speed = (ecto_input.VEL)u"m/s"
fluid = ecto_input.fluid
zenith_angle = (ecto_input.Z)u"°"
global_radiation = (ecto_input.QSOLR)u"W/m^2"
α_ground = ecto_input.alpha_sub
shade = ecto_input.SHADE
ϵ_ground = ecto_input.epsilon_sub
ϵ_sky = ecto_input.epsilon_sky
fO2 = ecto_input.O2gas / 100
fCO2 = ecto_input.CO2gas / 100
fN2 = ecto_input.N2gas / 100
convection_enhancement = ecto_input.conv_enhance

# organism morphology
mass = u"kg"((ecto_input.Ww_g)u"g")
ρ_flesh = (ecto_input.rho_body)u"kg/m^3"
shape_b = ecto_input.shape_b
shape_c = ecto_input.shape_c
eye_fraction = ecto_input.pct_eyes / 100
conduction_fraction = ecto_input.pct_cond / 100

F_sky = ecto_input.fatosk
F_ground = ecto_input.fatosb
α_body_dorsal = ecto_input.alpha
α_body_ventral = ecto_input.alpha
ϵ_body_dorsal = ecto_input.epsilon
ϵ_body_ventral = ecto_input.epsilon
Le = 0.025u"m"
#conduction_depth::D =    Param(2.5u"cm") this should be an environmental variable

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

# organism behaviour
if ecto_input.postur == 0.0
    solar_orientation = Intermediate()
elseif ecto_input.postur == 1.0
    solar_orientation = NormalToSun()
else
    solar_orientation = ParallelToSun()
end
pant = 1.0

# geometry
shape_pars = DesertIguana(mass, ρ_flesh) 
#shape_pars = Cylinder(mass, ρ_flesh, shape_b)
geometry = Body(shape_pars, Naked())
A_total = get_total_area(geometry)
A_dorsal = A_total * (1 - ventral_fraction)
A_ventral = A_total * ventral_fraction * (1 - conduction_fraction)
A_conduction = A_total * conduction_fraction
A_convection = A_total * (1 - conduction_fraction)
A_sil_normal = silhouette_area(shape_pars, NormalToSun())
A_sil_parallel = silhouette_area(shape_pars, ParallelToSun())
A_silhouette = silhouette_area(shape_pars, solar_orientation)
A_eff = A_convection * skin_wetness

# geometry test
@test A_total ≈ (ecto_output.AREA)u"m^2" rtol=1e-9
@test A_conduction ≈ (ecto_output.AV)u"m^2" rtol=1e-9
@test A_sil_normal ≈ (ecto_output.ASILN)u"m^2" rtol=1e-9
@test A_sil_parallel ≈ (ecto_output.ASILP)u"m^2" rtol=1e-9
@test A_eff ≈ (ecto_output.AEFF)u"m^2" rtol=1e-9
@test geometry.geometry.characteristic_dimension ≈ (ecto_output.AL)u"m" rtol=1e-9
@test geometry.geometry.length[1] ≈ (ecto_output.R1)u"m" rtol=1e-9
@test geometry.geometry.volume ≈ (ecto_output.VOL)u"m^3" rtol=1e-9

# organism state
T_core = u"K"((ecto_output.TC)u"°C")

# calculate heat fluxes

# metabolism|
Q_metab = metabolic_rate(AndrewsPough2(), mass, T_core; M1, M2, M3, M4, O2conversion=Typical())

# metabolism test
@test Q_metab ≈ (ecto_output.QMET)u"W" rtol=1e-9

# respiration
resp_out = respiration(; 
    T_lung = T_core, 
    Q_metab, 
    fO2_extract, 
    pant, 
    rq,
    mass, 
    T_air, 
    rh, 
    P_atmos, 
    fO2, 
    fCO2, 
    fN2,
    )
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
Q_gen_spec = Q_gen_net / geometry.geometry.volume

# compute skin and lung temperature
T_skin, T_lung = Tsurf_and_Tlung(;
    body = geometry, 
    k_flesh, 
    Q_gen_spec, 
    T_core)

# test lung temperature
@test T_lung ≈ u"K"((ecto_output.TLUNG)u"°C") rtol=1e-7

# solar radiation
diffuse_fraction = ecto_input.PDIF
solar_out = solar(;
    α_body_dorsal,
    α_body_ventral,
    A_silhouette,
    A_dorsal, 
    A_ventral, 
    F_ground, 
    F_sky, 
    α_ground, 
    shade, 
    zenith_angle, 
    global_radiation, 
    diffuse_fraction
    )
Q_solar = solar_out.Q_solar

# longwave radiation
ir_gain = radin(;
    A_dorsal,
    A_ventral, 
    F_sky, 
    F_ground, 
    ϵ_body_dorsal, 
    ϵ_body_ventral, 
    ϵ_ground, 
    ϵ_sky, 
    T_sky, 
    T_ground,
    )
Q_ir_in = ir_gain.Q_ir_in
ir_loss = radout(;
    T_dorsal = T_skin,
    T_ventral = T_skin,
    A_dorsal,
    A_ventral,
    F_sky,
    F_ground,
    ϵ_body_dorsal,
    ϵ_body_ventral,
    )
Q_ir_out = ir_loss.Q_ir_out

# conduction
Q_cond = conduction(;
    A_conduction, 
    L = Le, 
    T_surface = T_skin, 
    T_substrate, 
    k_substrate,
    )

# convection
conv_out = convection(; 
            body=geometry, 
            area=A_convection, 
            T_air, 
            T_surface=T_skin,
            wind_speed,
            P_atmos,
            fluid,
            fO2,
            fCO2,
            fN2,
            convection_enhancement,
            )
Q_conv = conv_out.Q_conv

# evaporation
evap_out = evaporation(;
            T_surface = T_skin,
            ψ_org,
            wetness = skin_wetness,
            area = A_convection,
            hd = conv_out.hd,
            hd_free = conv_out.hd_free,
            eye_fraction,
            T_air,
            rh,
            P_atmos,
            fO2,
            fCO2,
            fN2,
            )
Q_evap = evap_out.Q_evap

# energy balance test
@test Q_solar ≈ (ecto_output.QSOL)u"W" rtol=1e-9
@test Q_cond ≈ (ecto_output.QCOND)u"W" rtol=1e-5
@test Q_conv ≈ (ecto_output.QCONV)u"W" rtol=1e-4
@test Q_evap ≈ (ecto_output.QEVAP)u"W" rtol=1e-4

Q_in = Q_solar + Q_ir_in + Q_metab
Q_out = Q_ir_out + Q_conv + Q_evap + Q_resp + Q_cond
Q_in - Q_out

# This is broken
# (For it to work we either need an anonymous function 
# or a functor struct - you can make a struct into a function!
# then we have all the parameters attached to it.
#@test_broken T_core_s = find_zero(ectotherm, (T_air - 40K, T_air + 100K), Bisection())
#@test_broken T_core_C = (Unitful.ustrip(T_core_s) - 273.15)°C

# using structs to pass parameters

conduction_pars_external = ExternalConductionParameters(;
    conduction_fraction,
)

conduction_pars_internal = InternalConductionParameters(;
    k_flesh,
)

radiation_pars = RadiationParameters(;
    α_body_dorsal,
    α_body_ventral,
    ϵ_body_dorsal,
    ϵ_body_ventral,
    A_silhouette,
    A_dorsal, 
    A_ventral, 
    F_sky,
    F_ground,
    ventral_fraction,
)

convection_pars = ConvectionParameters(;
    convection_area = A_convection,
)

evaporation_pars = EvaporationParameters(;
    skin_wetness,
    eye_fraction,
)

hydraulic_pars = HydraulicParameters(;
    water_potential = ψ_org,
)

respiration_pars = RespirationParameters(;
    fO2_extract,
    pant,
    rq,
)

metabolism_pars = MetabolismParameters(;
    Q_metabolism = Q_metab,
    model = AndrewsPough2(),
)

traits = Traits(
    InsulationParameters(),
    conduction_pars_external,
    conduction_pars_internal,
    radiation_pars,
    convection_pars,
    evaporation_pars,
    hydraulic_pars,
    respiration_pars,
    metabolism_pars
)
# construct the Model which holds the parameters of the organism in the Organism concrete struct, of type AbstractOrganism
#lizard = Model(Organism(geometry, traits))
lizard = Organism(geometry, traits)

# get the environmental parameters
environment_pars = EnvironmentalPars(;
    α_ground,
    ϵ_ground,
    ϵ_sky,
    elevation,
    fluid,
    fN2,
    fO2,
    fCO2,
)

environment_vars = EnvironmentalVars(;
    T_air,
    T_sky,
    T_ground,
    T_substrate,
    rh,
    wind_speed,
    P_atmos,
    k_substrate,
    zenith_angle,
    global_radiation,
    diffuse_fraction,
    shade,
)

environment = (; environment_pars, environment_vars)
# define the method 'ectotherm' for passing to find_zero, which dispatches off 'lizard' 
T_air = environment_vars.T_air
ectotherm(T_air, lizard, environment)
T_core_s = find_zero(t -> ectotherm(t, lizard, environment), (T_air - 40u"K", T_air + 100u"K"), Bisection())
T_core_C = (Unitful.ustrip(T_core_s) - 273.15)u"°C"

# test core temperature calculation
@test T_core_C ≈ (ecto_output.TC)u"°C" rtol=1e-4

heat_balance_out = ectotherm(T_core_s, lizard, environment)

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