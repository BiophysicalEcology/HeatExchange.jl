using HeatExchange
using Unitful
using Test
using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R
using CSV, DataFrames
using ModelParameters
using Roots

testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

names_in = [
    :Ww_g, :alpha, :epsilon, :rho_body, :fatosk, :fatosb, :shape, :shape_a, :shape_b, :shape_c, 
    :conv_enhance, :pct_cond, :pct_touch, :postur, :orient, :k_flesh, :M_1, :M_2, 
    :M_3, :M_4, :Q_act, :pct_wet, :pct_eyes, :pct_mouth, :psi_body, :pantmax, :F_O2, :RQ, :delta_air, 
    :leaf, :g_vs_ab, :g_vs_ad, :elev, :alpha_sub, :epsilon_sub, :epsilon_sky, :pres, :fluid, :O2gas, 
    :CO2gas, :N2gas, :K_sub, :PDIF, :SHADE, :QSOLR, :Z, :TA, :TGRD, :TSUBST, :TSKY, :VEL, :RH
]
input_NMR_vec = CSV.read("$testdir/data/ectoR_input.csv", DataFrame)[:, 2]
input_NMR = (; zip(names_in, input_NMR_vec)...)

names_output_TC = [
    :TC, :TSKIN, :TLUNG
]
output_TC_vec = CSV.read("$testdir/data/temperatures.csv", DataFrame)[:, 2]
output_TC = (; zip(names_output_TC, output_TC_vec)...)

output_enbal_vec = CSV.read("$testdir/data/enbal.csv", DataFrame)
output_enbal = NamedTuple(output_enbal_vec[1, 2:end])
output_masbal_vec = CSV.read("$testdir/data/masbal.csv", DataFrame)
output_masbal = NamedTuple(output_masbal_vec[1, 2:end])
output_geom_vec = CSV.read("$testdir/data/geom.csv", DataFrame)
output_geom = (; zip(Symbol.(collect(output_geom_vec[:, 1])), collect(output_geom_vec[:, 2]))...)
output_solar_vec = CSV.read("$testdir/data/solar.csv", DataFrame)
output_solar = (; zip(Symbol.(collect(output_solar_vec[:, 1])), collect(output_solar_vec[:, 2]))...)
output_resp_vec = CSV.read("$testdir/data/resp.csv", DataFrame)
output_resp = (; zip(Symbol.(collect(output_resp_vec[:, 1])), collect(output_resp_vec[:, 2]))...)
output_conv_vec = CSV.read("$testdir/data/conv.csv", DataFrame)
output_conv = (; zip(Symbol.(collect(output_conv_vec[:, 1])), collect(output_conv_vec[:, 2]))...)
output_evap_vec = CSV.read("$testdir/data/evap.csv", DataFrame)
output_evap = (; zip(Symbol.(collect(output_evap_vec[:, 1])), collect(output_evap_vec[:, 2]))...)

# set the environmental variables
environmental_vars = EnvironmentalVars(
    T_air = K((input_NMR.TA)°C),
    T_sky = K((input_NMR.TSKY)°C),
    T_sub = K((input_NMR.TSUBST)°C),
    rh = input_NMR.RH,
    vel = (input_NMR.VEL)m/s,
    P_atmos = (input_NMR.pres)Pa,
    zen = (input_NMR.Z)°,
    k_sub = (input_NMR.K_sub)W/m/K,
    Q_sol = (input_NMR.QSOLR)W/m^2,
    Q_dir = (input_NMR.QSOLR)W/m^2 / cos((input_NMR.Z)°) - (input_NMR.QSOLR)W/m^2 * input_NMR.PDIF,
    Q_dif = (input_NMR.QSOLR)W/m^2 * input_NMR.PDIF,
)

# set the environmental parameters
environmental_params = EnvironmentalPars(
    α_sub = Param(input_NMR.alpha_sub, bounds=(0.0, 1.0)),
    ϵ_sub = Param(input_NMR.epsilon_sub, bounds=(0.2, 1.0)),
    ϵ_sky = Param(input_NMR.epsilon_sky, bounds=(0.2, 1.0)),
    elev = Param((input_NMR.elev), units=u"m"),
    fluid = Param(input_NMR.fluid),
    fN2 = Param(input_NMR.N2gas / 100),
    fO2 = Param(input_NMR.O2gas / 100),
    fCO2 = Param(input_NMR.CO2gas / 100),
)

# organism morphology
mass = kg((input_NMR.Ww_g)g)
ρ_body = (input_NMR.rho_body)kg/m^3
shapeb = input_NMR.shape_b
shapec = input_NMR.shape_c
shape_body = Ellipsoid(mass, ρ_body, shapeb, shapec)
 # construct a Body, which is naked - this constructor will apply the
 # 'geometry' function to the inputs and return a struct that has the
 # struct for the 'Shape' type, as well as the insulation and the geometry struct
geometric_traits = Body(shape_body, Naked())


morphology = MorphoPars(
    α_org_dorsal = Param(input_NMR.alpha, bounds=(0.2, 1.0)),
    α_org_ventral = Param(input_NMR.alpha, bounds=(0.2, 1.0)),
    ϵ_org_dorsal = Param(input_NMR.epsilon, bounds=(0.1, 1.0)),
    ϵ_org_ventral = Param(input_NMR.epsilon, bounds=(0.1, 1.0)),
    F_sky = Param(input_NMR.fatosk, bounds=(0.3, 0.5)),
    F_sub = Param(input_NMR.fatosb, bounds=(0.3, 0.5)),
    k_body = Param(input_NMR.k_flesh, bounds=(0.412, 2.8), units=u"W/m/K"),
    p_eyes = Param((input_NMR.pct_eyes) / 100, bounds=(0.0, 4e-4)),
    p_wet = Param(input_NMR.pct_wet/100, bounds=(0.0, 1.0)),
    p_cond = Param((input_NMR.pct_cond) / 100, bounds=(0.0, 1.0)),
)

physiology = PhysioPars(
    fO2_extract = Param(input_NMR.F_O2 / 100, bounds=(0.10, 0.30)),
    rq = Param(input_NMR.RQ, bounds=(0.7, 0.9)),
    M1 = Param(input_NMR.M_1, bounds=(0.01, 0.02)),
    M2 = Param(input_NMR.M_2, bounds=(0.7, 0.9)),
    M3 = Param(input_NMR.M_3, bounds=(0.02, 0.04)),
    pant = Param(input_NMR.pantmax, bounds=(1.0, 10.0)),
)

organism_vars = OrganismalVars(
    T_core = K((input_NMR.TA)°C),
    T_surf = K((input_NMR.TA)°C),
    T_lung = K((input_NMR.TA)°C),
    ψ_org = (input_NMR.psi_body)J/kg,
)

metab_out = metabolism(mass, u"K"((output_TC.TC)u"°C"), physiology.M1, physiology.M2, physiology.M3)
metab_out.Q_metab
output_enbal.QMET

resp_out = respiration(u"K"((output_TC.TC)u"°C"), metab_out.Q_metab, physiology.fO2_extract.val,
 physiology.pant.val, physiology.rq.val, environmental_vars.T_air, environmental_vars.rh,
  environmental_params.elev, environmental_vars.P_atmos, environmental_params.fO2.val,
   environmental_params.fCO2.val, environmental_params.fN2.val)
resp_out.Q_resp
output_resp.QRESP
A_tot = geometric_traits.geometry.area
A_cond = A_tot * morphology.p_cond.val
A_conv = A_tot * (1 - morphology.p_cond.val)
A_sil = calc_silhouette_area(geometric_traits, Z)
A_sil = calc_silhouette_area(geometric_traits, environmental_vars.zen)
solar_out = solar(morphology.α_org_dorsal.val, morphology.α_org_ventral.val, A_sil,
 geometric_traits.geometry.area, A_cond, morphology.F_sub.val, morphology.F_sky.val, 
 environmental_params.α_sub.val, environmental_vars.Q_sol, environmental_vars.Q_dir, 
 environmental_vars.Q_dif)
solar_out.Q_solar
output_solar.QSOLAR

# construct the Model which holds the parameters of the organism in the Organism concrete struct, of type AbstractOrganism
lizard = Model(Organism(geometric_traits, morphology, physiology))

# get the variables for both the organism and environment
variables = (organism=organism_vars, environment=environmental_vars)

# define the method 'heat_balance' for passing to find_zero, which dispatches off 'lizard' 
T_air = EnvironmentalVars().T_air
balance = heat_balance(T_air, lizard, environmental_params, variables)
T_core_s = find_zero(t -> heat_balance(t, lizard, environmental_params, variables), (T_air - 40K, T_air + 100K), Bisection())
T_core_C = °C(T_core_s)
heat_balance_out = heat_balance(T_core_s, lizard, environmental_params, variables)

heat_balance_out.enbal[1]
output_enbal[1]

u"°C"(heat_balance_out.T_core)
(output_TC.TC)u"°C"

@testset "heat budget" begin
    @test heat_balance_out.T_core ≈ u"K"((output_TC.TC)u"°C") atol=1e-6u"K"
    @test heat_balance_out.enbal ≈ (collect(output_enbal))u"W" atol=1e-6u"W"
    @test frog_out.geometry.volume ≈ (frog_out_NMR.VOL)u"m^3" atol=1e-6u"m^3"
    @test frog_silhouette_parallel ≈ (frog_out_NMR.ASILP)u"m^2" atol=1e-6u"m^2"
    @test frog_silhouette_normal ≈ (frog_out_NMR.ASILN)u"m^2" atol=1e-6u"m^2"
end