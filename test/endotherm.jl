using HeatExchange
using ModelParameters
using Unitful, UnitfulMoles
using FluidProperties
using Test
using DataFrames, CSV
#using Plots

# ellipsoid model
testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

# solvendo

testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

endo_input_vec = DataFrame(CSV.File("$testdir/data/endoR_input.csv"))[:, 2]
endo_input_names = Symbol.(DataFrame(CSV.File("$testdir/data/endoR_input_names.csv"))[:, 2])
treg_output_vec = first(Tables.rowtable(DataFrame(CSV.File("$testdir/data/endoR_treg.csv"))[:, 2:end]))
morph_output_vec = first(Tables.rowtable(DataFrame(CSV.File("$testdir/data/endoR_morph.csv"))[:, 2:end]))
enbal_output_vec = first(Tables.rowtable(DataFrame(CSV.File("$testdir/data/endoR_enbal.csv"))[:, 2:end]))
masbal_output_vec = first(Tables.rowtable(DataFrame(CSV.File("$testdir/data/endoR_masbal.csv"))[:, 2:end]))

endo_input = (; zip(endo_input_names, endo_input_vec)...)

# define shape
bodyshape = Ellipsoid((endo_input.AMASS)u"kg", (endo_input.ANDENS)u"kg/m^3", 
    (endo_input.SHAPE_B), (endo_input.SHAPE_C)) 

#bodyshape = Plate((endo_input.AMASS)u"kg", (endo_input.ANDENS)u"kg/m^3", 
#    (endo_input.SHAPE_B), (endo_input.SHAPE_C)) # define shape
#bodyshape = Cylinder((endo_input.AMASS)u"kg", (endo_input.ANDENS)u"kg/m^3", 
#    (endo_input.SHAPE_B)) # define shape
#bodyshape = Sphere((endo_input.AMASS)u"kg", (endo_input.ANDENS)u"kg/m^3") # define shape

environmental_vars = EnvironmentalVars(
    u"K"((endo_input.TA)u"°C"), # T_air
    u"K"((endo_input.TAREF)u"°C"), # T_air_reference
    u"K"((endo_input.TSKY)u"°C"), # T_sky
    u"K"((endo_input.TGRD)u"°C"), # T_ground
    u"K"((endo_input.TCONDSB)u"°C"), # T_substrate
    u"K"((endo_input.TBUSH)u"°C"), # T_bush
    u"K"((endo_input.TAREF)u"°C"), # T_vegetation
    endo_input.RH/100.0, # rh
    (endo_input.VEL)u"m/s", # wind_speed
    (endo_input.BP)u"Pa", # P_atmos
    (endo_input.Z)u"°", # zenith_angle
    (endo_input.KSUB)u"W/m/K", # k_substrate
    (endo_input.QSOLR)u"W/m^2", # solar_radiation
    endo_input.PDIF, # fraction_diffuse_radiation
)

environmental_pars = EnvironmentalPars(;
    α_ground = endo_input.ABSSB,
    ϵ_ground = 1.0,
    ϵ_sky = 1.0,
    elevation = (endo_input.ELEV)u"m",
    fluid = endo_input.FLTYPE,
    fN2 = endo_input.N2GAS / 100.0,
    fO2 = endo_input.O2GAS / 100.0,
    fCO2 = endo_input.CO2GAS / 100.0,
    convection_enhancement = endo_input.CONV_ENHANCE,
    )

body_pars = BodyPars(;
    mass = (endo_input.AMASS)u"kg",
    ρ_flesh = (endo_input.ANDENS)u"kg/m^3",
    ρ_fat = (endo_input.FATDEN)u"kg/m^3",
    shape_b = endo_input.SHAPE_B,
    shape_c = endo_input.SHAPE_C,
)

integument_pars = IntegumentPars(;
    α_body_dorsal = 1 - endo_input.REFLD,
    α_body_ventral = 1 - endo_input.REFLV,
    ϵ_body_dorsal = endo_input.EMISAN,
    ϵ_body_ventral = endo_input.EMISAN,
    F_sky = endo_input.FSKREF,
    F_ground = endo_input.FGDREF,
    F_vegetation = 0.0,
    F_bush = endo_input.FABUSH,
    eye_fraction = endo_input.PCTEYES / 100.0,
    bare_skin_fraction = endo_input.PCTBAREVAP / 100.0,
    ventral_fraction = endo_input.PVEN,
)

physio_pars = PhysioPars(;
    Q_minimum = (endo_input.QBASAL)u"W",
    q10 = endo_input.Q10,
    k_fat = (endo_input.AK2)u"W/m/K",
    fO2_extract = endo_input.EXTREF / 100,
    rq = endo_input.RQ,
    Δ_breath = (endo_input.DELTAR)u"K",
    rh_exit = endo_input.RELXIT / 100.0,
)

thermoreg_pars = ThermoregulationPars(;
    insulation_step = endo_input.PZFUR,
    shape_b_step = endo_input.UNCURL,
    shape_b_max = endo_input.SHAPE_B_MAX,
    T_core_max = u"K"((endo_input.TC_MAX)u"°C"),
    T_core_min = u"K"((endo_input.TC_MIN)u"°C"),
    T_core_step = (endo_input.TC_INC)u"K",
    k_flesh_step = (endo_input.AK1_INC)u"W/m/K",
    k_flesh_max = (endo_input.AK1_MAX)u"W/m/K",
    pant_step = endo_input.PANT_INC,
    pant_multiplier = endo_input.PANT_MULT,
    pant_max = endo_input.PANT_MAX,
    skin_wetness_step = endo_input.PCTWET_INC / 100.0,
    skin_wetness_max = endo_input.PCTWET_MAX / 100.0,
    )

if endo_input.ORIENT == 0.0
    solar_orientation = Intermedate()
elseif endo_input.ORIENT == 1.0
    solar_orientation = NormalToSun()
else
    solar_orientation = ParallelToSun()
end

thermoreg_vars = ThermoregulationVars(;
    insulation_depth_dorsal = (endo_input.ZFURD)u"m",
    insulation_depth_ventral = (endo_input.ZFURV)u"m",
    fat_fraction = endo_input.FATPCT / 100.0,
    conduction_fraction = endo_input.PCOND,
    k_flesh = (endo_input.AK1)u"W/m/K",
    T_core_target = u"K"((endo_input.TC)u"°C"),
    pant = endo_input.PANT,
    skin_wetness = endo_input.PCTWET / 100.0,
    insulation_wetness = endo_input.FURWET / 100.0,
    solar_orientation
    )

if endo_input.FURTHRMK == 0.0
    insulation_conductivity = nothing
else
    insulation_conductivity = (endo_input.FURTHRMK)u"W/m/K"
end
insulation_pars = InsulationPars(;
    insulation_conductivity_dorsal = insulation_conductivity,
    insulation_conductivity_ventral = insulation_conductivity,
    fibre_diameter_dorsal = (endo_input.DHAIRD)u"m",
    fibre_diameter_ventral = (endo_input.DHAIRV)u"m",
    fibre_length_dorsal = (endo_input.LHAIRD)u"m",
    fibre_length_ventral = (endo_input.LHAIRV)u"m",
    max_insulation_depth_dorsal = (endo_input.ZFURD_MAX)u"m",
    max_insulation_depth_ventral = (endo_input.ZFURV_MAX)u"m",
    fibre_density_dorsal = (endo_input.RHOD)u"1/m^2",
    fibre_density_ventral = (endo_input.RHOV)u"1/m^2",
    insulation_reflectance_dorsal = endo_input.REFLD,
    insulation_reflectance_ventral = endo_input.REFLV,
    insulation_depth_compressed = (endo_input.ZFURCOMP)u"m",
    fibre_conductivity = (endo_input.KHAIR)u"W/m/K",
    longwave_depth_fraction = endo_input.XR,
)

organism_vars = OrganismalVars(;
    T_core = u"K"((endo_input.TC)u"°C"),
    T_skin = u"K"((endo_input.TS)u"°C"),
    T_insulation = u"K"((endo_input.TFA)u"°C"),
    T_lung = u"K"((endo_input.TC)u"°C"),
    ψ_org = 0.0u"J/kg",
)

model_pars = EndoModelPars(
    thermoregulation_mode = endo_input.TREGMODE,
    thermoregulate = Bool(endo_input.THERMOREG),
    respire = Bool(endo_input.RESPIRE),
    torpor = Bool(endo_input.TORPOR),
    simulsol_tolerance = (endo_input.DIFTOL)u"K",
    resp_tolerance = endo_input.BRENTOL,
    )

endotherm_out = endotherm(; model_pars, bodyshape, body_pars, integument_pars, insulation_pars, 
    physio_pars, thermoreg_pars, thermoreg_vars, environmental_pars, organism_vars, environmental_vars)

thermoregulation = endotherm_out.thermoregulation
morphology = endotherm_out.morphology
energy_fluxes = endotherm_out.energy_fluxes
mass_fluxes = endotherm_out.mass_fluxes

(; insulation_test) = insulation_properties(; insulation=insulation_pars, 
    insulation_temperature = thermoregulation.T_insulation, 
    ventral_fraction = integument_pars.ventral_fraction, 
    insulation_depth_dorsal = thermoreg_vars.insulation_depth_dorsal, 
    insulation_depth_ventral = thermoreg_vars.insulation_depth_ventral)
@testset "endotherm thermoregulation comparisons" begin
    @test treg_output_vec.TC ≈ ustrip(u"°C", thermoregulation.T_core) rtol = 1e-3
    @test treg_output_vec.TLUNG ≈ ustrip(u"°C", thermoregulation.T_lung) rtol = 1e-3
    @test treg_output_vec.TSKIN_D ≈ ustrip(u"°C", thermoregulation.T_skin_dorsal) rtol = 1e-3
    @test treg_output_vec.TSKIN_V ≈ ustrip(u"°C", thermoregulation.T_skin_ventral) rtol = 1e-3
    @test treg_output_vec.TFA_D ≈ ustrip(u"°C", thermoregulation.T_insulation_dorsal) rtol = 1e-3
    @test treg_output_vec.TFA_V ≈ ustrip(u"°C", thermoregulation.T_insulation_ventral) rtol = 1e-3
    @test treg_output_vec.SHAPE_B ≈ thermoregulation.shape_b rtol = 1e-4
    @test treg_output_vec.PANT ≈ thermoregulation.pant rtol = 1e-4
    @test treg_output_vec.PCTWET / 100.0 ≈ thermoregulation.skin_wetness rtol = 1e-4
    @test treg_output_vec.K_FLESH ≈ ustrip(u"W/m/K", thermoregulation.k_flesh) rtol = 1e-4
    if insulation_test > 0.0u"m"
        @test treg_output_vec.K_FUR_D ≈ ustrip(u"W/m/K", thermoregulation.k_insulation_dorsal) rtol = 1e-4
        @test treg_output_vec.K_FUR_V ≈ ustrip(u"W/m/K", thermoregulation.k_insulation_ventral) rtol = 1e-4
    end
    #if isnothing(insulation_conductivity)
        @test treg_output_vec.K_FUR_EFF ≈ ustrip(u"W/m/K", thermoregulation.k_insulation_effective) rtol = 1e-4
        @test treg_output_vec.K_COMPFUR ≈ ustrip(u"W/m/K", thermoregulation.k_insulation_compressed) rtol = 1e-4
    #end
    @test treg_output_vec.Z_FUR_D ≈ ustrip(u"m", thermoregulation.insulation_depth_dorsal) rtol = 1e-4
    @test treg_output_vec.Z_FUR_V ≈ ustrip(u"m", thermoregulation.insulation_depth_ventral) rtol = 1e-4
end

@testset "endotherm morphology comparisons" begin
    @test morph_output_vec.AREA ≈ ustrip(u"m^2", morphology.area_total) rtol = 1e-7
    @test morph_output_vec.AREA_SKIN ≈ ustrip(u"m^2", morphology.area_skin) rtol = 1e-7
    @test morph_output_vec.AREA_SKIN_EVAP ≈ ustrip(u"m^2", morphology.area_evaporation) rtol = 1e-7
    @test morph_output_vec.AREA_CONV ≈ ustrip(u"m^2", morphology.area_convection) rtol = 1e-7
    @test morph_output_vec.AREA_COND ≈ ustrip(u"m^2", morphology.area_conduction) rtol = 1e-7
    #@test morph_output_vec.AREA_SIL ≈ ustrip(u"m^2", morphology.area_silhouette) rtol = 1e-7
    @test morph_output_vec.F_SKY ≈ morphology.F_sky rtol = 1e-7
    @test morph_output_vec.F_GROUND ≈ morphology.F_ground rtol = 1e-7
    @test morph_output_vec.VOLUME ≈ ustrip(u"m^3", morphology.volume) rtol = 1e-7
    @test morph_output_vec.FLESH_VOL ≈ ustrip(u"m^3", morphology.flesh_volume) rtol = 1e-7
    @test morph_output_vec.CHAR_DIM ≈ ustrip(u"m", morphology.characteristic_dimension) rtol = 1e-7
    @test morph_output_vec.MASS_FAT ≈ ustrip(u"kg", morphology.fat_mass) rtol = 1e-7
    @test morph_output_vec.LENGTH ≈ ustrip(u"m", morphology.a_semi_major * 2) rtol = 1e-7
    @test morph_output_vec.WIDTH ≈ ustrip(u"m", morphology.b_semi_minor * 2) rtol = 1e-7
    @test morph_output_vec.HEIGHT ≈ ustrip(u"m", morphology.c_semi_minor * 2) rtol = 1e-7
    @test morph_output_vec.FAT_THICK ≈ ustrip(u"m", morphology.fat) rtol = 1e-7
end

# check for near zero
QEVAP = enbal_output_vec.QEVAP < 1.0e-20 ? 0.0 : enbal_output_vec.QEVAP

@testset "endotherm energy flux comparisons" begin
    @test enbal_output_vec.QSOL ≈ ustrip(u"W", energy_fluxes.Q_solar) rtol = 1e-3
    @test enbal_output_vec.QIRIN ≈ ustrip(u"W", energy_fluxes.Q_longwave_in) rtol = 1e-3
    @test enbal_output_vec.QGEN ≈ ustrip(u"W", energy_fluxes.Q_gen) rtol = 1e-3
    @test QEVAP ≈ ustrip(u"W", energy_fluxes.Q_evaporation) rtol = 1e-2
    @test enbal_output_vec.QIROUT ≈ ustrip(u"W", energy_fluxes.Q_longwave_out) rtol = 1e-3
    @test enbal_output_vec.QCONV ≈ ustrip(u"W", energy_fluxes.Q_convection) rtol = 1e-3
    @test enbal_output_vec.QCOND ≈ ustrip(u"W", energy_fluxes.Q_conduction) rtol = 1e-3
    if !isnothing(energy_fluxes.balance)
        @test enbal_output_vec.ENB ≈ ustrip(u"W", energy_fluxes.balance) atol = 1e-3
    end
    @test enbal_output_vec.NTRY ≈ energy_fluxes.ntry
    @test Bool(enbal_output_vec.SUCCESS) ≈ energy_fluxes.success
end

@testset "endotherm mass flux comparisons" begin
    if model_pars.respire
        @test masbal_output_vec.AIR_L ≈ ustrip(u"L/hr", mass_fluxes.V_air) rtol = 1e-3
        @test masbal_output_vec.O2_L ≈ ustrip(u"L/hr", mass_fluxes.V_O2_STP) rtol = 1e-3
        @test masbal_output_vec.H2OResp_g ≈ ustrip(u"g/hr", mass_fluxes.m_resp) rtol = 1e-3
        @test masbal_output_vec.H2OCut_g ≈ ustrip(u"g/hr", mass_fluxes.m_sweat) rtol = 1e-3
        #@test masbal_output_vec.H2O_mol_in ≈ ustrip(u"mol/hr", mass_fluxes.J_H2O_in) rtol = 1e-8
        #@test masbal_output_vec.H2O_mol_out ≈ ustrip(u"mol/hr", mass_fluxes.J_H2O_out) rtol = 1e-8
        @test masbal_output_vec.O2_mol_in ≈ ustrip(u"mol/hr", mass_fluxes.J_O2_in) rtol = 1e-3
        @test masbal_output_vec.O2_mol_out ≈ ustrip(u"mol/hr", mass_fluxes.J_O2_out) rtol = 1e-3
        #@test masbal_output_vec.CO2_mol_in ≈ ustrip(u"mol/hr", mass_fluxes.J_CO2_in) rtol = 1e-8
        #@test masbal_output_vec.CO2_mol_out ≈ ustrip(u"mol/hr", mass_fluxes.J_CO2_out) rtol = 1e-8
        @test masbal_output_vec.N2_mol_in ≈ ustrip(u"mol/hr", mass_fluxes.J_N2_in) rtol = 1e-3
        @test masbal_output_vec.N2_mol_out ≈ ustrip(u"mol/hr", mass_fluxes.J_N2_out) rtol = 1e-3      
        @test masbal_output_vec.AIR_mol_in ≈ ustrip(u"mol/hr", mass_fluxes.J_air_in) rtol = 1e-3
        @test masbal_output_vec.AIR_mol_out ≈ ustrip(u"mol/hr", mass_fluxes.J_air_out) rtol = 1e-3 
    end      
end

