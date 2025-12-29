using HeatExchange
using Unitful, UnitfulMoles
using Test
using DataFrames, CSV

testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

ellipsoid_input_names = Symbol.(DataFrame(CSV.File("$testdir/data/ellipsoid_input.csv"))[:, 1])
ellipsoid_input_vals = DataFrame(CSV.File("$testdir/data/ellipsoid_input.csv"))[:, 2]
ellipsoid_input = (; zip(ellipsoid_input_names, ellipsoid_input_vals)...)
ellipsoid_output = DataFrame(CSV.File("$testdir/data/ellipsoid_output.csv"))[41, 2:end]
air_temperature = DataFrame(CSV.File("$testdir/data/ellipsoid_air_temperature.csv"))[41, 2]

ellipsoid_out = ellipsoid_endotherm(; 
    posture=ellipsoid_input.posture, 
    mass=(ellipsoid_input.mass)u"kg", 
    density=(ellipsoid_input.density)u"kg/m^3", 
    core_temperature=(ellipsoid_input.coreT)u"°C",
    insulation_depth=(ellipsoid_input.furdepth)u"mm", 
    insulation_conductivity=(ellipsoid_input.furcond)u"W/m/K",
    emissivity=0.95,
    oxygen_extraction_efficiency=ellipsoid_input.O2eff, 
    stress_factor=ellipsoid_input.stress,
    air_temperature=(air_temperature)u"°C", 
    wind_speed=(ellipsoid_input.windspd)u"m/s", 
    relative_humidity=ellipsoid_input.rh/100, 
    P_atmos = 101325.0u"Pa", 
    q10=ellipsoid_input.Q10,
    minimum_metabolic_rate=missing, 
    metabolic_multiplier=ellipsoid_input.basmult, 
    lethal_desiccation=0.15, 
    f_O2=0.2094)

rtol=1e-3
@testset "ellipsoid model" begin
    @test ellipsoid_out.Q_gen_required ≈ (ellipsoid_output.Qgen)u"W" rtol=rtol
    @test ellipsoid_out.Q_gen_final ≈ (ellipsoid_output.QgenFinal)u"W" rtol=rtol
    @test ellipsoid_out.O2_consumption_rate ≈ (ellipsoid_output.mlO2ph)u"ml/hr" rtol=rtol
    @test ellipsoid_out.basal_metabolic_rate_fraction ≈ ellipsoid_output.PctBasal/100 rtol=rtol
    @test ellipsoid_out.skin_temperature ≈ (ellipsoid_output.Tskin)u"°C" rtol=rtol
    @test ellipsoid_out.lower_critical_air_temperature ≈ (ellipsoid_output.LCT)u"°C" rtol=rtol
    @test ellipsoid_out.upper_critical_air_temperature ≈ (ellipsoid_output.UCT)u"°C" rtol=rtol
    @test ellipsoid_out.Q_respiration ≈ (ellipsoid_output.Qresp_W)u"W" rtol=rtol
    @test ellipsoid_out.Q_evap ≈ (ellipsoid_output.H2Oloss_W)u"W" rtol=rtol
    @test ellipsoid_out.respiratory_water_loss_rate ≈ (ellipsoid_output.Qresp_gph)u"g/hr" rtol=rtol
    @test ellipsoid_out.total_water_loss_rate ≈ (ellipsoid_output.H2O_gph)u"g/hr" rtol=rtol
    @test ellipsoid_out.fractional_mass_loss ≈ (ellipsoid_output.massph_percent/100)u"1/hr" rtol=rtol
end