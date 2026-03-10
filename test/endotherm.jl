using HeatExchange
using BiophysicalGeometry
using ModelParameters
using Unitful, UnitfulMoles
using FluidProperties
using Test
using DataFrames, CSV

testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

endo_input_names = Symbol.(DataFrame(CSV.File("$testdir/data/endoR_input_names.csv"))[:, 2])
shape_number = 4
furmult = 0
# loop through all shapes with and without fur
for shape_number in 1:4
    for furmult in 0:1
        endo_input_vec = DataFrame(
            CSV.File(
                "$testdir/data/endoR_input_" *
                string(shape_number) *
                "_" *
                string(furmult) *
                ".csv",
            ),
        )[
            :, 2
        ]
        treg_output_vec = first(
            Tables.rowtable(
                DataFrame(
                    CSV.File(
                        "$testdir/data/endoR_treg_" *
                        string(shape_number) *
                        "_" *
                        string(furmult) *
                        ".csv",
                    ),
                )[
                    :, 2:end
                ],
            ),
        )
        morph_output_vec = first(
            Tables.rowtable(
                DataFrame(
                    CSV.File(
                        "$testdir/data/endoR_morph_" *
                        string(shape_number) *
                        "_" *
                        string(furmult) *
                        ".csv",
                    ),
                )[
                    :, 2:end
                ],
            ),
        )
        enbal_output_vec = first(
            Tables.rowtable(
                DataFrame(
                    CSV.File(
                        "$testdir/data/endoR_enbal_" *
                        string(shape_number) *
                        "_" *
                        string(furmult) *
                        ".csv",
                    ),
                )[
                    :, 2:end
                ],
            ),
        )
        masbal_output_vec = first(
            Tables.rowtable(
                DataFrame(
                    CSV.File(
                        "$testdir/data/endoR_masbal_" *
                        string(shape_number) *
                        "_" *
                        string(furmult) *
                        ".csv",
                    ),
                )[
                    :, 2:end
                ],
            ),
        )

        endo_input = (; zip(endo_input_names, endo_input_vec)...)

        # define shape
        if endo_input.SHAPE == 1
            shape_pars = Cylinder(
                (endo_input.AMASS)u"kg", (endo_input.ANDENS)u"kg/m^3", (endo_input.SHAPE_B)
            ) # define shape
        end
        if endo_input.SHAPE == 2
            shape_pars = Sphere((endo_input.AMASS)u"kg", (endo_input.ANDENS)u"kg/m^3")
        end
        if endo_input.SHAPE == 3
            shape_pars = Plate(
                (endo_input.AMASS)u"kg",
                (endo_input.ANDENS)u"kg/m^3",
                (endo_input.SHAPE_B),
                (endo_input.SHAPE_C),
            )
        end
        if endo_input.SHAPE == 4
            shape_pars = Ellipsoid(
                (endo_input.AMASS)u"kg",
                (endo_input.ANDENS)u"kg/m^3",
                (endo_input.SHAPE_B),
                (endo_input.SHAPE_C),
            )
        end

        fat = Fat(endo_input.FATPCT / 100.0, (endo_input.FATDEN)u"kg/m^3")
        mean_insulation_depth =
            (endo_input.ZFURD * (1 - endo_input.PVEN) + endo_input.ZFURV * endo_input.PVEN)u"m"
        mean_fibre_diameter = (
            endo_input.DHAIRD * (1 - endo_input.PVEN) + endo_input.DHAIRV * endo_input.PVEN
        )u"m"
        mean_fibre_density =
            (endo_input.RHOD * (1 - endo_input.PVEN) + endo_input.RHOV * endo_input.PVEN)u"1/m^2"
        fur = Fur(mean_insulation_depth, mean_fibre_diameter, mean_fibre_density)
        geometry = Body(shape_pars, CompositeInsulation(fur, fat))

        environment_vars = EnvironmentalVars(;
            air_temperature=u"K"((endo_input.TA)u"°C"),
            reference_air_temperature=u"K"((endo_input.TAREF)u"°C"),
            sky_temperature=u"K"((endo_input.TSKY)u"°C"),
            ground_temperature=u"K"((endo_input.TGRD)u"°C"),
            substrate_temperature=u"K"((endo_input.TCONDSB)u"°C"),
            bush_temperature=u"K"((endo_input.TBUSH)u"°C"),
            vegetation_temperature=u"K"((endo_input.TAREF)u"°C"),
            relative_humidity=(endo_input.RH / 100.0),
            wind_speed=(endo_input.VEL)u"m/s",
            atmospheric_pressure=(endo_input.BP)u"Pa",
            zenith_angle=(endo_input.Z)u"°",
            substrate_conductivity=(endo_input.KSUB)u"W/m/K",
            global_radiation=(endo_input.QSOLR)u"W/m^2",
            diffuse_fraction=endo_input.PDIF,
            shade=(endo_input.SHADE / 100),
        )

        environment_pars = EnvironmentalPars(;
            ground_albedo=endo_input.ABSSB,
            ground_emissivity=1.0,
            sky_emissivity=1.0,
            elevation=(endo_input.ELEV)u"m",
            fluid=endo_input.FLTYPE,
            gas_fractions=GasFractions(
                endo_input.O2GAS / 100.0,
                endo_input.CO2GAS / 100.0,
                endo_input.N2GAS / 100.0,
            ),
            convection_enhancement=endo_input.CONV_ENHANCE,
        )

        conduction_pars_external = ExternalConductionParameters(;
            conduction_fraction=endo_input.PCOND
        )

        conduction_pars_internal = InternalConductionParameters(;
            fat_fraction=(endo_input.FATPCT / 100.0),
            flesh_conductivity=(endo_input.AK1)u"W/m/K",
            fat_conductivity=(endo_input.AK2)u"W/m/K",
            fat_density=(endo_input.FATDEN)u"kg/m^3",
        )

        if endo_input.ORIENT == 0.0
            solar_orientation = Intermediate()
        elseif endo_input.ORIENT == 1.0
            solar_orientation = NormalToSun()
        else
            solar_orientation = ParallelToSun()
        end
        radiation_pars = RadiationParameters(;
            body_absorptivity_dorsal=(1 - endo_input.REFLD),
            body_absorptivity_ventral=(1 - endo_input.REFLV),
            body_emissivity_dorsal=endo_input.EMISAN,
            body_emissivity_ventral=endo_input.EMISAN,
            sky_view_factor=endo_input.FSKREF,
            ground_view_factor=endo_input.FGDREF,
            bush_view_factor=endo_input.FABUSH,
            ventral_fraction=endo_input.PVEN,
            solar_orientation,
        )

        evaporation_pars = EvaporationParameters(;
            skin_wetness=(endo_input.PCTWET / 100.0),
            insulation_wetness=(endo_input.FURWET / 100.0),
            eye_fraction=(endo_input.PCTEYES / 100.0),
            bare_skin_fraction=(endo_input.PCTBAREVAP / 100.0),
            insulation_fraction=1.0,
        )

        hydraulic_pars = HydraulicParameters(;)

        respiration_pars = RespirationParameters(;
            oxygen_extraction_efficiency=(endo_input.EXTREF / 100.0),
            pant=endo_input.PANT,
            respiratory_quotient=endo_input.RQ,
            exhaled_temperature_offset=(endo_input.DELTAR)u"K",
            exhaled_relative_humidity=(endo_input.RELXIT / 100.0),
        )

        metabolism_pars = MetabolismParameters(;
            core_temperature=u"K"((endo_input.TC)u"°C"),
            metabolic_heat_flow=(endo_input.QBASAL)u"W",
            q10=endo_input.Q10,
            model=Kleiber(),
        )

        insulation_pars = InsulationParameters(;
            dorsal=FibreProperties(;
                diameter=(endo_input.DHAIRD)u"m",
                length=(endo_input.LHAIRD)u"m",
                density=(endo_input.RHOD)u"1/m^2",
                depth=(endo_input.ZFURD)u"m",
                reflectance=endo_input.REFLD,
                conductivity=(endo_input.KHAIR)u"W/m/K",
            ),
            ventral=FibreProperties(;
                diameter=(endo_input.DHAIRV)u"m",
                length=(endo_input.LHAIRV)u"m",
                density=(endo_input.RHOV)u"1/m^2",
                depth=(endo_input.ZFURV)u"m",
                reflectance=endo_input.REFLV,
                conductivity=(endo_input.KHAIR)u"W/m/K",
            ),
            depth_compressed=(endo_input.ZFURCOMP)u"m",
            longwave_depth_fraction=endo_input.XR,
        )

        options = SolveMetabolicRateOptions(;
            respire=Bool(endo_input.RESPIRE),
            temperature_tolerance=(endo_input.DIFTOL)u"K",
            resp_tolerance=endo_input.BRENTOL,
        )

        traits = HeatExchangeTraits(
            shape_pars,
            insulation_pars,
            conduction_pars_external,
            conduction_pars_internal,
            radiation_pars,
            ConvectionParameters(),
            evaporation_pars,
            hydraulic_pars,
            respiration_pars,
            metabolism_pars,
            options,
        )

        mammal = Organism(geometry, traits)

        environment = (; environment_pars, environment_vars)

        # initial conditions
        skin_temperature = u"K"((endo_input.TS)u"°C")
        insulation_temperature = u"K"((endo_input.TFA)u"°C")

        endotherm_out = solve_metabolic_rate(
            mammal, environment, skin_temperature, insulation_temperature
        )

        thermoregulation = endotherm_out.thermoregulation
        morphology = endotherm_out.morphology
        energy_flows = endotherm_out.energy_flows
        mass_flows = endotherm_out.mass_flows

        (; insulation_test) = insulation_properties(
            insulation_pars, thermoregulation.insulation_temperature, radiation_pars.ventral_fraction
        )

        rtol = 1e-3

        @testset "endotherm thermoregulation comparisons" begin
            @test treg_output_vec.TC ≈ ustrip(u"°C", thermoregulation.core_temperature) rtol = rtol
            @test treg_output_vec.TLUNG ≈ ustrip(u"°C", thermoregulation.lung_temperature) rtol = rtol
            @test treg_output_vec.TSKIN_D ≈ ustrip(u"°C", thermoregulation.skin_temperature_dorsal) rtol =
                rtol
            @test treg_output_vec.TSKIN_V ≈ ustrip(u"°C", thermoregulation.skin_temperature_ventral) rtol =
                rtol
            @test treg_output_vec.TFA_D ≈
                ustrip(u"°C", thermoregulation.insulation_temperature_dorsal) rtol = rtol
            @test treg_output_vec.TFA_V ≈
                ustrip(u"°C", thermoregulation.insulation_temperature_ventral) rtol = rtol
            if insulation_test > 0.0u"m"
                @test treg_output_vec.K_FUR_D ≈
                    ustrip(u"W/m/K", thermoregulation.insulation_conductivity_dorsal) rtol = rtol
                @test treg_output_vec.K_FUR_V ≈
                    ustrip(u"W/m/K", thermoregulation.insulation_conductivity_ventral) rtol = rtol
            end
            @test treg_output_vec.K_FUR_EFF ≈
                ustrip(u"W/m/K", thermoregulation.insulation_conductivity_effective) rtol = rtol
            @test treg_output_vec.K_COMPFUR ≈
                ustrip(u"W/m/K", thermoregulation.insulation_conductivity_compressed) rtol = rtol
        end

        fat = morphology.fat < 1.0e-10u"m" ? 0.0u"m" : morphology.fat
        rtol = 1e-6
        @testset "endotherm morphology comparisons" begin
            @test morph_output_vec.AREA ≈ ustrip(u"m^2", morphology.total_area) rtol = rtol
            @test morph_output_vec.AREA_SKIN ≈ ustrip(u"m^2", morphology.area_skin) rtol =
                rtol
            @test morph_output_vec.AREA_SKIN_EVAP ≈
                ustrip(u"m^2", morphology.area_evaporation) rtol = rtol
            @test morph_output_vec.AREA_CONV ≈ ustrip(u"m^2", morphology.area_convection) rtol =
                rtol
            @test morph_output_vec.AREA_COND ≈ ustrip(u"m^2", morphology.area_conduction) rtol =
                rtol
            @test morph_output_vec.AREA_SIL ≈ ustrip(u"m^2", morphology.area_silhouette) rtol =
                rtol
            @test morph_output_vec.F_SKY ≈ morphology.sky_view_factor rtol = rtol
            @test morph_output_vec.F_GROUND ≈ morphology.ground_view_factor rtol = rtol
            @test morph_output_vec.VOLUME ≈ ustrip(u"m^3", morphology.volume) rtol = rtol
            @test morph_output_vec.FLESH_VOL ≈ ustrip(u"m^3", morphology.volume_flesh) rtol =
                rtol
            @test morph_output_vec.CHAR_DIM ≈
                ustrip(u"m", morphology.characteristic_dimension) rtol = rtol
            @test morph_output_vec.MASS_FAT ≈ ustrip(u"kg", morphology.fat_mass) rtol = rtol
            if mammal.body.shape isa Cylinder
                @test morph_output_vec.LENGTH ≈ ustrip(u"m", morphology.length_fur) rtol =
                    rtol
                @test morph_output_vec.WIDTH ≈ ustrip(u"m", morphology.radius_fur * 2) rtol =
                    rtol
            end
            if mammal.body.shape isa Sphere
                @test morph_output_vec.LENGTH ≈ ustrip(u"m", morphology.radius_fur * 2) rtol =
                    rtol
                @test morph_output_vec.WIDTH ≈ ustrip(u"m", morphology.radius_fur * 2) rtol =
                    rtol
            end
            if mammal.body.shape isa Plate
                @test morph_output_vec.LENGTH ≈ ustrip(u"m", morphology.length_fur) rtol =
                    rtol
                @test morph_output_vec.WIDTH ≈ ustrip(u"m", morphology.width_fur) rtol =
                    rtol
            end
            if mammal.body.shape isa Ellipsoid
                @test morph_output_vec.LENGTH ≈
                    ustrip(u"m", morphology.a_semi_major_fur * 2) rtol = rtol
                @test morph_output_vec.WIDTH ≈ ustrip(u"m", morphology.b_semi_minor_fur * 2) rtol =
                    rtol
                @test morph_output_vec.HEIGHT ≈
                    ustrip(u"m", morphology.c_semi_minor_fur * 2) rtol = rtol
            end
            @test morph_output_vec.FAT_THICK ≈ ustrip(u"m", fat) rtol = rtol
        end

        # check for near zero
        QEVAP = enbal_output_vec.QEVAP < 1.0e-20 ? 0.0 : enbal_output_vec.QEVAP

        rtol = 1e-3
        @testset "endotherm energy flow comparisons" begin
            @test enbal_output_vec.QSOL ≈ ustrip(u"W", energy_flows.solar_flow) rtol = rtol
            @test enbal_output_vec.QIRIN ≈ ustrip(u"W", energy_flows.longwave_flow_in) rtol =
                rtol
            @test enbal_output_vec.QGEN ≈ ustrip(u"W", energy_flows.generated_heat_flow) rtol = rtol * 10
            @test QEVAP ≈ ustrip(u"W", energy_flows.evaporation_heat_flow) rtol = rtol
            @test enbal_output_vec.QIROUT ≈ ustrip(u"W", energy_flows.longwave_flow_out) rtol =
                rtol
            @test enbal_output_vec.QCONV ≈ ustrip(u"W", energy_flows.convection_heat_flow) rtol =
                rtol
            @test enbal_output_vec.QCOND ≈ ustrip(u"W", energy_flows.conduction_flow) rtol =
                rtol * 100 # could be because it's a very small number
            if !isnothing(energy_flows.balance)
                @test enbal_output_vec.ENB ≈ ustrip(u"W", energy_flows.balance) atol = 1e-3
            end
            @test enbal_output_vec.NTRY ≈ energy_flows.ntry
            @test Bool(enbal_output_vec.SUCCESS) ≈ energy_flows.success
        end

        rtol = 1e-3
        @testset "endotherm mass flow comparisons" begin
            if options.respire
                @test masbal_output_vec.AIR_L ≈ ustrip(u"L/hr", mass_flows.air_flow) rtol =
                    rtol
                @test masbal_output_vec.O2_L ≈ ustrip(u"L/hr", mass_flows.oxygen_flow_standard) rtol =
                    rtol
                @test masbal_output_vec.H2OResp_g ≈ ustrip(u"g/hr", mass_flows.respiration_mass) rtol =
                    rtol
                @test masbal_output_vec.H2OCut_g ≈ ustrip(u"g/hr", mass_flows.m_sweat) rtol =
                    rtol
                #@test masbal_output_vec.H2O_mol_in ≈ ustrip(u"mol/hr", mass_flows.molar_fluxes_in.water) rtol = rtol
                #@test masbal_output_vec.H2O_mol_out ≈ ustrip(u"mol/hr", mass_flows.molar_fluxes_out.water) rtol = rtol
                @test masbal_output_vec.O2_mol_in ≈ ustrip(u"mol/hr", mass_flows.molar_fluxes_in.oxygen) rtol =
                    rtol
                @test masbal_output_vec.O2_mol_out ≈ ustrip(u"mol/hr", mass_flows.molar_fluxes_out.oxygen) rtol =
                    rtol
                #@test masbal_output_vec.CO2_mol_in ≈ ustrip(u"mol/hr", mass_flows.molar_fluxes_in.carbon_dioxide) rtol = rtol
                #@test masbal_output_vec.CO2_mol_out ≈ ustrip(u"mol/hr", mass_flows.molar_fluxes_out.carbon_dioxide) rtol = rtol
                @test masbal_output_vec.N2_mol_in ≈ ustrip(u"mol/hr", mass_flows.molar_fluxes_in.nitrogen) rtol =
                    rtol
                @test masbal_output_vec.N2_mol_out ≈ ustrip(u"mol/hr", mass_flows.molar_fluxes_out.nitrogen) rtol =
                    rtol
                @test masbal_output_vec.AIR_mol_in ≈ ustrip(u"mol/hr", mass_flows.molar_fluxes_in.air) rtol =
                    rtol
                @test masbal_output_vec.AIR_mol_out ≈
                    ustrip(u"mol/hr", mass_flows.molar_fluxes_out.air) rtol = rtol
            end
        end
    end
end
