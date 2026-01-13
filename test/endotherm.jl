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
            T_air=u"K"((endo_input.TA)u"°C"),
            T_air_reference=u"K"((endo_input.TAREF)u"°C"),
            T_sky=u"K"((endo_input.TSKY)u"°C"),
            T_ground=u"K"((endo_input.TGRD)u"°C"),
            T_substrate=u"K"((endo_input.TCONDSB)u"°C"),
            T_bush=u"K"((endo_input.TBUSH)u"°C"),
            T_vegetation=u"K"((endo_input.TAREF)u"°C"),
            rh=(endo_input.RH / 100.0),
            wind_speed=(endo_input.VEL)u"m/s",
            P_atmos=(endo_input.BP)u"Pa",
            zenith_angle=(endo_input.Z)u"°",
            k_substrate=(endo_input.KSUB)u"W/m/K",
            global_radiation=(endo_input.QSOLR)u"W/m^2",
            diffuse_fraction=endo_input.PDIF,
            shade=(endo_input.SHADE / 100),
        )

        environment_pars = EnvironmentalPars(;
            α_ground=endo_input.ABSSB,
            ϵ_ground=1.0,
            ϵ_sky=1.0,
            elevation=(endo_input.ELEV)u"m",
            fluid=endo_input.FLTYPE,
            gasfrac=GasFractions(
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
            k_flesh=(endo_input.AK1)u"W/m/K",
            k_fat=(endo_input.AK2)u"W/m/K",
            ρ_fat=(endo_input.FATDEN)u"kg/m^3",
        )

        if endo_input.ORIENT == 0.0
            solar_orientation = Intermediate()
        elseif endo_input.ORIENT == 1.0
            solar_orientation = NormalToSun()
        else
            solar_orientation = ParallelToSun()
        end
        radiation_pars = RadiationParameters(;
            α_body_dorsal=(1 - endo_input.REFLD),
            α_body_ventral=(1 - endo_input.REFLV),
            ϵ_body_dorsal=endo_input.EMISAN,
            ϵ_body_ventral=endo_input.EMISAN,
            F_sky=endo_input.FSKREF,
            F_ground=endo_input.FGDREF,
            F_bush=endo_input.FABUSH,
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
            fO2_extract=(endo_input.EXTREF / 100.0),
            pant=endo_input.PANT,
            rq=endo_input.RQ,
            Δ_breath=(endo_input.DELTAR)u"K",
            rh_exit=(endo_input.RELXIT / 100.0),
        )

        metabolism_pars = MetabolismParameters(;
            T_core=u"K"((endo_input.TC)u"°C"),
            Q_metabolism=(endo_input.QBASAL)u"W",
            q10=endo_input.Q10,
            model=Kleiber(),
        )

        if endo_input.FURTHRMK == 0.0
            insulation_conductivity = nothing
        else
            insulation_conductivity = (endo_input.FURTHRMK)u"W/m/K"
        end
        insulation_pars = InsulationParameters(;
            insulation_conductivity_dorsal=insulation_conductivity,
            insulation_conductivity_ventral=insulation_conductivity,
            fibre_diameter_dorsal=(endo_input.DHAIRD)u"m",
            fibre_diameter_ventral=(endo_input.DHAIRV)u"m",
            fibre_length_dorsal=(endo_input.LHAIRD)u"m",
            fibre_length_ventral=(endo_input.LHAIRV)u"m",
            insulation_depth_dorsal=(endo_input.ZFURD)u"m",
            insulation_depth_ventral=(endo_input.ZFURV)u"m",
            fibre_density_dorsal=(endo_input.RHOD)u"1/m^2",
            fibre_density_ventral=(endo_input.RHOV)u"1/m^2",
            insulation_reflectance_dorsal=endo_input.REFLD,
            insulation_reflectance_ventral=endo_input.REFLV,
            insulation_depth_compressed=(endo_input.ZFURCOMP)u"m",
            fibre_conductivity=(endo_input.KHAIR)u"W/m/K",
            longwave_depth_fraction=endo_input.XR,
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
        )

        mammal = Organism(geometry, traits)

        environment = (; environment_pars, environment_vars)

        options = EndothermMetabolicRateOptions(;
            respire=Bool(endo_input.RESPIRE),
            simulsol_tolerance=(endo_input.DIFTOL)u"K",
            resp_tolerance=endo_input.BRENTOL,
        )

        # initial conditions
        T_skin = u"K"((endo_input.TS)u"°C")
        T_insulation = u"K"((endo_input.TFA)u"°C")

        endotherm_out = solve_metabolic_rate(
            T_skin, T_insulation, mammal, environment, options
        )

        thermoregulation = endotherm_out.thermoregulation
        morphology = endotherm_out.morphology
        energy_fluxes = endotherm_out.energy_fluxes
        mass_fluxes = endotherm_out.mass_fluxes

        (; insulation_test) = insulation_properties(;
            insulation=insulation_pars,
            insulation_temperature=thermoregulation.T_insulation,
            ventral_fraction=radiation_pars.ventral_fraction,
        )

        rtol = 1e-3

        @testset "endotherm thermoregulation comparisons" begin
            @test treg_output_vec.TC ≈ ustrip(u"°C", thermoregulation.T_core) rtol = rtol
            @test treg_output_vec.TLUNG ≈ ustrip(u"°C", thermoregulation.T_lung) rtol = rtol
            @test treg_output_vec.TSKIN_D ≈ ustrip(u"°C", thermoregulation.T_skin_dorsal) rtol =
                rtol
            @test treg_output_vec.TSKIN_V ≈ ustrip(u"°C", thermoregulation.T_skin_ventral) rtol =
                rtol
            @test treg_output_vec.TFA_D ≈
                ustrip(u"°C", thermoregulation.T_insulation_dorsal) rtol = rtol
            @test treg_output_vec.TFA_V ≈
                ustrip(u"°C", thermoregulation.T_insulation_ventral) rtol = rtol
            if insulation_test > 0.0u"m"
                @test treg_output_vec.K_FUR_D ≈
                    ustrip(u"W/m/K", thermoregulation.k_insulation_dorsal) rtol = rtol
                @test treg_output_vec.K_FUR_V ≈
                    ustrip(u"W/m/K", thermoregulation.k_insulation_ventral) rtol = rtol
            end
            if isnothing(insulation_conductivity)
                @test treg_output_vec.K_FUR_EFF ≈
                    ustrip(u"W/m/K", thermoregulation.k_insulation_effective) rtol = rtol
                @test treg_output_vec.K_COMPFUR ≈
                    ustrip(u"W/m/K", thermoregulation.k_insulation_compressed) rtol = rtol
            end
        end

        fat = morphology.fat < 1.0e-10u"m" ? 0.0u"m" : morphology.fat
        rtol = 1e-6
        @testset "endotherm morphology comparisons" begin
            @test morph_output_vec.AREA ≈ ustrip(u"m^2", morphology.area_total) rtol = rtol
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
            @test morph_output_vec.F_SKY ≈ morphology.F_sky rtol = rtol
            @test morph_output_vec.F_GROUND ≈ morphology.F_ground rtol = rtol
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
        @testset "endotherm energy flux comparisons" begin
            @test enbal_output_vec.QSOL ≈ ustrip(u"W", energy_fluxes.Q_solar) rtol = rtol
            @test enbal_output_vec.QIRIN ≈ ustrip(u"W", energy_fluxes.Q_longwave_in) rtol =
                rtol
            @test enbal_output_vec.QGEN ≈ ustrip(u"W", energy_fluxes.Q_gen) rtol = rtol * 10
            @test QEVAP ≈ ustrip(u"W", energy_fluxes.Q_evaporation) rtol = rtol
            @test enbal_output_vec.QIROUT ≈ ustrip(u"W", energy_fluxes.Q_longwave_out) rtol =
                rtol
            @test enbal_output_vec.QCONV ≈ ustrip(u"W", energy_fluxes.Q_convection) rtol =
                rtol
            @test enbal_output_vec.QCOND ≈ ustrip(u"W", energy_fluxes.Q_conduction) rtol =
                rtol * 100 # could be because it's a very small number
            if !isnothing(energy_fluxes.balance)
                @test enbal_output_vec.ENB ≈ ustrip(u"W", energy_fluxes.balance) atol = 1e-3
            end
            @test enbal_output_vec.NTRY ≈ energy_fluxes.ntry
            @test Bool(enbal_output_vec.SUCCESS) ≈ energy_fluxes.success
        end

        rtol = 1e-3
        @testset "endotherm mass flux comparisons" begin
            if model_pars.respire
                @test masbal_output_vec.AIR_L ≈ ustrip(u"L/hr", mass_fluxes.V_air) rtol =
                    rtol
                @test masbal_output_vec.O2_L ≈ ustrip(u"L/hr", mass_fluxes.V_O2_STP) rtol =
                    rtol
                @test masbal_output_vec.H2OResp_g ≈ ustrip(u"g/hr", mass_fluxes.m_resp) rtol =
                    rtol
                @test masbal_output_vec.H2OCut_g ≈ ustrip(u"g/hr", mass_fluxes.m_sweat) rtol =
                    rtol
                #@test masbal_output_vec.H2O_mol_in ≈ ustrip(u"mol/hr", mass_fluxes.J_H2O_in) rtol = rtol
                #@test masbal_output_vec.H2O_mol_out ≈ ustrip(u"mol/hr", mass_fluxes.J_H2O_out) rtol = rtol
                @test masbal_output_vec.O2_mol_in ≈ ustrip(u"mol/hr", mass_fluxes.J_O2_in) rtol =
                    rtol
                @test masbal_output_vec.O2_mol_out ≈ ustrip(u"mol/hr", mass_fluxes.J_O2_out) rtol =
                    rtol
                #@test masbal_output_vec.CO2_mol_in ≈ ustrip(u"mol/hr", mass_fluxes.J_CO2_in) rtol = rtol
                #@test masbal_output_vec.CO2_mol_out ≈ ustrip(u"mol/hr", mass_fluxes.J_CO2_out) rtol = rtol
                @test masbal_output_vec.N2_mol_in ≈ ustrip(u"mol/hr", mass_fluxes.J_N2_in) rtol =
                    rtol
                @test masbal_output_vec.N2_mol_out ≈ ustrip(u"mol/hr", mass_fluxes.J_N2_out) rtol =
                    rtol
                @test masbal_output_vec.AIR_mol_in ≈ ustrip(u"mol/hr", mass_fluxes.J_air_in) rtol =
                    rtol
                @test masbal_output_vec.AIR_mol_out ≈
                    ustrip(u"mol/hr", mass_fluxes.J_air_out) rtol = rtol
            end
        end
    end
end
