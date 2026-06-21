function example_environment_vars(;
    air_temperature=u"K"((20.0)u"°C"),
    relative_humidity=0.05,
    wind_speed=0.1u"m/s",
    atmospheric_pressure=101325.0u"Pa",
    zenith_angle=20.0u"°",
    substrate_conductivity=2.79u"W/m/K",
    global_radiation=0.0u"W/m^2",
    diffuse_fraction=0.1,
    shade=0,
)
    EnvironmentalVars(;
        air_temperature,
        reference_air_temperature = air_temperature,
        sky_temperature           = air_temperature,
        ground_temperature        = air_temperature,
        substrate_temperature     = air_temperature,
        bush_temperature          = air_temperature,
        vegetation_temperature    = air_temperature,
        relative_humidity,
        wind_speed,
        atmospheric_pressure,
        zenith_angle,
        substrate_conductivity,
        global_radiation,
        diffuse_fraction,
        shade,
    )
end

function example_environment_pars(;
    ground_albedo=0.8,
    ground_emissivity=1.0,
    sky_emissivity=1.0,
    elevation=0.0u"m",
    fluid=Air(),
    gas_fractions=GasFractions(),
    convection_enhancement=1.0,
)
    EnvironmentalPars(;
        ground_albedo,
        ground_emissivity,
        sky_emissivity,
        elevation,
        fluid,
        gas_fractions,
        convection_enhancement,
    )
end

function example_ellipsoid_shape_pars(;
    mass=65.0u"kg",
    ρ_flesh=1000.0u"kg/m^3",
    axis_ratio_b=1.1,
    axis_ratio_c=1.1,
)
    Ellipsoid(mass, ρ_flesh, axis_ratio_b, axis_ratio_c)
end

example_shape_pars(; kwargs...) = example_ellipsoid_shape_pars(; kwargs...)

function example_conduction_pars_external(;
    conduction_fraction=0.0,
)
    ExternalConductionParameters(;
        conduction_fraction,
    )
end

function example_conduction_pars_internal(;
    fat_fraction=0.0,
    flesh_conductivity=0.9u"W/m/K",
    fat_conductivity=0.23u"W/m/K",
    fat_density=901.0u"kg/m^3",
)
    InternalConductionParameters(;
        fat_fraction,
        flesh_conductivity,
        fat_conductivity,
        fat_density,
    )
end

function example_radiation_pars(;
    body_absorptivity_dorsal=0.8,
    body_absorptivity_ventral=0.8,
    body_emissivity_dorsal=0.99,
    body_emissivity_ventral=0.99,
    sky_view_factor=0.5,
    ground_view_factor=0.5,
    bush_view_factor=0.0,
    ventral_fraction=0.5,
    solar_orientation=Intermediate(),
)
    RadiationParameters(;
        body_absorptivity_dorsal,
        body_absorptivity_ventral,
        body_emissivity_dorsal,
        body_emissivity_ventral,
        sky_view_factor,
        ground_view_factor,
        bush_view_factor,
        ventral_fraction,
        solar_orientation,
    )
end

function example_evaporation_pars(;
    skin_wetness=0.005,
    insulation_wetness=0.0,
    eye_fraction=0.0,
    bare_skin_fraction=0.0,
    insulation_fraction=1.0,
)
    AnimalEvaporationParameters(;
        skin_wetness,
        insulation_wetness,
        eye_fraction,
        bare_skin_fraction,
        insulation_fraction,
    )
end

function example_hydraulic_pars(;
    water_potential=0.0u"J/kg",
    hydraulic_conductance=0.0u"kg / (m^2 * s * (J/kg))",
    specific_hydration=0.000304u"m^3 / (m^3 * (J/kg))",
)
    HydraulicParameters(;
        water_potential,
        hydraulic_conductance,
        specific_hydration,
    )
end

function example_respiration_pars(;
    oxygen_extraction_efficiency=0.2,
    pant=1.0,
    respiratory_quotient=0.8,
    exhaled_temperature_offset=0.0u"K",
    exhaled_relative_humidity=1.0,
)
    RespirationParameters(;
        oxygen_extraction_efficiency,
        pant,
        respiratory_quotient,
        exhaled_temperature_offset,
        exhaled_relative_humidity,
    )
end

function example_metabolism_pars(;
    core_temperature=u"K"((37.0)u"°C"),
    metabolic_heat_flow=77.61842u"W",
    q10=2.0,
    model=Kleiber(),
)
    MetabolismParameters(;
        core_temperature,
        metabolic_heat_flow,
        q10,
        model,
    )
end

function example_insulation_pars(;
    fibre_diameter_dorsal=30e-06u"m",
    fibre_diameter_ventral=30e-06u"m",
    fibre_length_dorsal=23.9e-03u"m",
    fibre_length_ventral=23.9e-03u"m",
    insulation_depth_dorsal=2e-03u"m",
    insulation_depth_ventral=2e-03u"m",
    fibre_density_dorsal=3000e+04u"1/m^2",
    fibre_density_ventral=3000e+04u"1/m^2",
    insulation_reflectance_dorsal=0.2,
    insulation_reflectance_ventral=0.2,
    insulation_depth_compressed=2e-03u"m",
    fibre_conductivity=0.209u"W/m/K",
    longwave_depth_fraction=1.0,
)
    InsulationParameters(;
        dorsal=FibreProperties(;
            diameter=fibre_diameter_dorsal,
            length=fibre_length_dorsal,
            density=fibre_density_dorsal,
            depth=insulation_depth_dorsal,
            reflectance=insulation_reflectance_dorsal,
            conductivity=fibre_conductivity,
        ),
        ventral=FibreProperties(;
            diameter=fibre_diameter_ventral,
            length=fibre_length_ventral,
            density=fibre_density_ventral,
            depth=insulation_depth_ventral,
            reflectance=insulation_reflectance_ventral,
            conductivity=fibre_conductivity,
        ),
        depth_compressed=insulation_depth_compressed,
        longwave_depth_fraction,
    )
end

function example_convection_pars(;
    convection_area=0.0u"m^2",
    characteristic_dimension_formula=VolumeCubeRoot(),
)
    ConvectionParameters(;
        convection_area,
        characteristic_dimension_formula,
    )
end

function example_leaf_evaporation_pars(;
    abaxial_vapour_conductance=0.3u"mol/m^2/s",
    adaxial_vapour_conductance=0.0u"mol/m^2/s",
    cuticular_conductance=0.01u"mol/m^2/s",
)
    LeafEvaporationParameters(;
        abaxial_vapour_conductance,
        adaxial_vapour_conductance,
        cuticular_conductance,
    )
end

function example_metabolic_rate_options(;
    respire=true,
    temperature_error_tolerance=1e-3u"K",
    resp_tolerance=1e-5,
)
    SolveMetabolicRateOptions(;
        respire,
        temperature_error_tolerance,
        resp_tolerance,
    )
end

"""
    example_heat_exchange_traits(; kwargs...)

Create example `HeatExchangeTraits` with sensible defaults.
"""
function example_heat_exchange_traits(;
    shape_pars=example_shape_pars(),
    insulation_pars=example_insulation_pars(),
    conduction_pars_external=example_conduction_pars_external(),
    conduction_pars_internal=example_conduction_pars_internal(),
    convection_pars=example_convection_pars(),
    radiation_pars=example_radiation_pars(),
    evaporation_pars=example_evaporation_pars(),
    hydraulic_pars=example_hydraulic_pars(),
    respiration_pars=example_respiration_pars(),
    metabolism_pars=example_metabolism_pars(),
    options=example_metabolic_rate_options(),
)
    HeatExchangeTraits(
        shape_pars,
        insulation_pars,
        conduction_pars_external,
        conduction_pars_internal,
        radiation_pars,
        convection_pars,
        evaporation_pars,
        hydraulic_pars,
        respiration_pars,
        metabolism_pars,
        options,
    )
end

# ── Ectotherm variants ────────────────────────────────────────────────────────

"""
    example_ectotherm_conduction_pars_external(; conduction_fraction=0.1) → ExternalConductionParameters

Create example `ExternalConductionParameters` for an ectotherm with NicheMapR-style defaults.
`conduction_fraction=0.1` corresponds to `pct_cond=10` in NicheMapR's ectotherm model.
"""
function example_ectotherm_conduction_pars_external(;
    conduction_fraction=0.1,
)
    ExternalConductionParameters(; conduction_fraction)
end

"""
    example_ectotherm_metabolism_pars(; kwargs...) → MetabolismParameters

Create example `MetabolismParameters` for an ectotherm.

Defaults to `AndrewsPough2()` (Andrews & Pough 1985 squamate metabolic rate equation).
`metabolic_heat_flow` is not used directly — metabolic heat is computed from the model
given body mass and temperature.
"""
function example_ectotherm_metabolism_pars(;
    metabolic_heat_flow=0.0u"W",
    model=AndrewsPough2(),
)
    MetabolismParameters(; metabolic_heat_flow, model)
end

"""
    example_ectotherm_radiation_pars(; kwargs...) → RadiationParameters

Create example `RadiationParameters` for an ectotherm with NicheMapR defaults.

| Parameter                    | Default        | NicheMapR equivalent |
|------------------------------|----------------|----------------------|
| `body_absorptivity_dorsal`   | 0.85           | alpha_max = 0.85     |
| `body_absorptivity_ventral`  | 0.85           | alpha_min = 0.85     |
| `body_emissivity_dorsal`     | 0.95           | epsilon = 0.95       |
| `body_emissivity_ventral`    | 0.95           | epsilon = 0.95       |
| `sky_view_factor`            | 0.4            | fatosk = 0.4         |
| `ground_view_factor`         | 0.4            | fatosb = 0.4         |
| `bush_view_factor`           | 0.0            | (none)               |
| `ventral_fraction`           | 0.5            | (none)               |
| `solar_orientation`          | `Intermediate()` | postur = 0         |
"""
function example_ectotherm_radiation_pars(;
    body_absorptivity_dorsal  = 0.85,
    body_absorptivity_ventral = 0.85,
    body_emissivity_dorsal    = 0.95,
    body_emissivity_ventral   = 0.95,
    sky_view_factor           = 0.4,
    ground_view_factor        = 0.4,
    bush_view_factor          = 0.0,
    ventral_fraction          = 0.5,
    solar_orientation         = Intermediate(),
)
    RadiationParameters(;
        body_absorptivity_dorsal,
        body_absorptivity_ventral,
        body_emissivity_dorsal,
        body_emissivity_ventral,
        sky_view_factor,
        ground_view_factor,
        bush_view_factor,
        ventral_fraction,
        solar_orientation,
    )
end

"""
    example_ectotherm_evaporation_pars(; kwargs...) → AnimalEvaporationParameters

Create example `AnimalEvaporationParameters` for an ectotherm with NicheMapR defaults.

| Parameter            | Default  | NicheMapR equivalent         |
|----------------------|----------|------------------------------|
| `skin_wetness`       | 0.001    | pct_wet = 0.1 (%)            |
| `eye_fraction`       | 0.0003   | pct_eyes = 0.03 (%)          |
| `bare_skin_fraction` | 1.0      | naked skin (no insulation)   |
| `insulation_wetness` | 0.0      | (no insulation)              |
| `insulation_fraction`| 0.0      | (no insulation)              |
"""
function example_ectotherm_evaporation_pars(;
    skin_wetness        = 0.001,
    eye_fraction        = 0.0003,
    bare_skin_fraction  = 1.0,
    insulation_wetness  = 0.0,
    insulation_fraction = 0.0,
)
    AnimalEvaporationParameters(;
        skin_wetness,
        insulation_wetness,
        eye_fraction,
        bare_skin_fraction,
        insulation_fraction,
    )
end

"""
    example_ectotherm_respiration_pars(; kwargs...) → RespirationParameters

Create example `RespirationParameters` for an ectotherm with NicheMapR defaults.

| Parameter                       | Default  | NicheMapR equivalent                             |
|---------------------------------|----------|--------------------------------------------------|
| `oxygen_extraction_efficiency`  | 0.2      | F_O2 = 20 (%)                                    |
| `respiratory_quotient`          | 0.8      | RQ = 0.8                                         |
| `pant`                          | 1.0      | pantmax = 1 (no panting)                         |
| `exhaled_temperature_offset`    | 0.1 K    | delta_air = 0.1 (°C above Tair)                  |
| `exhaled_relative_humidity`     | 1.0      | (saturated exhaled air)                          |
| `mouth_fraction`                | 0.05     | pct_mouth/100: added to skin_wetness when panting|
"""
function example_ectotherm_respiration_pars(;
    oxygen_extraction_efficiency = 0.2,
    respiratory_quotient         = 0.8,
    pant                         = 1.0,
    exhaled_temperature_offset   = 0.1u"K",
    exhaled_relative_humidity    = 1.0,
    mouth_fraction               = 0.05,
)
    RespirationParameters(; oxygen_extraction_efficiency, respiratory_quotient, pant,
                            exhaled_temperature_offset, exhaled_relative_humidity, mouth_fraction)
end

"""
    example_ectotherm_hydraulic_pars(; kwargs...) → HydraulicParameters

Create example `HydraulicParameters` for an ectotherm with NicheMapR defaults.

| Parameter              | Default        | NicheMapR equivalent          |
|------------------------|----------------|-------------------------------|
| `water_potential`      | -707.0 J/kg    | psi_body = -707               |
| `hydraulic_conductance`| 0.0 …          | (no liquid exchange by default)|
| `specific_hydration`   | 0.000304 …     | (standard value)              |
"""
function example_ectotherm_hydraulic_pars(;
    water_potential      = -707.0u"J/kg",
    hydraulic_conductance = 0.0u"kg / (m^2 * s * (J/kg))",
    specific_hydration   = 0.000304u"m^3 / (m^3 * (J/kg))",
)
    HydraulicParameters(; water_potential, hydraulic_conductance, specific_hydration)
end

"""
    example_ectotherm_conduction_pars_internal(; kwargs...) → InternalConductionParameters

Create example `InternalConductionParameters` for an ectotherm with NicheMapR defaults.

| Parameter             | Default       | NicheMapR equivalent |
|-----------------------|---------------|----------------------|
| `flesh_conductivity`  | 0.5 W/m/K     | k_flesh = 0.5        |
| `fat_fraction`        | 0.0           | rinsul = 0           |
"""
function example_ectotherm_conduction_pars_internal(;
    flesh_conductivity = 0.5u"W/m/K",
    fat_fraction       = 0.0,
    fat_conductivity   = 0.23u"W/m/K",
    fat_density        = 901.0u"kg/m^3",
)
    InternalConductionParameters(; fat_fraction, flesh_conductivity, fat_conductivity, fat_density)
end

"""
    example_ectotherm_heat_exchange_traits(; kwargs...) → HeatExchangeTraits

Create example `HeatExchangeTraits` for an ectotherm using NicheMapR defaults throughout.
All component parameter functions can be overridden individually.
"""
function example_ectotherm_heat_exchange_traits(;
    shape_pars               = DesertIguana(40.0u"g", 1000.0u"kg/m^3"),
    conduction_pars_external = example_ectotherm_conduction_pars_external(),
    conduction_pars_internal = example_ectotherm_conduction_pars_internal(),
    convection_pars          = ConvectionParameters(),
    radiation_pars           = example_ectotherm_radiation_pars(),
    evaporation_pars         = example_ectotherm_evaporation_pars(),
    hydraulic_pars           = example_ectotherm_hydraulic_pars(),
    respiration_pars         = example_ectotherm_respiration_pars(),
    metabolism_pars          = example_ectotherm_metabolism_pars(),
    options                  = example_metabolic_rate_options(),
    insulation_pars          = example_insulation_pars(),
)
    HeatExchangeTraits(
        shape_pars,
        insulation_pars,
        conduction_pars_external,
        conduction_pars_internal,
        radiation_pars,
        convection_pars,
        evaporation_pars,
        hydraulic_pars,
        respiration_pars,
        metabolism_pars,
        options,
    )
end
