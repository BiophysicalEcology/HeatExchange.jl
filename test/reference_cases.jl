# A grid of endotherm solves wider than the NicheMapR comparisons, to pin the current results
# while the solve is rewritten. `generate_reference.jl` writes them to data/reference_solve.csv.

using HeatExchange, BiophysicalGeometry, Unitful

const REFERENCE_SHAPES = (;
    cylinder = Cylinder(; mass = 2.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 3.0),
    sphere = Sphere(; mass = 2.0u"kg", density = 1000.0u"kg/m^3"),
    plate = Plate(; mass = 2.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 2.0, axis_ratio_c = 5.0),
    ellipsoid = Ellipsoid(; mass = 2.0u"kg", density = 1000.0u"kg/m^3", axis_ratio_b = 2.0, axis_ratio_c = 2.0),
)

const REFERENCE_CASES = [
    (; shape, fur, fat, sun, wind, contact, asymmetric)
    for shape in keys(REFERENCE_SHAPES), fur in (false, true), fat in (0.0, 0.2), sun in (false, true),
        wind in (0.1, 3.0), contact in (0.0, 0.3), asymmetric in (false, true)
    if fur || !asymmetric
]

function reference_solve(case)
    dorsal_depth = case.fur ? (case.asymmetric ? 10e-3u"m" : 5e-3u"m") : 0.0u"m"
    ventral_depth = case.fur ? (case.asymmetric ? 2e-3u"m" : 5e-3u"m") : 0.0u"m"
    insulation_pars = example_insulation_pars(;
        insulation_depth_dorsal = dorsal_depth, insulation_depth_ventral = ventral_depth,
        insulation_depth_compressed = ventral_depth)
    shape = REFERENCE_SHAPES[case.shape]
    fur = FibrousLayer((dorsal_depth + ventral_depth) / 2, 30e-6u"m", 3000e4u"1/m^2")
    fat = FatLayer(case.fat, 901.0u"kg/m^3")
    traits = example_heat_exchange_traits(;
        shape_pars = shape,
        insulation_pars,
        conduction_pars_external = example_conduction_pars_external(; conduction_fraction = case.contact),
        conduction_pars_internal = example_conduction_pars_internal(; fat_fraction = case.fat),
    )
    organism = Organism(Body(shape, CompositeInsulation(fur, fat)), traits)
    environment_vars = example_environment_vars(;
        wind_speed = case.wind * u"m/s",
        global_radiation = case.sun ? 800.0u"W/m^2" : 0.0u"W/m^2",
        air_temperature = u"K"(15.0u"°C"),
    )
    environment = (; environment_pars = example_environment_pars(), environment_vars)
    solve_metabolic_rate(organism, environment, u"K"(33.0u"°C"), u"K"(20.0u"°C"))
end

# The results of a solve as plain numbers in SI units.
function reference_row(out)
    t = out.thermoregulation
    e = out.energy_flows
    m = out.mass_flows
    g = out.morphology
    si(x) = Float64(ustrip(upreferred(x)))
    sides = map((:dorsal, :ventral)) do side
        ts = getproperty(t, side)
        es = getproperty(e, side)
        NamedTuple{Symbol.(string(side, "_") .* ("skin_temperature", "insulation_temperature", "solar",
            "convection", "conduction", "net_metabolic", "skin_evaporation", "longwave"))}(
            si.((ts.skin_temperature, ts.insulation_temperature, es.solar, es.convection, es.conduction,
                es.net_metabolic, es.skin_evaporation, es.longwave)))
    end
    merge((;
        core_temperature = si(t.core_temperature), skin_temperature = si(t.skin_temperature),
        insulation_temperature = si(t.insulation_temperature), lung_temperature = si(t.lung_temperature),
        solar_flow = si(e.solar_flow), longwave_flow_in = si(e.longwave_flow_in),
        longwave_flow_out = si(e.longwave_flow_out), metabolic_heat_flow = si(e.metabolic_heat_flow),
        evaporation_heat_flow = si(e.evaporation_heat_flow), respiration_heat_flow = si(e.respiration_heat_flow),
        convection_heat_flow = si(e.convection_heat_flow), conduction_flow = si(e.conduction_flow),
        air_flow = si(m.air_flow), m_evap = si(m.m_evap), m_sweat = si(m.m_sweat),
        total_area = si(g.total_area), area_convection = si(g.area_convection),
        area_silhouette = si(g.area_silhouette), characteristic_dimension = si(g.characteristic_dimension),
    ), sides...)
end
