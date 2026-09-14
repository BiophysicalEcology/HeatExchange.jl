using HeatExchange
using BiophysicalGeometry
using Test
using Unitful
using UnitfulMoles

@testset "leaf evaporation" begin
    # Representative mass transfer coefficients (m/s) for a leaf-sized body.
    # These are supplied from convection() in practice; we construct them directly here
    # to keep the test independent of the pre-existing convection() gas_fractions keyword bug.
    mass_windy = TransferCoefficients(0.005u"m/s", 0.001u"m/s", 0.004u"m/s")
    mass_calm  = TransferCoefficients(0.001u"m/s", 0.001u"m/s", 0.0u"m/s")

    air_temperature      = u"K"(30.0u"°C")
    surface_temperature  = u"K"(35.0u"°C")
    area                 = 0.01u"m^2"
    atmos = AtmosphericConditions(0.2, 2.0u"m/s", 101325.0u"Pa")

    leaf_pars = LeafEvaporationParameters(;
        abaxial_vapour_conductance = 0.3u"mol/m^2/s",
        adaxial_vapour_conductance = 0.0u"mol/m^2/s",
        cuticular_conductance      = 0.01u"mol/m^2/s",
    )

    out = evaporation(leaf_pars, mass_windy, atmos, area, surface_temperature, air_temperature)

    # Positive heat loss when leaf is warmer and drier than air
    @test out.evaporation_heat_flow > 0.0u"W"
    @test out.transpiration_mass_flow > 0.0u"g/s"

    # Stomata fully closed → much less transpiration than open stomata
    leaf_pars_closed = LeafEvaporationParameters(;
        abaxial_vapour_conductance = 0.0u"mol/m^2/s",
        adaxial_vapour_conductance = 0.0u"mol/m^2/s",
        cuticular_conductance      = 0.01u"mol/m^2/s",
    )
    out_closed = evaporation(leaf_pars_closed, mass_windy, atmos, area, surface_temperature, air_temperature)
    @test out_closed.transpiration_mass_flow < out.transpiration_mass_flow
    @test out_closed.evaporation_heat_flow > 0.0u"W"

    # Less boundary layer transport (calmer conditions) → less evaporation
    out_calm = evaporation(leaf_pars, mass_calm, atmos, area, surface_temperature, air_temperature)
    @test out_calm.evaporation_heat_flow < out.evaporation_heat_flow

    # Negative water potential reduces surface humidity → less evaporation
    out_dry = evaporation(leaf_pars, mass_windy, atmos, area, surface_temperature, air_temperature;
        water_potential = -1e6u"J/kg")
    @test out_dry.evaporation_heat_flow < out.evaporation_heat_flow
end

@testset "SimpleForcedCorrelation" begin
    body = Body(Ellipsoid(0.2u"g" |> u"kg", 700.0u"kg/m^3", 1.0, 100.0), Naked())
    common = (; body, area=0.01u"m^2", air_temperature=u"K"(25.0u"°C"),
        surface_temperature=u"K"(30.0u"°C"), atmospheric_pressure=101325.0u"Pa", fluid=Air())

    windy = convection(; common..., wind_speed=2.0u"m/s", correlation=SimpleForcedCorrelation())
    @test isfinite(ustrip(u"W", windy.convection_flow))
    @test windy.convection_flow > 0.0u"W"

    # No free-convection term: zero wind gives zero convective flow.
    calm = convection(; common..., wind_speed=0.0u"m/s", correlation=SimpleForcedCorrelation())
    @test calm.convection_flow == 0.0u"W"

    # Nusselt correlation, not shape, drives forced convection: same
    # Reynolds number gives the same Nu regardless of shape argument.
    @test nusselt_forced(SimpleForcedCorrelation(), Plate(1.0u"kg", 1.0u"kg/m^3", 1.0, 1.0), 500.0) ==
        nusselt_forced(SimpleForcedCorrelation(), Ellipsoid(1.0u"kg", 1.0u"kg/m^3", 1.0, 1.0), 500.0)
    @test nusselt_forced(SimpleForcedCorrelation(), Plate(1.0u"kg", 1.0u"kg/m^3", 1.0, 1.0), 500.0) ≈ 0.6 * sqrt(500.0)
    @test nusselt_free(SimpleForcedCorrelation(), Plate(1.0u"kg", 1.0u"kg/m^3", 1.0, 1.0), 1.0e5, 0.71) == 0.0

    # ShapeCorrelation (default) includes free convection -- differs from
    # SimpleForcedCorrelation even at the same wind speed.
    shaped = convection(; common..., wind_speed=2.0u"m/s")
    @test shaped.convection_flow != windy.convection_flow
end
