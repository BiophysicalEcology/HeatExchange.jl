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

    # The same air as a vapour density gives the same evaporation
    vapour_density = HeatExchange.wet_air_properties(air_temperature, 0.2, 101325.0u"Pa").vapour_density
    vapour_atmos = VapourDensityAtmosphericConditions(vapour_density, 2.0u"m/s", 101325.0u"Pa")
    @test evaporation(leaf_pars, mass_windy, vapour_atmos, area, surface_temperature, air_temperature).transpiration_mass_flow ≈
        out.transpiration_mass_flow
    # Supersaturated air condenses onto a leaf at air temperature, more the more supersaturated
    saturated = HeatExchange.wet_air_properties(air_temperature, 1.0, 101325.0u"Pa").vapour_density
    dew(excess) = evaporation(leaf_pars, mass_windy, VapourDensityAtmosphericConditions(saturated * excess, 2.0u"m/s",
        101325.0u"Pa"), area, air_temperature, air_temperature).transpiration_mass_flow
    @test dew(1.05) < 0.0u"g/s"
    @test dew(1.1) < dew(1.05)
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

@testset "TemperatureDifferenceModel" begin
    body = Body(Ellipsoid(0.2u"g" |> u"kg", 700.0u"kg/m^3", 1.0, 100.0), Naked())
    air_temperature = u"K"(25.0u"°C")
    width = 0.1u"K"
    free(ΔT, temperature_difference_model) = convection(; body, area=0.01u"m^2", air_temperature,
        surface_temperature=air_temperature + ΔT, wind_speed=0.5u"m/s", atmospheric_pressure=101325.0u"Pa",
        fluid=Air(), temperature_difference_model).heat_transfer_coefficient.free
    exact, smooth = ExactTemperatureDifference(), SmoothTemperatureDifference(width)

    # Exact is the default
    @test convection(; body, area=0.01u"m^2", air_temperature, surface_temperature=air_temperature + 1.0u"K",
        wind_speed=0.5u"m/s", atmospheric_pressure=101325.0u"Pa", fluid=Air()).heat_transfer_coefficient.free ==
        free(1.0u"K", exact)
    # At ΔT = 0, smooth free convection is exact free convection at ΔT = width
    @test free(0.0u"K", smooth) ≈ free(width, exact)
    # Smooth has a flat top at ΔT = 0, where exact has a cusp (beyond its 1e-6 K AD floor)
    δ = 1.0e-4u"K"
    @test abs(free(δ, smooth) - free(0.0u"K", smooth)) < 1.0e-6 * free(0.0u"K", smooth)
    @test free(δ, exact) - free(0.0u"K", exact) > 100 * (free(δ, smooth) - free(0.0u"K", smooth))
    # Far from ΔT = 0 the two agree
    @test free(10.0u"K", smooth) ≈ free(10.0u"K", exact) rtol = 1.0e-4
end
