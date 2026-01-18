using HeatExchange
using Unitful
using Test

@testset "Metabolic Rate Equations" begin
    @testset "AndrewsPough2 (squamates)" begin
        eq = AndrewsPough2()
        mass = 100.0u"g"
        T = 30.0u"°C"

        # Basic output
        Q = metabolic_rate(eq, mass, T)
        @test Q isa Unitful.Power
        @test Q > 0.0u"W"

        # Temperature capping (1-50°C range)
        @test metabolic_rate(eq, mass, 0.0u"°C") == metabolic_rate(eq, mass, 1.0u"°C")
        @test metabolic_rate(eq, mass, 60.0u"°C") == metabolic_rate(eq, mass, 50.0u"°C")

        # Metabolic rate increases with temperature and mass
        @test metabolic_rate(eq, mass, 30.0u"°C") > metabolic_rate(eq, mass, 20.0u"°C")
        @test metabolic_rate(eq, 100.0u"g", T) > metabolic_rate(eq, 10.0u"g", T)

        # Resting (M4=1) > standard (M4=0)
        @test metabolic_rate(AndrewsPough2(M4=1.0), mass, T) > metabolic_rate(AndrewsPough2(M4=0.0), mass, T)

        # Accepts Kelvin
        @test metabolic_rate(eq, mass, 303.15u"K") ≈ metabolic_rate(eq, mass, 30.0u"°C")
    end

    @testset "Kleiber (mammals)" begin
        Q = metabolic_rate(Kleiber(), 1.0u"kg")
        @test Q isa Unitful.Power
        @test Q > 0.0u"W"

        # Metabolic rate increases with mass (allometric scaling)
        @test metabolic_rate(Kleiber(), 10.0u"kg") > metabolic_rate(Kleiber(), 1.0u"kg")

        # T_body argument is ignored
        @test metabolic_rate(Kleiber(), 1.0u"kg", 37.0u"°C") == metabolic_rate(Kleiber(), 1.0u"kg")
    end

    @testset "McKechnieWolf (birds)" begin
        Q = metabolic_rate(McKechnieWolf(), 100.0u"g")
        @test Q isa Unitful.Power
        @test Q > 0.0u"W"

        # Metabolic rate increases with mass
        @test metabolic_rate(McKechnieWolf(), 1000.0u"g") > metabolic_rate(McKechnieWolf(), 100.0u"g")

        # T_body argument is ignored
        @test metabolic_rate(McKechnieWolf(), 100.0u"g", 40.0u"°C") == metabolic_rate(McKechnieWolf(), 100.0u"g")
    end
end

@testset "O2-Joules Conversion" begin
    Q = 1.0u"W"
    V = 100.0u"ml/hr"
    rq = 0.8

    @testset "Typical" begin
        # Round-trip conversion
        V_converted = Joules_to_O2(Typical(), Q, rq)
        Q_back = O2_to_Joules(Typical(), V_converted, rq)
        @test Q ≈ Q_back

        # Default (no equation specified)
        @test Joules_to_O2(Q) ≈ Joules_to_O2(Typical(), Q, 1.0)
        @test O2_to_Joules(V) ≈ O2_to_Joules(Typical(), V, 1.0)
    end

    @testset "Kleiber1961" begin
        # Round-trip conversion
        V_converted = Joules_to_O2(Kleiber1961(), Q, rq)
        Q_back = O2_to_Joules(Kleiber1961(), V_converted, rq)
        @test Q ≈ Q_back
    end
end
