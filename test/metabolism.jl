using HeatExchange
using Unitful
using Test

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
