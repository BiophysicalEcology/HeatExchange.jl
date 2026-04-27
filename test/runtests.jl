using HeatExchange
using Aqua
using SafeTestsets
using Test

@testset "Quality assurance" begin
    Aqua.test_unbound_args(HeatExchange)
    Aqua.test_stale_deps(HeatExchange)
    Aqua.test_undefined_exports(HeatExchange)
    Aqua.test_project_extras(HeatExchange)
    Aqua.test_deps_compat(HeatExchange)
end

@safetestset "biophysics" begin include("biophysics.jl") end
@safetestset "ectotherm" begin include("ectotherm.jl") end
@safetestset "leaf" begin include("leaf.jl") end
@safetestset "endotherm" begin include("endotherm.jl") end
@safetestset "heat_balance" begin include("heat_balance.jl") end
@safetestset "ellipsoid" begin include("ellipsoid.jl") end
@safetestset "metabolism" begin include("metabolism.jl") end
