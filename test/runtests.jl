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
@safetestset "geometry" begin include("geometry.jl") end
@safetestset "environment" begin include("environment.jl") end
@safetestset "organism" begin include("geometry.jl") end
