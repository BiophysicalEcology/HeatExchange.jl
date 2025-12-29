using HeatExchange
using SafeTestsets
using Aqua

Aqua.test_unbound_args(HeatExchange)
Aqua.test_stale_deps(HeatExchange)
Aqua.test_undefined_exports(HeatExchange)
Aqua.test_project_extras(HeatExchange)
Aqua.test_deps_compat(HeatExchange)
Aqua.test_project_toml_formatting(HeatExchange)

@safetestset "ectotherm" begin include("ectotherm.jl") end
@safetestset "endotherm" begin include("endotherm.jl") end
@safetestset "ellipsoid" begin include("ellipsoid.jl") end
