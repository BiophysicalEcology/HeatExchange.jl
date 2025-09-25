using HeatExchange
using SafeTestsets
using Aqua

Aqua.test_unbound_args(HeatExchange)
Aqua.test_stale_deps(HeatExchange)
Aqua.test_undefined_exports(HeatExchange)
Aqua.test_project_extras(HeatExchange)
Aqua.test_deps_compat(HeatExchange)
Aqua.test_project_toml_formatting(HeatExchange)

#@safetestset "biophysics" begin include("biophysics.jl") end
#@safetestset "geometry" begin include("geometry.jl") end
#@safetestset "environment" begin include("environment.jl") end
#@safetestset "organism" begin include("geometry.jl") end
@safetestset "plate shape comparisons" begin include("test_geometry_NicheMapR.jl") end
@safetestset "cylinder shape comparisons" begin include("test_geometry_NicheMapR.jl") end
@safetestset "ellipsoid shape comparisons" begin include("test_geometry_NicheMapR.jl") end
@safetestset "iguana shape comparisons" begin include("test_geometry_NicheMapR.jl") end
@safetestset "frog shape comparisons" begin include("test_geometry_NicheMapR.jl") end
