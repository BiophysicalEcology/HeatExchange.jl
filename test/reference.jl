# The solve reproduces data/reference_solve.csv, recorded by generate_reference.jl.

using CSV, DataFrames, Test
include("reference_cases.jl")

reference = CSV.read(joinpath(@__DIR__, "data", "reference_solve.csv"), DataFrame)

for (case, expected) in zip(REFERENCE_CASES, eachrow(reference))
    @testset "$case" begin
        @test all(k -> string(getproperty(expected, k)) == string(getproperty(case, k)), keys(case))
        row = reference_row(reference_solve(case))
        for k in keys(row)
            @test row[k] ≈ expected[k] rtol = 1e-6 atol = 1e-12
        end
    end
end
