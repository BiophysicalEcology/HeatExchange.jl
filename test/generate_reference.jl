# Writes data/reference_solve.csv from the current solve. Run from this directory with the
# test environment active, only when the reference is meant to change.

using CSV, DataFrames
include("reference_cases.jl")

rows = map(REFERENCE_CASES) do case
    merge(case, reference_row(reference_solve(case)))
end
CSV.write(joinpath(@__DIR__, "data", "reference_solve.csv"), DataFrame(rows))
