using Base: Float64
using SheetModel, Test
const S = SheetModel

# ----------------------#
# Test the model output #
# ----------------------#

# read the reference ϕ and h fields
include("test_references.jl")

# run the cases
include(joinpath(@__DIR__, "../examples/SHMIP_cases.jl"))

ϕ_test = Dict()
h_test = Dict()
δx     = Dict()
δy     = Dict()

# A1 test case: sqrt geometry
inputs, outputs = run_SHMIP(test_case="A1", nx=128, ny=128, dt = 1e9,
                            γ_ϕ= 0.9, γ_h=0.8, dτ_ϕ_=1.0, dτ_h_= 6e-6)
ϕ_test["A1"] = outputs.ϕ
h_test["A1"] = outputs.h
δx["A1"]     = fld((inputs.input_params.nx-3), (size(ϕ_ref["A1"], 1) - 1))
δy["A1"]     = fld((inputs.input_params.ny-3), (size(ϕ_ref["A1"], 2) - 1))

# E1 test case: valley geometry
#inputs, outputs = run_SHMIP(test_case="E1", nx=256, ny=256, dt = 1e9,
#                            γ_ϕ= 0.7, γ_h=0.5, dτ_ϕ_= 1.0, dτ_h_= 7e-4)
#ϕ_test["E1"] = outputs.ϕ
#h_test["E1"] = outputs.h
#δx["E1"]     = fld((inputs.input_params.nx-3), (size(ϕ_ref["A1"], 1) - 1))
#δy["E1"]     = fld((inputs.input_params.ny-3), (size(ϕ_ref["A1"], 2) - 1))

# test whether the runs agree with the references
@testset "Model runs" begin
    for test_case in keys(ϕ_test)
        @test all((ϕ_test[test_case][2:δx[test_case]:end-1, 2:δy[test_case]:end-1] .- ϕ_ref[test_case]) ./ (ϕ_ref[test_case] .+ eps(Float64)) .< 1e-2) # add an eps(Float64) in case ϕ_ref is zero
        @test all((h_test[test_case][2:δx[test_case]:end-1, 2:δy[test_case]:end-1] .- h_ref[test_case]) ./ (h_ref[test_case] .+ eps(Float64)) .< 1e-2)
    end
end
