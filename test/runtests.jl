using SheetModel, Test
const S = SheetModel

# -------------------------#
# Test the small functions #
# -------------------------#

@testset "Helper functions" begin
    @test S.calc_q(0.5, 0.0, 0.0, 1.0, 1.25, 1.5, 1e-25) == 0.0 # shouldn't produce NaNs
    @test S.calc_q(0.5, 1.0, 1.0, 1.0, 1.25, 1.5, 1e-25) == -0.35355339059327373
    @test S.calc_pw(1.0, 1.0, 1.0, 0.5) == 0.5
    @test S.calc_N(0.2, 0.9, 1.0, 1.0, 0.5, 0.1) == 0.35
    @test S.calc_vc(0.2, 0.5, 0.9, 1.0, 1.0, 0.5, 0.1, 3.0, 1.0) == 0.0015879629629629627
    @test S.calc_vo(1.2, 1.0, 1.0, 1.0) == 0.0
    @test S.calc_vo(0.6, 1.0, 1.0, 1.0) == 0.4
end


# ----------------------#
# Test the model output #
# ----------------------#

# read the reference ϕ and h fields
include("test_references.jl")

# run the cases
include(joinpath(@__DIR__, "../examples/SHMIP_cases.jl"))
nx, ny = 64, 32
ϕ_test = Dict()
h_test = Dict()
nit_test = Dict()

for test_case in keys(ϕ_ref)
    inputs, outputs = run_SHMIP(test_case, Nx=nx, Ny=ny, dt = 5e7)
    ϕ_test[test_case] = outputs.ϕ
    h_test[test_case] = outputs.h
end

# test whether the runs agree with the references
@testset "Model runs" begin
    for test_case in keys(ϕ_ref)
        @test ϕ_test[test_case][1:20:end, 1:10:end] ≈ ϕ_ref[test_case]
        @test h_test[test_case][1:20:end, 1:10:end] ≈ h_ref[test_case]
    end
end
