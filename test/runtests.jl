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

# repr(ϕ[1:20:end, 1:20:end]) and repr(h[1:20:end, 1:20:end])
# ttot = ..., dt = ...
ϕ_ref = Dict("A1"   =>  [],
             "A2"   =>  [],
             "A3"   =>  [],
             "A4"   =>  [],
             "A5"   =>  [],
             "A6"   =>  [])

include(joinpath(@__DIR__, "../examples/SHMIP_cases.jl"))
nx, ny = 64, 32
ϕ_test = Dict()
h_test = Dict()
nit_test = Dict()
test_cases = ["A1", "A2", "A3", "A4", "A5", "A6"] #,
              #"D1", "D2", "D3", "D4", "D5",
              #"E1", "E2", "E3", "E4", "E5",
              #"F1", "F2", "F3", "F4", "F5"]
for test_case in test_cases
    nit, ϕ, h, rest = run_SHMIP(test_case, Nx=nx, Ny=ny)
    ϕ_test[test_case] = ϕ
    h_test[test_case] = h
end

@testset "Model runs" begin
    for test_case in test_cases
        @test ϕ_test[test_case][1:20:end, 1:20:end] ≈ ϕ_ref[test_case]["ϕ"]
        @test h_test[test_case][1:20:end, 1:20:end] ≈ ϕ_ref[test_case]["h"]
    end
end
