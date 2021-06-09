using SheetModel, Test
const S = SheetModel

# -------------------------#
# Test the small functions #
# -------------------------#

@testset "Helper functions" begin
    @test S.upstream.([1,2], [1,2]', 3, [1 2; 3 -4]) == [1 4; 2 6]
    @test S.calc_q(0.5, 0.0, 1.0, 1.25, 1.5, 1e-25) == 0.0 # shouldn't produce NaNs
    @test S.calc_q(0.5, 1.0, 1.0, 1.25, 1.5, 1e-25) == -0.42044820762685725
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
onetimestep = Dict("ϕ" => [0.0 0.0;
                           4.050057825824443e6 4.050057825824443e6;
                           7.964048824264514e6 7.964048824264514e6;
                           1.0035820159457687e7 1.0035820159457687e7],
                 "h"   => [0.0 0.0;
                           0.015873015873015872 0.015874689059231187;
                           0.031746031746031744 0.03174747527941701;
                           0.047619047619047616 0.0476127369421051])

include(joinpath(@__DIR__, "../run_the_sheet.jl"))
lx = 100e3
nx = 64
@testset "Model runs" begin
    @test xc == LinRange(0.0, lx, nx)
    @test ϕ[1:20:end, 1:20:end] ≈ onetimestep["ϕ"]
    @test h[1:20:end, 1:20:end] ≈ onetimestep["h"]
end
