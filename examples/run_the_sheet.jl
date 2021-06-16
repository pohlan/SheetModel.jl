## RUN SHMIP TEST CASES

using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters
const S = SheetModel

# water input, source term
# Suite A: use different steady and spatially uniform water inputs
runoff = (A1 = (x, y) -> 7.93e-11,
          A2 = (x, y) -> 1.59e-9,
          A3 = (x, y) -> 5.79e-9,
          A4 = (x, y) -> 2.5e-8,
          A5 = (x, y) -> 4.5e-8,
          A6 = (x, y) -> 5.79e-7,
          )

# surface elevation
surf = (sqrt   = (x, y) -> (6 *( sqrt(x+5e3) - sqrt(5e3) ) + 1 ),
        # valley = ...
)

# bed elevation
bed = (sqrt   = (x, y) -> 0.0,
       # valley = ...
)

# Define numerical domain and input parameters
input_params = S.Para(
    lx = 100e3,  # domain length in x-direction, m
    ly = 20e3,  # domain length in y-direction, m
    nx = 64,
    ny = 32,

    calc_zs = surf.sqrt,    # surface elevation, m
    calc_zb = bed.sqrt,     # bed elevation, m
    calc_m = runoff.A1,     # source term, m/s

    ttot = 25000.0,
    dt = 2500.0
)

# Initial condition
@unpack xc, yc, lx, ly = input_params
ϕ0, h0 = S.initial_conditions(
    xc,
    yc,
    calc_ϕ = (x, y) -> 1e6/lx * x,
    #calc_ϕ = (x, y) -> exp(- 1e-2*(x-Lx/2)^2) * exp(-1e-2*(yc-Ly/2)^2),
    #calc_ϕ = (x, y) -> rand(),
    calc_h = (x, y) -> 0.04
)

@time ϕ, h = S.runthemodel(input_params, ϕ0, h0, printit=1000);
S.plot_output(xc, yc, ϕ, h)
