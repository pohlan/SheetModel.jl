using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters
const S = SheetModel


# Define numerical domain and input parameters
input_params = S.Para(
    # SHMIP A3
    lx = 100e3,  # domain length in x-direction, m
    ly = 20e3,  # domain length in y-direction, m
    nx = 64,
    ny = 32,

    calc_H = (x, y)  -> (6 *( sqrt(x+5e3) - sqrt(5e3) ) + 1 ), # ice thickness, m
    calc_zb = (x, y) -> 0.0,        # bed elevation, m
    calc_m = (x, y) -> 5.79e-09,    # source term, m/s

    ttot = 100.0,
    dt = 100.0
)

# Initial condition
@unpack xc, yc, lx, ly = input_params
ϕ0, h0 = S.initial_conditions(
    xc,
    yc,
    calc_ϕ = (x, y) -> 1e6/lx * x,
    #calc_ϕ = (x, y) -> exp(- 1e-2*(x-Lx/2)^2) * exp(-1e-2*(yc-Ly/2)^2),
    #calc_ϕ = (x, y) -> rand(),
    calc_h = (x, y) -> 0.05/lx * x
)

@time ϕ, h = S.runthemodel(input_params, ϕ0, h0, printit=100);
S.plot_output(xc, yc, ϕ, h)
