 # GlaDS SHEET model run

using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters
const S = SheetModel

s_per_day = 3600*24.0
# Define numerical domain and input parameters
input_params = S.Para(
    lx = 60e3,  # domain length in x-direction, m
    ly = 20e3,  # domain length in y-direction, m
    nx = 64,
    ny = 32,

    calc_zs = (x, y)  -> 6.124 * sqrt(x), # surface elevation, m
    calc_zb = (x, y) -> 0.0,              # bed elevation, m
    calc_m_xyt = (x, y, t) -> 0.14/s_per_day - 0.10/(1e-3 *  s_per_day) * 6.124 * sqrt(x) >= 0.0 ? 0.14/s_per_day - 0.10/(1e-3 * s_per_day) * 6.124 * sqrt(x) : 0.0,    # source term, m/s

    ttot = 100.0,
    dt = 100.0
)

@unpack xc, yc, lx, ly = input_params
ϕ0, h0 = S.initial_conditions(
    xc,
    yc,
    #calc_ϕ = (x, y) -> 1e6/lx * x,
    #calc_ϕ = (x, y) -> exp(- 1e-2*(x-Lx/2)^2) * exp(-1e-2*(yc-Ly/2)^2),
    calc_ϕ = (x, y) -> rand(),
    calc_h = (x, y) -> 0.05
)

@time ϕ, h = S.runthemodel(input_params, ϕ0, h0);
S.plot_output(xc, yc, ϕ, h)