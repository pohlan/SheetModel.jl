## RUN SHMIP TEST CASES

using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters
const S = SheetModel


function run_SHMIP(test_case; Nx, Ny, t_tot, make_plot=false)

    # suite A: use different steady and spatially uniform water inputs
    runoff = Dict("A1" => (x, y, t) -> 7.93e-11,
                  "A2" => (x, y, t) -> 1.59e-9,
                  "A3" => (x, y, t) -> 5.79e-9,
                  "A4" => (x, y, t) -> 2.5e-8,
                  "A5" => (x, y, t) -> 4.5e-8,
                  "A6" => (x, y, t) -> 5.79e-7,
              )

    # suite D: seasonally changing water input with constant shift DT
    DT = Dict("D1" => -4.0,
              "D2" => -2.0,
              "D3" => 0.0,
              "D4" => 2.0,
              "D5" => 4.0
              )

    # water input function for suite D
    function make_runoff_fct(zs, DT)
        year  = 31536000.0   # number of seconds per year
        lapse = -0.0075      # lapse rate, K/m
        DDF   = 0.01 / 86400 # degree day factor, m/(K*s)
        basal = 7.93e-11     # basal melt rate, m/s, equal to the source of scenario A1
        temp(t) = -16*cos(2*pi/year*t)- 5 + DT
        runoff(zs,t) = max(0.0, (zs * lapse .+ temp(t)) * DDF) .+ basal
        return (x, y, t) -> runoff(zs(x, y), t)
    end

    # geometry: "sqrt" -> ice-sheet margin; "valley" -> valley glacier
    geom  = Dict("sqrt" => (Lx   = 100e3,
                            Ly   = 20e3,
                            surf = (x, y) -> (6 *( sqrt(x+5e3) - sqrt(5e3) ) + 1 ),
                            bed  = (x, y) -> 0.0
                            ),
                 "valley" => (Lx = 6e3,
                             Ly = 1e3,
                             # surf = ...
                             # bed = ...
                             )
    )

    if any(startswith.(test_case, ["A", "D"]))
        topo = geom["sqrt"]
    else
        topo = geom["valley"]
    end

    if startswith(test_case, "A")
        water_input = runoff[test_case]
    elseif startswith(test_case, "D")
        water_input = make_runoff_fct(topo.surf, DT[test_case])
    end

    # Define numerical domain and input parameters
    input_params = S.Para(
        lx = topo.Lx,  # domain length in x-direction, m
        ly = topo.Ly,  # domain length in y-direction, m
        nx = Nx,
        ny = Ny,

        calc_zs = topo.surf,    # surface elevation, m
        calc_zb = topo.bed,     # bed elevation, m
        calc_m_xyt  = water_input,     # source term, m/s

        ttot = t_tot,
        dt = 2500.0 #  TODO: Adaptive time stepping, in the end it shouldn't be specified as input
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

    @time N, h = S.runthemodel(input_params, ϕ0, h0, printit=1000);

    if make_plot
        S.plot_output(xc, yc, N, h)
    end

    return nothing
end

run_SHMIP("D1", Nx=64, Ny=32, t_tot=2500.0, make_plot=true);
