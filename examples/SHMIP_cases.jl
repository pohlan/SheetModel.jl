## RUN SHMIP TEST CASES

using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters
using PyPlot
const S = SheetModel

s_per_day  = 24 * 60 * 60
s_per_year = 365 * s_per_day

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

# for "valley" geometry
const para_bench = 0.05
surface_val(x, y) = 100(x+200)^(1/4) + 1/60*x - 2e10^(1/4) + 1
f(x, para) = (surface_val(6e3,0) - para*6e3)/6e3^2 * x^2 + para*x
g(y) = 0.5e-6 * abs(y)^3
r(x, para) = (-4.5*x/6e3 + 5) * (surface_val(x, 0) - f(x, para)) /
               (surface_val(x, 0) - f(x, para_bench) + eps())
bed_val(x,y, para) = f(x,para) + g(y) * r(x, para)


function run_SHMIP(test_case; Nx, Ny, make_plot=false, printtime=10^5, printit=10^5,
                   γ_ϕ=0.8, γ_h=0.9, dτ_ϕ_=1.0, dτ_h_=7e-6)      # parameters for pseudo-transient time stepping

    # suite A: use different steady and spatially uniform water inputs
    runoff = Dict("A1" => (x, y, t) -> 7.93e-11,
                  "A2" => (x, y, t) -> 1.59e-9,
                  "A3" => (x, y, t) -> 5.79e-9,
                  "A4" => (x, y, t) -> 2.5e-8,
                  "A5" => (x, y, t) -> 4.5e-8,
                  "A6" => (x, y, t) -> 5.79e-7,
                  # suite E
                  "E1"  => (x, y, t) -> 1.158e-6,
                  "E2"  => (x, y, t) -> 1.158e-6,
                  "E3"  => (x, y, t) -> 1.158e-6,
                  "E4"  => (x, y, t) -> 1.158e-6,
                  "E5"  => (x, y, t) -> 1.158e-6
              )

    # suite D: seasonally changing water input with constant shift DT
    DT = Dict(# suite D (sqrt)
              "D1" => -4.0,
              "D2" => -2.0,
              "D3" => 0.0,
              "D4" => 2.0,
              "D5" => 4.0,
              # suite F (valley)
              "F1" => -6.0,
              "F2" => -3.0,
              "F3" => 0.0,
              "F4" => 3.0,
              "F5" => 6.0
              )

    para = Dict(# suite E
                "E1" => 0.05,
                "E2" => 0.0,
                "E3" => -0.1,
                "E4" => -0.5,
                "E5" => -0.7,
                # suite F
                "F1" => 0.05,
                "F2" => 0.05,
                "F3" => 0.05,
                "F4" => 0.05,
                "F5" => 0.05
                )

    # geometry: "sqrt" -> ice-sheet margin; "valley" -> valley glacier
    geom  = Dict("sqrt" => (xrange   = (0.0, 100e3),
                            yrange   = (0.0, 20e3),
                            surf = (x, y) -> (6 *( sqrt(x+5e3) - sqrt(5e3) ) + 1 ),
                             bed  = (x, y) -> 0.0
                            # bed = (x, y) -> x * 1e-3 * ((y-10e3)*2e-4)^2 # alternative bed topography varying in y-direction; no convergence
                            ),
                 "valley" => (xrange = (0.0, 6e3),
                              yrange = (-500.0, 500.0),
                              surf = (x, y) -> surface_val(x, y),
                              bed = (x, y) -> bed_val(x, y, para[test_case])
                              )
    )

    if any(startswith.(test_case, ["A", "D"]))
        topo = geom["sqrt"]
    elseif any(startswith.(test_case, ["E", "F"]))
        topo = geom["valley"]
    end

    if any(startswith.(test_case, ["A", "E"]))
        water_input = runoff[test_case]
    elseif any(startswith.(test_case, ["D", "F"]))
        water_input = make_runoff_fct(topo.surf, DT[test_case])
    end

    # Define numerical domain and input parameters
    input_params = S.Para(
        xrange = topo.xrange,  # domain length in x-direction, m
        yrange = topo.yrange,  # domain length in y-direction, m
        nx = Nx,
        ny = Ny,

        calc_zs = topo.surf,    # surface elevation, m
        calc_zb = topo.bed,     # bed elevation, m
        calc_m_xyt  = water_input,     # source term, m/s

        #ttot = 1000.0,
        ttot = 2000, #*s_per_day,
        dt   = 2000, #*s_per_day, #  TODO: Adaptive time stepping, in the end it shouldn't be specified as input

        itMax = 5*10^4,
        γ_ϕ  = γ_ϕ,  # damping parameter for ϕ
        γ_h  = γ_h,  # damping parameter for h
        dτ_ϕ_ = dτ_ϕ_, # scaling factor for dτ_ϕ
        dτ_h_ = dτ_h_  # scaling factor for dτ_h
    )

    # Initial condition
    @unpack xc, yc, H = input_params
    ϕ0, h0 = S.initial_conditions(
        xc,
        yc,
        H,
        calc_ϕ = (x, y) -> 100.0,
        #calc_ϕ = (x, y) -> 1e6/lx * x,
        #calc_ϕ = (x, y) -> exp(- 1e-2*(x-Lx/2)^2) * exp(-1e-2*(yc-Ly/2)^2),
        #calc_ϕ = (x, y) -> rand(),
        calc_h = (x, y) -> 0.04
    )


    model_output = S.runthemodel(input_params, ϕ0, h0, printtime=printtime, printit=printit);
    @unpack N, ϕ, h, qx, qy, ittot, Err_ϕ, Err_h, errs_ϕ, errs_h, Res_ϕ, Res_h, qx_ice, qy_ice = model_output

    if make_plot
        S.plot_output(xc, yc, H, N, h, qx, qy, qx_ice, qy_ice, Err_ϕ, Err_h, errs_h, errs_ϕ)
    end

    return (;input_params, ϕ0, h0), model_output
end
