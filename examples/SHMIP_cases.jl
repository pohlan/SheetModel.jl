## RUN SHMIP TEST CASES

using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters
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
para_bench = 0.05
surface_val(x, y) = 100(x+200)^(1/4) + 1/60*x - 2e10^(1/4) + 1
f(x, para) = (surface_val(6e3,0) - para*6e3)/6e3^2 * x^2 + para*x
g(y) = 0.5e-6 * abs(y)^3
r(x, para) = (-4.5*x/6e3 + 5) * (surface_val(x, 0) - f(x, para)) /
               (surface_val(x, 0) - f(x, para_bench) + eps())
bed_val(x,y, para) = f(x,para) + g(y) * r(x, para)


function run_SHMIP(test_case; Nx, Ny, make_plot=false, printtime=10^5,
                   γ_ϕ=0.29, γ_h=0.32, dτ_ϕ_=2e6, dτ_h_=2.1)      # parameters for pseudo-transient time stepping

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
    geom  = Dict("sqrt" => (Lx   = 100e3,
                            Ly   = 20e3,
                            surf = (x, y) -> (6 *( sqrt(x+5e3) - sqrt(5e3) ) + 1 ),
                            bed  = (x, y) -> 0.0
                            ),
                 "valley" => (Lx = 6e3,
                             Ly = 1e3,
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
        lx = topo.Lx,  # domain length in x-direction, m
        ly = topo.Ly,  # domain length in y-direction, m
        nx = Nx,
        ny = Ny,

        calc_zs = topo.surf,    # surface elevation, m
        calc_zb = topo.bed,     # bed elevation, m
        calc_m_xyt  = water_input,     # source term, m/s

        ttot = 10 * s_per_day,
        dt   = 3000.0, #  TODO: Adaptive time stepping, in the end it shouldn't be specified as input

        γ_ϕ  = γ_ϕ,  # damping parameter for ϕ
        γ_h  = γ_h,  # damping parameter for h
        dτ_ϕ_ = dτ_ϕ_, # scaling factor for dτ_ϕ
        dτ_h_ = dτ_h_  # scaling factor for dτ_h
    )

    # Initial condition
    @unpack xc, yc, lx, ly, H = input_params
    ϕ0, h0 = S.initial_conditions(
        xc,
        yc,
        H,
        calc_ϕ = (x, y) -> 0.0,
        #calc_ϕ = (x, y) -> exp(- 1e-2*(x-Lx/2)^2) * exp(-1e-2*(yc-Ly/2)^2),
        #calc_ϕ = (x, y) -> rand(),
        calc_h = (x, y) -> 0.04
    )


    N, ϕ, h, nit = S.runthemodel(input_params, ϕ0, h0, printtime=printtime);

    if make_plot
        S.plot_output(xc, yc, N, h)
    end

    return nit, ϕ
end
