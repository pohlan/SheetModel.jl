## RUN SHMIP TEST CASES (to be run with "run_SHMIP.jl")

using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters, Infiltrator
using PyPlot
const S = SheetModel

const day  = 24 * 60 * 60
const year = 365day

# water input function for suite D
function make_runoff_fct(zs, DT)
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

# Return arrays of initial conditions for ϕ and h
function initial_conditions(xc, yc; calc_ϕ = (x,y) -> 0.0, calc_h = (x,y) -> 0.01)
    ϕ_init = calc_ϕ.(xc, yc')
    h_init = calc_h.(xc, yc')
    return ϕ_init, h_init
end


# Turn all negative numbers into 0.0
pos(x) = x > 0.0 ? x : 0.0

# Determine dx or dy from Lx and dx or Ly and dy
grid_size(Lx, nx) = Lx / (nx-3)

# Set boundaries (ghostpoints) to zero.
function ghost(A)
    A[[1, end], :] .= 0.
    A[:, [1, end]] .= 0.
    return A
end

function plot_output(xc, yc, H, ϕ, h, qx, qy, Res_ϕ, Res_h, iters, errs_h, errs_ϕ)
    pygui(true)

    # (I) ϕ and h
    x_plt = [xc[1]; xc .+ (xc[2]-xc[1])]
    y_plt = [yc[1]; yc .+ (yc[2]-yc[1])]
    ϕ[H .== 0.0] .= NaN
    h[H .== 0.0] .= NaN

    figure()
    # (Ia) pcolor of ϕ and h fields
    subplot(2, 2, 1)
    pcolor(x_plt, y_plt, h')#, edgecolors="black")
    colorbar()
    title("h")
    subplot(2, 2, 2)
    pcolor(x_plt, y_plt, ϕ')#, edgecolors="black")
    colorbar()
    title("ϕ")
    # (Ib) cross-sections of ϕ and h
    subplot(2, 2, 3)
    ind = size(ϕ ,2)÷2
    plot(xc, h[:, ind])
    title(join(["h at y = ", string(round(yc[ind], digits=1))]))
    subplot(2, 2, 4)
    plot(xc, ϕ[:, ind])
    title(join(["ϕ at y = ", string(round(yc[ind], digits=1))]))

    # (II) fluxes
    # don't show any value outside of glacier domain
    # qx[H[1:end-1, :] .== 0.] .= NaN
    # qx[H[2:end, :]   .== 0.] .= NaN
    # qy[H[:, 1:end-1] .== 0.] .= NaN
    # qy[H[:, 2:end]   .== 0.] .= NaN

    # figure()
    # subplot(1, 2, 1)
    # pcolor(qx')
    # colorbar()
    # title("qx (m/s)")
    # subplot(1, 2, 2)
    # pcolor(qy')
    # colorbar()
    # title("qy (m/s)")

    # (III) residual fields
    # Res_ϕ[H .== 0.0] .= NaN
    # Res_h[H .== 0.0] .= NaN

    # figure()
    # subplot(1, 2, 1)
    # pcolormesh(Res_h')
    # colorbar()
    # title("err_h")
    # subplot(1, 2, 2)
    # pcolormesh(Res_ϕ')
    # colorbar()
    # title("err_ϕ")

    # (IV) iteration vs. error
    figure()
    semilogy(iters, errs_ϕ, label=L"\mathrm{Res}_ϕ", color="gold")
    semilogy(iters, errs_h, label=L"\mathrm{Res}_h", color="deepskyblue")
    title("errors: norm(...) / length(..)")

    xlabel(L"# iterations $i$")
    ylabel("error")
    legend()
end

function run_SHMIP(;test_case="A1", nx=64, ny=32, itMax=10^7, tol=1e-6, make_plot=false, do_print=true, warmup=0,
                   dt=1e20day, tsteps=1, γ_ϕ= 0.8, γ_h=0.8, dτ_ϕ_=1.0, dτ_h=6e-6)      # parameters for pseudo-transient time stepping

    # suite A: use different steady and spatially uniform water inputs
    runoff = Dict("A1" => (x, y, t) -> 7.93e-11,
                  "A2" => (x, y, t) -> 1.59e-9,
                  "A3" => (x, y, t) -> 5.79e-9,
                  "A4" => (x, y, t) -> 2.5e-8,
                  "A5" => (x, y, t) -> 4.5e-8,
                  "A6" => (x, y, t) -> 5.79e-7,
                  # suite E
                  "E1"  => (x, y, t) -> 1.158e-6, # 5e-10
                  "E2"  => (x, y, t) -> 1.158e-6, # 5e-10
                  "E3"  => (x, y, t) -> 1.158e-6, # 5e-10
                  "E4"  => (x, y, t) -> 1.158e-6, # 5e-10
                  "E5"  => (x, y, t) -> 1.158e-6  # 5e-10
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
                            ),
                 "valley" => (xrange = (0.0, 6e3),
                              yrange = (-500.0, 500.0),
                              surf = (x, y) -> surface_val(x, y),
                              bed = (x, y) -> bed_val(x, y, para[test_case])
                              )
    )

    if any(startswith.(test_case, ["A", "D"]))
        topo = geom["sqrt"]
        ev = 0
    elseif any(startswith.(test_case, ["E", "F"]))
        topo = geom["valley"]
        ev = 1e-3
    end

    if any(startswith.(test_case, ["A", "E"]))
        water_input = runoff[test_case]
    elseif any(startswith.(test_case, ["D", "F"]))
        water_input = make_runoff_fct(topo.surf, DT[test_case])
    end

    # input parameters
    x1, xend          = topo.xrange
    y1, yend          = topo.yrange
    Lx                = xend - x1
    Ly                = yend - y1
    dx                = grid_size(Lx, nx)
    dy                = grid_size(Ly, ny)
    xc                = LinRange(x1-dx, xend+dx, nx)         # including ghost points
    yc                = LinRange(y1-dy, yend+dy, ny)
    zb                = ghost(topo.bed.(xc, yc'))
    H                 = ghost(pos.(topo.surf.(xc, yc') .- zb))
    ttot              = tsteps * dt

    # definition of function calc_m (calculating the source term)
    # (doing calc_m(ix, iy, t) =  water_input(xc[ix], yc[iy], t) directly gives an error as xc and yc are not accessible anymore later)
    function make_calc_m(xc, yc, water_input)
        calc_m(ix, iy, t) = water_input(xc[ix], yc[iy], t)
        return calc_m
    end
    calc_m = make_calc_m(xc, yc, water_input)

    # initial conditions
    ϕ_init = 100. * ones(nx, ny)
    h_init = 0.04 * ones(nx, ny)

    # masks for ice domain and boundary conditions
    ice_mask    = H .> 0.
    bc_diric    = falses(nx, ny); bc_diric[2, :] .= true
    bc_no_xflux = diff(ice_mask, dims=1) .!= 0
    bc_no_yflux = diff(ice_mask, dims=2) .!= 0

    # struct of input parameters
    params_struct = S.model_input(;H, zb, Lx, Ly, dx, dy, ttot, dt, itMax, tol, γ_ϕ, γ_h, dτ_ϕ_, dτ_h, ev)

    # call the SheetModel
    input = (;params_struct, ϕ_init, h_init, calc_m, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux, do_print, warmup)
    output = runthemodel(;input...);

    # plotting
    @unpack N, ϕ, h, qx, qy,
            ittot, iters, Res_ϕ, Res_h, errs_ϕ, errs_h = output

    if make_plot && !any(isnan.([errs_ϕ[end], errs_h[end]]))
        plot_output(xc, yc, H, ϕ, h, qx, qy, Res_ϕ, Res_h, iters, errs_h, errs_ϕ)
    end

    return (;input..., test_case), output
end
