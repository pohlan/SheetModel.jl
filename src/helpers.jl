using Base: Float64, Int64
using LinearAlgebra, Parameters, Statistics, PyPlot

# misc fudge factors to avoid dividing by zero
small = eps(Float32)

"""
Turn all negative numbers into 0.0
"""
pos(x) = x > 0.0 ? x : 0.0

"""
Creates a vector of lenght nx where boundary points are zeros and interior points ones.
Used to achieve zero ice thickness (H=0) at ghost points
"""
ghostp(nx) = [0.0; ones(nx-2); 0.0]

"""
Create struct including all model parameters, physical and numerical
"""
@with_kw struct Para{F1, F2} @deftype Float64
    # Scalars (one global value)
    g     = 9.81              # gravitational acceleration, m/s^2
    ρw    = 1000.0            # water density, kg/m^3
    ρi    = 910.0             # ice density, kg/m^3
    α     = 1.25              # first sheet flow exponent
    β     = 1.5               # second sheet flow exponent
    k     = 0.005             # sheet conductivity, m^(7/4)kg^(-1/2)
    n     = 3.0               # Glen's flow law exponent
    A     = 3.375e-24         # ice flow constant, Pa^(-n)s^(-1)
    ev    = 0.0               # englacial void ratio; SHMIP: 0 for ice-sheet, 1e-3 for valley glacier
    lr    = 2.0               # horizontal cavity spacing, m
    hr    = 0.1               # bedrock bump height, m
    ub    = 1e-6              # basal sliding speed, m/s

    # Numerical domain
    xrange::Tuple{Float64, Float64}   # domain size
    yrange::Tuple{Float64, Float64}
    nx::Int64                         # number of grid points, including ghost points where ice thickness = 0
    ny::Int64
    dx = (xrange[2]-xrange[1]) / (nx-3)          # grid size
    dy = (yrange[2]-yrange[1]) / (ny-3)
    xc::LinRange{Float64} = LinRange(xrange[1]-dx, xrange[2]+dx, nx) # vector of x-coordinates
    yc::LinRange{Float64} = LinRange(yrange[1]-dy, yrange[2]+dy, ny) # vector of y-coordinates

    # Field parameters (defined on every grid point)
    calc_zs::Function
    calc_zb::Function
    calc_m_xyt::F1 # f(x, y, t)

    H::Matrix{Float64} = (ghostp(nx) * ghostp(ny)') .* ( pos.(calc_zs.(xc, yc') .- calc_zb.(xc, yc')) )  # ice thickness, m
    zb::Matrix{Float64} = calc_zb.(xc, yc')                    # bed elevation, m
    calc_m_t::F2 = t -> calc_m_xyt.(xc, yc', t)                # source term, m/s, f(t)

    # Physical time stepping
    ttot        # total simulation time
    dt          # physical time step

    # Pseudo-time iteration
    tol    = 1e-6       # tolerance
    itMax  = 5*10^3     # max number of iterations
    γ_ϕ    = 1e-3       # damping parameter for ϕ update
    γ_h    = 0.8        # damping parameter for h update
    dτ_ϕ_   = 1e6       # scaling factor for dτ_ϕ
    dτ_h_   = 50.0      # scaling factor for dτ_h

    # Dimensionless numbers
    Σ   = NaN
    Γ   = NaN
    Λ   = NaN
    r_ρ = ρw / ρi
end
Broadcast.broadcastable(p::Para) = Ref(p)

@with_kw struct model_output
    N::Matrix{Float64}
    ϕ::Matrix{Float64}
    h::Matrix{Float64}
    qx::Matrix{Float64}
    qy::Matrix{Float64}
    Err_ϕ::Matrix{Float64}
    Err_h::Matrix{Float64}
    Res_ϕ::Matrix{Float64}
    Res_h::Matrix{Float64}
    ittot::Int64
    iters::Vector{Int64}
    errs_ϕ::Vector{Float64}
    errs_h::Vector{Float64}
    errs_ϕ_rel::Vector{Float64}
    errs_h_rel::Vector{Float64}
    errs_ϕ_res::Vector{Float64}
    errs_h_res::Vector{Float64}
    errs_ϕ_resrel::Vector{Float64}
    errs_h_resrel::Vector{Float64}
end
Broadcast.broadcastable(out::model_output) = Ref(out)

"""
Pre-allocate arrays
"""
function array_allocation(nu::Para)
    @unpack nx, ny = nu
    Δϕ     = zeros(nx, ny)
    Δh     = zeros(nx, ny)
    qx     = zeros(nx-1, ny)
    qy     = zeros(nx, ny-1)
    d_eff  = zeros(nx-2, ny-2)
    m      = zeros(nx, ny)
    N      = zeros(nx, ny)
    dϕ_dτ  = zeros(nx, ny)
    dh_dτ  = zeros(nx, ny)
    Res_ϕ  = zeros(nx, ny)
    Res_h  = zeros(nx, ny)
    return Δϕ, Δh, qx, qy, d_eff, m, N, dϕ_dτ, dh_dτ, Res_ϕ, Res_h
end

"""
Convert input parameters to non-dimensional quantities
"""
function scaling(p::Para, ϕ0, h0)
    @unpack g, ρw, k, A, lr, hr,
            H, zb, calc_m_t, ub,
            dx, dy, xc, yc,
            ttot, dt,
            r_ρ, α, β, n,
            Σ, Γ, Λ = p

    # Scaling factors
    g_ = g
    ρ_ = ρw
    k_ = k
    A_ = A
    lr_ = lr
    h_ = hr
    xy_ = max(xc[end]-xc[1], yc[end]-yc[1])

    H_ = mean(H)
    zb_ = H_ / r_ρ
    m_ = mean(calc_m_t(0.0)) # for time-dependent input: temporal peak
    ub_ = ub

    ϕ_ = g_ * H_ * ρ_ / r_ρ
    N_ = ϕ_
    q_ = k_ * h_^α * ( ϕ_  / xy_ )^(β-1)
    t_ = xy_ * h_ / q_
    vo_ = ub_ * h_ / lr_
    vc_ = A_ * h_ * N_ ^n

    if any(.!isnan.([Σ, Γ, Λ]))
        @warn "Σ, Γ and Λ have already been assigned."
    end

    # Dimensionless parameters
    scaled_params = reconstruct_to(p, Para,
        # Scalars (one global value)
        g = g / g_,
        ρw = ρw / ρ_,
        ρi = ρw / ρ_ / r_ρ,
        k = k / k_,
        A = A / A_,
        lr = lr / lr_,
        hr = hr / h_,
        ub = ub / ub_,

        # Field parameters (defined on every grid point)
        H = H ./ H_,
        zb = zb ./ zb_,
        calc_m_t = t -> calc_m_t(t) ./ m_,

        # Numerical domain
        dx = dx ./ xy_,
        dy = dy ./ xy_,
        xc = xc ./ xy_,
        yc = yc ./ xy_,

        ttot = ttot / t_,
        dt = dt / t_,

        Σ = vo_ * xy_ / q_,
        Γ = vc_ * xy_ / q_,
        Λ = m_ * xy_ / q_
        )::Para

    # variables
    ϕ0 = ϕ0 ./ ϕ_
    h0 = h0 ./ h_

    return scaled_params, ϕ0, h0, ϕ_, N_, h_, q_
end

function descaling(output::model_output, N_, ϕ_, h_, q_)
    @unpack N, ϕ, h, qx, qy,
            Err_ϕ, Err_h , Res_ϕ, Res_h = output
    output_descaled = model_output(output,
        N = N .* N_,
        ϕ = ϕ .* ϕ_,
        h = h .* h_,
        qx = qx .* q_,
        qy = qy .* q_,
        Err_ϕ = Err_ϕ .* ϕ_,
        Err_h = Err_h .* h_,
        Res_ϕ = Res_ϕ .* ϕ_,
        Res_h = Res_h .* h_,
        )::model_output
    return output_descaled
end


#----------------------------------#
# Functions used in SHMIP_cases.jl #
#----------------------------------#
"""
Return arrays of initial conditions for ϕ and h
"""
function initial_conditions(xc, yc, H; calc_ϕ = (x,y) -> 0.0, calc_h = (x,y) -> 0.01)
    ϕ0 = calc_ϕ.(xc, yc')
    h0 = calc_h.(xc, yc')
    return ϕ0, h0
end

function plot_output(xc, yc, H, N, h, qx, qy, Err_ϕ, Err_h,
                     iters, errs_h, errs_ϕ, errs_ϕ_rel, errs_h_rel,
                     errs_ϕ_res, errs_h_res, errs_ϕ_resrel, errs_h_resrel)
    x_plt = [xc[1]; xc .+ (xc[2]-xc[1])]
    y_plt = [yc[1]; yc .+ (yc[2]-yc[1])]
    N[H .== 0.0] .= NaN
    h[H .== 0.0] .= NaN
    pygui(true)
    # pcolor of ϕ and h fields
    figure()
    subplot(2, 2, 1)
    pcolor(x_plt, y_plt, h')#, edgecolors="black")
    colorbar()
    title("h")
    subplot(2, 2, 2)
    pcolor(x_plt, y_plt, N')#, edgecolors="black")
    colorbar()
    title("N")
    # cross-sections of ϕ and h
    subplot(2, 2, 3)
    ind = size(N,2)÷2
    plot(xc, h[:, ind])
    title(join(["h at y = ", string(round(yc[ind], digits=1))]))
    subplot(2, 2, 4)
    plot(xc, N[:, ind])
    title(join(["N at y = ", string(round(yc[ind], digits=1))]))

    # don't show any value outside of glacier domain
    qx_plot = copy(qx)
    qy_plot = copy(qy)
    qx_plot[H[1:end-1, :] .== 0.] .= NaN
    qx_plot[H[2:end, :]   .== 0.] .= NaN
    qy_plot[H[:, 1:end-1] .== 0.] .= NaN
    qy_plot[H[:, 2:end]   .== 0.] .= NaN

    figure()
    subplot(1, 2, 1)
    pcolor(qx_plot')
    colorbar()
    title("qx (m/s)")
    subplot(1, 2, 2)
    pcolor(qy_plot')
    colorbar()
    title("qy (m/s)")

    Err_ϕ[H .== 0.0] .= NaN
    Err_h[H .== 0.0] .= NaN
    figure()
    subplot(1, 2, 1)
    pcolormesh(Err_h')
    colorbar()
    title("err_h")
    subplot(1, 2, 2)
    pcolormesh(Err_ϕ')
    colorbar()
    title("err_ϕ")

    figure()
    semilogy(iters, errs_ϕ, label="err_ϕ", color="darkorange")
    semilogy(iters, errs_h, label="err_h", color="darkblue")
    semilogy(iters, errs_ϕ_rel, label="relative err_ϕ", color="darkorange", linestyle=":")
    semilogy(iters, errs_h_rel, label="relative err_h", color="darkblue", linestyle=":")
    semilogy(iters, errs_ϕ_res, label="err_ϕ_res", color="gold")
    semilogy(iters, errs_h_res, label="err_h_res", color="deepskyblue")
    semilogy(iters, errs_ϕ_resrel, label="relative err_ϕ_res", color="gold", linestyle=":")
    semilogy(iters, errs_h_resrel, label="relative err_h_res", color="deepskyblue", linestyle=":")

    xlabel("# iterations")
    ylabel("error")
    legend()
end