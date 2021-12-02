using Base: Float64, Int64
using LinearAlgebra, Parameters, Statistics, PyPlot, LaTeXStrings

# misc fudge factors to avoid dividing by zero
small = eps(Float64)

"""
Turn all negative numbers into 0.0
"""
pos(x) = x > 0.0 ? x : 0.0

"""
Create a vector of lenght nx where boundary points are zeros and interior points ones.
Used to achieve zero ice thickness (H=0) at ghost points
"""
ghostp(nx) = [0.0; ones(nx-2); 0.0]

"""
Create struct including all model parameters, physical and numerical
"""
@with_kw struct Para{T, F1, F2, F3} @deftype Float64
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
    @assert nx % 32 == 0 && ny % 8 == 0 "nx and ny must be multiples of 32 and 8, respectively."
    dx = (xrange[2]-xrange[1]) / (nx-3)          # grid size
    dy = (yrange[2]-yrange[1]) / (ny-3)
    xc::LinRange{Float64} = LinRange(xrange[1]-dx, xrange[2]+dx, nx) # vector of x-coordinates
    yc::LinRange{Float64} = LinRange(yrange[1]-dy, yrange[2]+dy, ny) # vector of y-coordinates

    # Field parameters (defined on every grid point)
    calc_zs::Function
    calc_zb::Function
    calc_m_phys::F1 # calculate m on the physical domain, f(x, y, t)

    H::T = (ghostp(nx) * ghostp(ny)') .* ( pos.(calc_zs.(xc, yc') .- calc_zb.(xc, yc')) )  # ice thickness, m
    zb::T = calc_zb.(xc, yc')                                                              # bed elevation, m
    calc_m_num::F2 = (ix, iy, t) -> calc_m_phys(xc[ix], yc[iy], t)                         # source term, m/s, calculated from matrix indices f(ix, iy, t)
    # directly create the kernel function as the function calc_m_num cannot be given to kernels at a later stage
    calc_Λ_m!::F3 = @parallel_indices (ix,iy) function calc_Λ_m!(Λ_m, Λ, t)
                        if (ix <= size(Λ_m, 1) && iy <= size(Λ_m, 2))
                            Λ_m[ix, iy] = Λ * calc_m_num(ix, iy, t)
                        end
                        return
                    end

    # Physical time stepping
    ttot        # total simulation time
    dt          # physical time step

    # Pseudo-time iteration
    tol    = 1e-6       # tolerance
    itMax  = 10^6       # max number of iterations
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

"""
Create struct containing output parameters of the model: N, ϕ, h, qx, qy and different errors
"""
@with_kw struct model_output{T}
    N::T
    ϕ::T
    h::T
    qx::T
    qy::T
    Res_ϕ::T
    Res_h::T
    ittot::Int64
    iters::Vector{Int64}
    errs_ϕ::Vector{Float64}
    errs_h::Vector{Float64}
    time_tot::Float64
    T_eff::Float64
end
Broadcast.broadcastable(out::model_output) = Ref(out)

"""
Pre-allocate arrays
"""
function array_allocation(nu::Para)
    @unpack nx, ny = nu
    qx     = @zeros(nx-1, ny)
    qy     = @zeros(nx, ny-1)
    d_eff  = @zeros(nx, ny)
    Λ_m    = @zeros(nx, ny)
    N      = @zeros(nx, ny)
    dϕ_dτ  = @zeros(nx, ny)
    dh_dτ  = @zeros(nx, ny)
    Res_ϕ  = @zeros(nx, ny)
    Res_h  = @zeros(nx, ny)
    return qx, qy, d_eff, Λ_m, N, dϕ_dτ, dh_dτ, Res_ϕ, Res_h
end

"""
Convert input parameters to non-dimensional quantities.
Convert arrays of H and zb to correct types depending on whether CPU or GPU is used.
"""
function format(p::Para, ϕ_init, h_init)
    @unpack g, ρw, k, A, lr, hr,
            H, zb, calc_m_phys, ub,
            dx, dy, xc, yc, xrange, yrange,
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
    xy_ = max(xrange[2]-xrange[1], yrange[2]-yrange[1])

    H_ = mean(H)
    zb_ = H_ / r_ρ
    m_ = mean(calc_m_phys(xc, yc', 0.0)) # for time-dependent input: temporal peak
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
    scaled_params = reconstruct(Para, p,
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
        H = Data.Array(H ./ H_),
        zb = Data.Array(zb ./ zb_),
        calc_m_num = (ix, iy, t) -> calc_m_num(ix, iy, t) / m_,

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
    ϕ_init = Data.Array(ϕ_init ./ ϕ_)
    h_init = Data.Array(h_init ./ h_)

    return scaled_params, ϕ_init, h_init, ϕ_, N_, h_, q_
end

mat(x) = Matrix{Float64}(x) # shorter notation for use in reformat() function

"""
Convert unitless output parameters back to dimensional quantities.
Convert all arrays to Matrix{Float64} type.
"""
function reformat(output::model_output, N_, ϕ_, h_, q_)
    @unpack N, ϕ, h, qx, qy, Res_ϕ, Res_h = output
    output_descaled = reconstruct(model_output{Matrix{Float64}}, output,
        N  = mat(N .* N_),
        ϕ  = mat(ϕ .* ϕ_),
        h  = mat(h .* h_),
        qx = mat(qx .* q_),
        qy = mat(qy .* q_),
        Res_ϕ = mat(Res_ϕ .* ϕ_),
        Res_h = mat(Res_h .* h_),
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
    ϕ_init = calc_ϕ.(xc, yc')
    h_init = calc_h.(xc, yc')
    return ϕ_init, h_init
end

function plot_output(xc, yc, H, N, h, qx, qy, Res_ϕ, Res_h, iters, errs_h, errs_ϕ)
    pygui(true)

    # (I) ϕ and h
    x_plt = [xc[1]; xc .+ (xc[2]-xc[1])]
    y_plt = [yc[1]; yc .+ (yc[2]-yc[1])]
    N[H .== 0.0] .= NaN
    h[H .== 0.0] .= NaN

    figure()
    # (Ia) pcolor of ϕ and h fields
    subplot(2, 2, 1)
    pcolor(x_plt, y_plt, h')#, edgecolors="black")
    colorbar()
    title("h")
    subplot(2, 2, 2)
    pcolor(x_plt, y_plt, N')#, edgecolors="black")
    colorbar()
    title("N")
    # (Ib) cross-sections of ϕ and h
    subplot(2, 2, 3)
    ind = size(N,2)÷2
    plot(xc, h[:, ind])
    title(join(["h at y = ", string(round(yc[ind], digits=1))]))
    subplot(2, 2, 4)
    plot(xc, N[:, ind])
    title(join(["N at y = ", string(round(yc[ind], digits=1))]))

    # (II) fluxes
    # don't show any value outside of glacier domain
    qx[H[1:end-1, :] .== 0.] .= NaN
    qx[H[2:end, :]   .== 0.] .= NaN
    qy[H[:, 1:end-1] .== 0.] .= NaN
    qy[H[:, 2:end]   .== 0.] .= NaN

    figure()
    subplot(1, 2, 1)
    pcolor(qx')
    colorbar()
    title("qx (m/s)")
    subplot(1, 2, 2)
    pcolor(qy')
    colorbar()
    title("qy (m/s)")

    # (III) residual fields
    Res_ϕ[H .== 0.0] .= NaN
    Res_h[H .== 0.0] .= NaN

    figure()
    subplot(1, 2, 1)
    pcolormesh(Res_h')
    colorbar()
    title("err_h")
    subplot(1, 2, 2)
    pcolormesh(Res_ϕ')
    colorbar()
    title("err_ϕ")

    # (IV) iteration vs. error
    figure()
    semilogy(iters, errs_ϕ, label=L"\mathrm{Res}_ϕ", color="gold")
    semilogy(iters, errs_h, label=L"\mathrm{Res}_h", color="deepskyblue")
    title("errors: norm(...) / length(..)")

    xlabel(L"# iterations $i$")
    ylabel("error")
    legend()
end