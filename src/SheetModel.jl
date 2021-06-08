module SheetModel

using Base: Float64, Int64
using LinearAlgebra, Parameters, Statistics, PyPlot

export Para

@with_kw struct Para @deftype Float64
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
    lx         # domain size
    ly
    nx::Int64         # number of grids
    ny::Int64
    dx = lx/nx      # grid size
    dy = ly/ny
    xc::LinRange = LinRange(0.0, lx, nx) # vector of x-coordinates
    yc::LinRange = LinRange(0.0, ly, ny) # vector of y-coordinates

    # Field parameters (defined on every grid point)
    calc_H::Function
    calc_zb::Function
    calc_m::Function
    H::Array = calc_H.(xc, yc')  # ice thickness, m
    zb::Array = calc_zb.(xc, yc')  # bed elevation, m
    m::Array = calc_m.(xc, yc')  # source term, m/s

    # Physical time stepping
    ttot        # total simulation time
    dt          # physical time step

    # Pseudo-time iteration
    tol    = 1e-6       # tolerance
    itMax  = 60       # max number of iterations
    damp   = 1-41/nx    # damping (this is a tuning parameter, dependent on e.g. grid resolution) # TODO: define how?
    dτ     = (1.0/(dx^2/k/2.1) + 1.0/dt)^-1 # pseudo-time step; TODO: define how?

    # Dimensionless numbers
    Σ   = NaN
    Γ   = NaN
    Λ   = NaN
    r_ρ = ρw / ρi
end
Broadcast.broadcastable(p::Para) = Ref(p)

"""
Converts input parameters to non-dimensional quantities
"""
function scaling(p::Para, ϕ0, h0)
    @unpack g, ρw, k, A, lr, hr,
            H, zb, m, ub,
            lx, ly, dx, dy, xc, yc,
            ttot, dt, dτ,
            r_ρ, α, β, n,
            Σ, Γ, Λ = p

    # Scaling factors
    g_ = g
    ρ_ = ρw
    k_ = k
    A_ = A
    lr_ = lr
    h_ = hr
    xy_ = max(lx, ly)

    H_ = mean(H)
    zb_ = H_
    m_ = mean(m)
    ub_ = ub

    ϕ_ = g_ * H_ * ρ_ / r_ρ
    q_ = k_ * h_^α * ( ϕ_  / xy_ )^(β-1)
    t_ = xy_ * h_ / q_
    vo_ = ub_ * h_ / lr_
    vc_ = A_ * h_ * ϕ_ ^n

    if any(.!isnan.([Σ, Γ, Λ]))
        @warn "Σ, Γ and Λ have already been assigned."
    end

    # Dimensionless parameters
    scaled_params = Para(p,
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
        m = m ./ m_,

        # Numerical domain
        lx = lx ./ xy_,
        ly = ly ./ xy_,
        dx = dx ./ xy_,
        dy = dy ./ xy_,
        xc = xc ./ xy_,
        yc = yc ./ xy_,

        ttot = ttot / t_,
        dt = dt / t_,
        dτ = dτ / t_,

        Σ = vo_ * xy_ / q_,
        Γ = vc_ * xy_ / q_,
        Λ = m_ * xy_ / q_
        )

    # variables
    ϕ0 = ϕ0 ./ ϕ_
    h0 = h0 ./ h_

    return scaled_params, ϕ0, h0
end

"""
Returns arrays of initial conditions for ϕ and h
"""
function initial_conditions(xc, yc; calc_ϕ = (x,y) -> 0.0, calc_h = (x,y) -> 0.0)
    ϕ0 = calc_ϕ.(xc, yc')
    h0 = calc_h.(xc, yc')
    return ϕ0, h0
end

"""
Allocates zeros()
"""
function array_allocation(nu::Para)
    @unpack nx, ny = nu
    vo     = zeros(nx, ny)
    vc     = zeros(nx, ny)
    dϕ_dx  = zeros(nx-1, ny)
    dϕ_dy  = zeros(nx, ny-1)
    qx     = zeros(nx-1, ny)
    qy     = zeros(nx, ny-1)
    dϕ_dτ  = zeros(nx-2, ny-2)
    dh_dτ  = zeros(nx-2, ny-2)
    Res_ϕ  = zeros(nx-2, ny-2)
    Res_h  = zeros(nx-2, ny-2)
    return vo, vc, dϕ_dx, dϕ_dy, qx, qy, dϕ_dτ, dh_dτ, Res_ϕ, Res_h
end


"""
Calculates discharge
"""
function calc_q(h, dϕ_du, p::Para) # u can be x or y
    @unpack k, α, β = p
    q = - k * h^α * abs(dϕ_du)^(β-2) * dϕ_du
    if isnan(q)
        q = 0.0     # change NaNs to zero; (if dϕ_du is zero, (dϕ_du)^(β-2) is Inf -> Inf*0 -> NaN)
    end
    return q
end

"""
Calculates water pressure
"""
function calc_pw(ϕ, ρw, g, zb)
    #@unpack ρw, g, zb = p
    return ϕ - ρw * g * zb
end

"""
Calculates effective pressure
"""
function calc_N(ϕ, ρi, ρw, g, H, zb)
    #@unpack ρi, g, H = p
    pw = calc_pw(ϕ, ρw, g, zb)
    return ρi * g * H - pw
end

"""
Calculates closure rate
"""
function calc_vc(ϕ, h, ρi, ρw, g, H, zb, n, A)
    #@unpack n, A = p
    N = calc_N(ϕ, ρi, ρw, g, H, zb)
    return 2 / n^n * A * h * abs(N)^(n-1) * N
end

"""
Calculates opening rate
"""
function calc_vo(h, ub, hr, lr)
    #@unpack ub, hr, lr = p
    if h < hr
        vo = ub * (hr .- h[h .< hr]) ./ lr
    else
        vo = 0.0
    end
    return vo
end



@views function runthemodel(input::Para, ϕ0, h0)

    params, ϕ0, h0 = scaling(input, ϕ0, h0)
    @unpack ev, g, ρw, ρi, n, A, Σ, Γ, Λ, m, dx, dy, xc, yc, H, zb, ub, hr, lr, dt, ttot, tol, itMax, damp, dτ = params

    # Array allocation
    vo, vc, dϕ_dx, dϕ_dy, qx, qy, dϕ_dτ, dh_dτ, Res_ϕ, Res_h = array_allocation(params)

    ϕ0[1, :]  = ρw .* g .* zb[1, :] # boundary condition, zero water pressure
    ϕ_old   = copy(ϕ0)
    ϕ       = copy(ϕ0)
    h_old   = copy(h0)
    h       = copy(h0)

    t = 0.0; it = 0; ittot = 0

    # Physical time loop
    while t<ttot
        iter = 0; err_ϕ = 2*tol; err_h = 2*tol;

        # Pseudo-transient iteration
        while max(err_ϕ, err_h) >tol && iter<itMax
            dϕ_dx   .= diff(ϕ, dims=1) ./ dx    # hydraulic gradient
            dϕ_dy   .= diff(ϕ, dims=2) ./ dy
            qx      .= calc_q.(h[1:end-1, :], dϕ_dx, params)  # TODO: which h values to take?!
            qy      .= calc_q.(h[:, 1:end-1], dϕ_dy, params)
            qx[end, :]      .= 0.0  # boundary condition
            qy[:, [1, end]] .= 0.0  # boundary condition

            vo     .= calc_vo.(h, ub, hr, lr)     # opening rate
            vc     .= calc_vc.(ϕ, h, ρi, ρw, g, H, zb, n, A)  # closure rate

            # calculate residuals
            Res_ϕ       .=      .- ev/(ρw*g) * (ϕ[2:end-1, 2:end-1] .- ϕ_old[2:end-1, 2:end-1])/dt .-         # dhe/dt
                                (diff(qx, dims=1)[:, 2:end-1]/dx .+ diff(qy, dims=2)[2:end-1, :]/dy) .-      # div(q)
                                (Σ * vo[2:end-1, 2:end-1] .- Γ * vc[2:end-1, 2:end-1])            .+         # dh/dt
                                Λ * m[2:end-1, 2:end-1]                                                                        # source term
            Res_h       .=    - (h[2:end-1, 2:end-1] .- h_old[2:end-1, 2:end-1]) / dt  .+
                                (Σ * vo[2:end-1, 2:end-1] .- Γ * vc[2:end-1, 2:end-1])

            # damped rate of change
            dϕ_dτ      .= Res_ϕ .+ damp .* dϕ_dτ
            dh_dτ      .= Res_h .+ damp .* dh_dτ

            # update fields
            ϕ[2:end-1, 2:end-1]  .= ϕ[2:end-1, 2:end-1] .+ dτ * dϕ_dτ   # update ϕ
            h[2:end-1, 2:end-1]  .= h[2:end-1, 2:end-1] .+ dτ * dh_dτ   # updat h

            err_ϕ = norm(Res_ϕ)/length(Res_ϕ)
            err_h   = norm(Res_h)/length(Res_h)
            iter += 1
        end
        ittot += iter; it += 1; t += dt

        ϕ_old .= ϕ
        h_old .= h
    end

    # TODO: convert back to dimensional quantities

    return ϕ, h
end

function plot_output(xc, yc, ϕ, h)
    pygui(true)
    figure()
    subplot(1, 2, 1)
    pcolor(xc, yc, h')
    colorbar()
    title("h")
    subplot(1, 2, 2)
    pcolor(xc, yc, ϕ')
    colorbar()
    title("ϕ")
end



end # module
