module SheetModel

#using Core: ScalarIndex
using Base: Float64
using LinearAlgebra, Parameters, Statistics, PyPlot

export Para

@with_kw struct Para{F<:Float64, I<:Int64, Arr<:Array{Float64}, Lin<:LinRange{Float64}}
    # Scalars (one global value)
    g::F     = 9.81              # gravitational acceleration, m/s^2
    ρw::F    = 1000.0            # water density, kg/m^3
    ρi::F    = 910.0             # ice density, kg/m^3
    α::F     = 1.25              # first sheet flow exponent
    β::F     = 1.5               # second sheet flow exponent
    k::F     = 0.005             # sheet conductivity, m^(7/4)kg^(-1/2)
    n::F     = 3.0               # Glen's flow law exponent
    A::F     = 3.375e-24         # ice flow constant, Pa^(-n)s^(-1)
    ev::F    = 0.0               # englacial void ratio; SHMIP: 0 for ice-sheet, 1e-3 for valley glacier
    lr::F    = 2.0               # horizontal cavity spacing, m
    hr::F    = 0.1               # bedrock bump height, m
    ub::F    = 1e-6              # basal sliding speed, m/s

    # Field parameters (defined on every grid point)
    H::Arr  # ice thickness, m
    zb::Arr # bed elevation, m
    m::Arr  # source term, m/s

    # Numerical domain
    lx::F         # domain size
    ly::F
    nx::I         # number of grids
    ny::I
    dx::F = lx/nx      # grid size
    dy::F = ly/ny
    xc::Lin = LinRange(0.0, lx, nx) # vector of x-coordinates
    yc::Lin = LinRange(0.0, ly, ny) # vector of y-coordinates

    # Physical time stepping
    ttot::F        # total simulation time
    dt::F          # physical time step

    # Pseudo-time iteration
    tol::F    = 1e-6       # tolerance
    itMax::I  = 60       # max number of iterations
    damp::F   = 1-41/nx    # damping (this is a tuning parameter, dependent on e.g. grid resolution) # TODO: define how?
    dτ::F     = (1.0/(dx^2/k/2.1) + 1.0/dt)^-1 # pseudo-time step; TODO: define how?

    # Dimensionless numbers
    Σ::F   = NaN
    Γ::F   = NaN
    Λ::F   = NaN
    r_ρ::F = ρw / ρi
end

function scaling(p::Para, ϕ0, h0)
    @unpack g, ρw, k, A, lr, hr,
            H, zb, m, ub,
            lx, ly, dx, dy, xc, yc,
            ttot, dt, dτ,
            r_ρ, α, β, n = p

    # Scaling factors
    g_ = g
    ρ_ = ρw
    k_ = k
    A_ = A
    lr_ = lr
    h_ = hr
    xy_ = max(lx, ly) # ?

    H_ = mean(H)
    zb_ = H_ / r_ρ # ?
    m_ = mean(m)
    ub_ = ub

    ϕ_ = g_ * H_ * ρ_ / r_ρ
    q_ = k_ * h_^α * ( ϕ_  / xy_ )^(β-1)
    t_ = xy_ * h_ / q_
    vo_ = ub_ * h_ / lr_
    vc_ = A_ * h_ * ϕ_ ^n

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
        #dτ = dτ / t_, # ?

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
Allocates zeros()
"""
function array_allocation(nu::Para)
    @unpack nx, ny = nu
    qx     = zeros(nx-1, ny-1)
    qy     = zeros(nx-1, ny-1)
    dϕ_dτ  = zeros(nx-2, ny-2)
    dh_dτ  = zeros(nx-2, ny-2)
    Res_ϕ  = zeros(nx-2, ny-2)
    Res_h  = zeros(nx-2, ny-2)
    return qx, qy, dϕ_dτ, dh_dτ, Res_ϕ, Res_h
end


"""
Calculates discharge
"""
function calc_q(h, dϕ_du, p::Para) # u can be x or y
    @unpack k, α, β = p
    q = - k .* h.^α .* abs.(dϕ_du).^(β-2) .* dϕ_du
    q[isnan.(q)] .= 0.0  # if dϕ_du is zero, Inf *
    return q
end

"""
Calculates water pressure
"""
function calc_pw(ϕ, p::Para)
    @unpack ρw, g, zb = p
    return ϕ .- ρw .* g .* zb
end

"""
Calculates effective pressure
"""
function calc_N(ϕ, p::Para)
    @unpack ρi, g, H = p
    pw = calc_pw(ϕ, p)
    return ρi * g * H .- pw
end

"""
Calculates closure rate
"""
function calc_vc(ϕ, h, p::Para)
    @unpack n, A = p
    N = calc_N(ϕ, p)
    return 2 / n^n * A * h .* abs.(N).^(n-1) .* N
end

"""
Calculates opening rate
"""
function calc_vo(h, p::Para)
    @unpack ub, hr, lr = p
    vo = zeros(size(h))
    vo[h .< hr] = ub * (hr .- h[h .< hr]) ./ lr
    vo[h .>= hr] .= 0.0
    return vo
end



function runthemodel(input::Para, ϕ0, h0)

    params, ϕ0, h0 = scaling(input, ϕ0, h0)
    @unpack ev, g, ρw, Σ, Γ, Λ, m, dx, dy, xc, yc, dt, ttot, tol, itMax, damp, dτ = params

    # Array allocation
    qx, qy, dϕ_dτ, dh_dτ, Res_ϕ, Res_h = array_allocation(params)

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
            dϕ_dx, dϕ_dy = diff(ϕ, dims=1) ./ dx, diff(ϕ, dims=2) ./ dy   # hydraulic gradient
            qx, qy     = calc_q(h[1:end-1, :], dϕ_dx, params), calc_q(h[:, 1:end-1], dϕ_dy, params) # TODO: which h values to take?!

            vo     = calc_vo(h, params)                   # opening rate
            vc     = calc_vc(ϕ, h, params)  # closure rate

            # calculate residuals
            Res_ϕ        =      ev/(ρw*g) * (ϕ[2:end-1, 2:end-1] .- ϕ_old[2:end-1, 2:end-1])/dt .+    # dhe/dt
                                (diff(qx, dims=1)[:, 2:end-1]/dx .+ diff(qy, dims=2)[2:end-1, :]/dy) .+    # div(q)
                                Σ * vo[2:end-1, 2:end-1] .- Γ * vc[2:end-1, 2:end-1]            .-    # dh/dt
                                Λ * m                                                                    # source term
            Res_h          =    (h[2:end-1, 2:end-1] .- h_old[2:end-1, 2:end-1]) / dt  .-
                                Σ * vo[2:end-1, 2:end-1] .+ Γ * vc[2:end-1, 2:end-1]

            # damped rate of change
            dϕ_dτ      = Res_ϕ .+ damp * dϕ_dτ
            dh_dτ        = Res_h .+ damp * dh_dτ

            # update fields
            ϕ[2:end-1, 2:end-1]  = ϕ[2:end-1, 2:end-1] + dτ * dϕ_dτ   # update rule, sets the BC as ϕ[1]=ϕ[end]=0
            h[2:end-1, 2:end-1]    = h[2:end-1, 2:end-1] + dτ * dh_dτ

            err_ϕ = norm(Res_ϕ)/length(Res_ϕ)
            err_h   = norm(Res_h)/length(Res_h)
            iter += 1

        end
        ittot += iter; it += 1; t += dt
        print(iter)

        ϕ_old .= ϕ
        h_old .= h
    end


    return xc, yc, ϕ0, ϕ
end

function plot_output(xc, yc, ϕ0, ϕ)
    pygui(true)
    figure()
    subplot(1, 2, 1)
    pcolor(xc, yc, ϕ0')
    subplot(1, 2, 2)
    pcolor(xc, yc, ϕ')
end



end # module
