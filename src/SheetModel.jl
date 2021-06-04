module SheetModel

using Base: Float64
using PyPlot, Printf, LinearAlgebra, Parameters

export Params

function scaling(p::Params)
    @unpack H, ub, m, lx, lr, g, A, hr, ρw, r_ρ, α, β, n = p

    # Scaling factors
    H_ = mean(H)
    ub_ = ub
    m_ = mean(m)
    x_ = lx
    lr_ = lr
    g_ = g
    A_ = A
    h_ = hr
    ρ_ = ρw
    t_ = x_ * h_ / q_
    ϕ_ = g_ * H_ * ρ_ / r_ρ
    q_ = k_ * h_^α * (ϕ_/x_)^(β-1)
    vo_ = ub_ * h_ / lr_
    vc_ = A_ * h_ * ϕ_^n

    # Dimensionless parameters
    Scaled_Params = Params(p,
    g = g / g_,
    ρw = ρw / ρ_,

    )

end



@with_kw struct Params{F<:Float64, I<:Int64, Arr<:Array{Float64}, Lin<:LinRange{Float64}}
    # Scalars (one global value)
    g::F     = 9.81              # gravitational acceleration, m/s^2
    ρw::F    = 1000.0            # water density, kg/m^3
    ρi::F    = 910.0             # ice density, kg/m^3
    α::F     = 1.25              # first sheet flow exponent
    β::F     = 1.5               # second sheet flow exponent
    k::F     = 0.01              # sheet conductivity, m^(7/4)kg^(-1/2)
    n::F     = 3.0               # Glen's flow law exponent
    A::F     = 5e-25             # ice flow constant, Pa^(-n)s^(-1)
    ev::F    = 1e-3              # englacial void ratio
    lr::F    = 2.0               # horizontal cavity spacing, m
    hr::F    = 0.1               # bedrock bump height, m
    ub::F    = 1e-6              # basal sliding speed, m/s

    # Numerical domain
    lx::F         # domain size
    ly::F
    nx::I         # number of grids
    ny::I
    dx::F = lx/nx      # grid size
    dy::F = ly/ny
    xc::Lin = LinRange(dx/2, lx-dx/2, nx) # vector of x-coordinates
    yc::Lin = LinRange(dy/2, ly-dy/2, ny) # vector of y-coordinates

    # Field parameters (defined on every grid point)
    H::Arr  # ice thickness
    zb::Arr # bed elevation
    m::Arr  # source term

    # Dimensionless numbers
    Σ::F   = nothing
    Γ::F   = nothing
    Λ::F   = nothing
    r_ρ::F = ρw / ρi

    # Physical time stepping
    ttot::F        # total simulation time
    dt::F          # physical time step

    # Pseudo-time iteration
    tol::F    = 1e-6       # tolerance
    itMax::I  = 10^3       # max number of iterations
    damp::F   = 1-41/nx    # damping (this is a tuning parameter, dependent on e.g. grid resolution) # TODO: define how?
    dτ::F     = (1.0/(dx^2/k/2.1) + 1.0/dt)^-1 # pseudo-time step; TODO: define how?
end



# @with_kw struct Para{F<:Float64, I<:Int64, Arr<:Array{Float64}, Lin<:LinRange{Float64}}
#     # Scalars (one global value)
#     g::F     = 1.0              # gravitational acceleration
#     ρw::F    = 1.0              # water density
#     ρi::F    = 910.0 / 1000.0
#     α::F     = 1.25
#     β::F     = 1.5
#     k::F     = 1.0
#     n::F     = 3.0
#     A::F     = 1.0
#     ev::F    = 0.5
#     lr::F    = 1.0              # horizontal cavity spacing, spatially uniform
#     hr::F    = 1.0              # bedrock bump height, spatially uniform
#     # Numerical domain
#     lx::F = 12.0       # domain size
#     ly::F = 10.0
#     nx::I = 64         # number of grids
#     ny::I = 32
#     dx::F = lx/nx      # grid size
#     dy::F = ly/ny
#     xc::Lin = LinRange(dx/2, lx-dx/2, nx) # vector of x-coordinates
#     yc::Lin = LinRange(dy/2, ly-dy/2, ny) # vector of y-coordinates
#     # Field parameters (defined on every grid point)
#     zs::Arr = 1/lx * xc * ones(ny)' # surface elevation
#     zb::Arr = zeros(nx, ny) # bed elevation
#     m::Arr  = zeros(nx-2, ny-2)  # source term
#     ub::Arr = zeros(nx, ny) # basal sliding speed
#     # Dimensionless numbers
#     Σ::F   = 1e2
#     Γ::F   = 1e5
#     Λ::F   = 0.5
#     r_ρ::F = 1000 / 910 # ρw / ρi
#     # Physical time stepping
#     ttot::F   = 2.0        # total simulation time
#     dt::F     = 0.1        # physical time step
#     # Pseudo-time iteration
#     tol::F    = 1e-6       # tolerance
#     itMax::I  = 10^3        # max number of iterations
#     damp::F   = 1-41/nx    # damping (this is a tuning parameter, dependent on e.g. grid resolution) # TODO: define how?
#     dτ::F     = (1.0/(dx^2/k/2.1) + 1.0/dt)^-1 # pseudo-time step; TODO: define how?
# end



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
    @unpack ρi, g, zs, zb = p
    pw = calc_pw(ϕ, p)
    return ρi * g * (zs .- zb) .- pw
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
    return ub .* (hr .- h) ./ lr
end



function runthemodel()


    Params = Para()
    @unpack ev, g, ρw, Σ, Γ, Λ, m, lx, ly, nx, ny, dx, dy, xc, yc, dt, ttot, tol, itMax, damp, dτ = Params

    # Array allocation
    qx, qy, dϕ_dτ, dh_dτ, Res_ϕ, Res_h = array_allocation(Params)

    # Initial condition
    #ϕ0    = exp.(- 1e-2*(xc.-lx/2).^2) * exp.(-1e-2*(yc.-ly/2).^2)'
     ϕ0 = 1/lx * xc * ones(ny)'
    # ϕ0 = rand(nx, ny)
    ϕ_old   = copy(ϕ0)
    ϕ      = copy(ϕ0)

    h0 = 0.5 * ones(nx, ny)
    h_old = copy(h0)
    h = copy(h0)

    t = 0.0; it = 0; ittot = 0

    # Physical time loop
    while t<ttot
        iter = 0; err_ϕ = 2*tol; err_h = 2*tol;

        # Pseudo-transient iteration
        while max(err_ϕ, err_h) >tol && iter<itMax
            dϕ_dx, dϕ_dy = diff(ϕ, dims=1) ./ dx, diff(ϕ, dims=2) ./ dy   # hydraulic gradient
            qx, qy     = calc_q(h[1:end-1, :], dϕ_dx, Params), calc_q(h[:, 1:end-1], dϕ_dy, Params) # TODO: which h values to take?!

            vo     = calc_vo(h, Params)                   # opening rate
            vc     = calc_vc(ϕ, h, Params)  # closure rate

            # calculate residuals
            Res_ϕ        =    - ev/(ρw*g) * (ϕ[2:end-1, 2:end-1] .- ϕ_old[2:end-1, 2:end-1])/dt .-    # dhe/dt
                                (diff(qx, dims=1)[:, 2:end-1]/dx .+ diff(qy, dims=2)[2:end-1, :]/dy)    .+    # div(q)
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

        ϕ_old .= ϕ
        h_old .= h
    end


    return xc, yc, ϕ0, ϕ
end

function plot_output(xc, yc, ϕ0, ϕ)
    figure()
    subplot(1, 2, 1)
    pcolor(xc, yc, ϕ0')
    subplot(1, 2, 2)
    pcolor(xc, yc, ϕ')
end



end # module
