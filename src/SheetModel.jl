__precompile__(false)
module SheetModel

using LinearAlgebra: size
using Base: Float64, Int64
using LinearAlgebra, LazyArrays, Parameters, Statistics, Printf, PyPlot

export Para

"""
All model parameters; physical and numerical
"""
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
    xc::LinRange{Float64} = LinRange(0.0, lx, nx) # vector of x-coordinates
    yc::LinRange{Float64} = LinRange(0.0, ly, ny) # vector of y-coordinates

    # Field parameters (defined on every grid point)
    calc_zs::Function
    calc_zb::Function
    calc_m::Function
    H::Matrix{Float64} = calc_zs.(xc, yc') .- calc_zb.(xc, yc')  # ice thickness, m
    zb::Matrix{Float64} = calc_zb.(xc, yc')  # bed elevation, m
    m::Matrix{Float64} = calc_m.(xc, yc')  # source term, m/s

    # Physical time stepping
    ttot        # total simulation time
    dt          # physical time step

    # misc fudge factors
    small = eps(Float32) # maybe a Float64 is needed here

    # Pseudo-time iteration
    tol    = 1e-6       # tolerance
    itMax  = 10^5       # max number of iterations
    damp   = 0.9        # damping (this is a tuning parameter, dependent on e.g. grid resolution) # TODO: define how?
    #dτ     = (1.0/(min(dx, dy)^2/k/4.1) + 1.0/dt)^-1 # pseudo-time step; TODO: define how?

    # Dimensionless numbers
    Σ   = NaN
    Γ   = NaN
    Λ   = NaN
    r_ρ = ρw / ρi
end
Broadcast.broadcastable(p::Para) = Ref(p)

"""
Convert input parameters to non-dimensional quantities
"""
function scaling(p::Para, ϕ0, h0)
    @unpack g, ρw, k, A, lr, hr,
            H, zb, m, ub,
            lx, ly, dx, dy, xc, yc,
            ttot, dt, #dτ,
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
        #dτ = dτ / t_,

        Σ = vo_ * xy_ / q_,
        Γ = vc_ * xy_ / q_,
        Λ = m_ * xy_ / q_
        )::Para

    # variables
    ϕ0 = ϕ0 ./ ϕ_
    h0 = h0 ./ h_

    return scaled_params, ϕ0, h0, ϕ_, h_
end

"""
Return arrays of initial conditions for ϕ and h
"""
function initial_conditions(xc, yc; calc_ϕ = (x,y) -> 0.0, calc_h = (x,y) -> 0.01)
    ϕ0 = calc_ϕ.(xc, yc')
    h0 = calc_h.(xc, yc')
    return ϕ0, h0
end

"""
Pre-allocate arrays
"""
function array_allocation(nu::Para)
    @unpack nx, ny = nu
    vo     = zeros(nx, ny)
    vc     = zeros(nx, ny)
    dϕ_dx  = zeros(nx-1, ny)
    dϕ_dy  = zeros(nx, ny-1)
    qx     = zeros(nx-1, ny)
    qy     = zeros(nx, ny-1)
    ix     = zeros(nx-1, ny)
    iy     = zeros(nx, ny-1)
    dϕ_dτ  = zeros(nx-2, ny-2)
    dh_dτ  = zeros(nx, ny)
    Res_ϕ  = zeros(nx-2, ny-2)
    Res_h  = zeros(nx, ny)
    Err_ϕ  = zeros(nx, ny)
    d_eff  = zeros(nx, ny)
    dτ     = zeros(nx, ny)
    return vo, vc, dϕ_dx, dϕ_dy, qx, qy, ix, iy, dϕ_dτ, dh_dτ, Res_ϕ, Res_h, Err_ϕ, d_eff, dτ
end

"""
Determine global index for upstream scheme
"""
function upstream(ix, iy, nx, dϕ_du; dims=1)
    if dϕ_du >= 0
        di = 0
    else
        di = 1
    end
    if dims == 1
        ix = ix + di
    elseif dims == 2
        iy = iy + di
    end
    return ix + (iy-1) * nx # global index
end

"""
Calculate discharge
"""
function calc_q(h, dϕ_du, k, α, β, small) # u can be x or y
    #@unpack k, α, β, small = p
    return - k * h^α * (abs(dϕ_du) + small)^(β-2) * dϕ_du
end

"""
Calculate water pressure

# Example
```jldoctest
julia> calc_pw(1.0, 1.0, 1.0, 0.5)
0.5
```
"""
function calc_pw(ϕ, ρw, g, zb)
    #@unpack ρw, g, zb = p
    return ϕ - ρw * g * zb
end

"""
Calculate effective pressure
"""
function calc_N(ϕ, ρi, ρw, g, H, zb)
    #@unpack ρi, g, H = p
    pw = calc_pw(ϕ, ρw, g, zb)
    return ρi * g * H - pw
end

"""
Calculate closure rate
"""
function calc_vc(ϕ, h, ρi, ρw, g, H, zb, n, A)
    #@unpack n, A = p
    N = calc_N(ϕ, ρi, ρw, g, H, zb)
    return 2 / n^n * A * h * abs(N)^(n-1) * N
end

"""
Calculate opening rate
"""
function calc_vo(h, ub, hr, lr)
    #@unpack ub, hr, lr = p
    if h < hr
        vo = ub * (hr - h) / lr
    else
        vo = 0.0
    end
    return vo
end

"""
Scale the parameters and call the model run function
"""
function runthemodel(input::Para, ϕ0, h0;
                    printit=10^5)         # error is printed after printit iterations of pseudo-transient time
    params, ϕ0, h0, ϕ_, h_ = scaling(input, ϕ0, h0)
    ϕ, h = runthemodel_scaled(params::Para, ϕ0, h0, printit)
    ϕ .= ϕ .* ϕ_
    h .= h .* h_
    return ϕ, h
end

"""
Run the model with scaled parameters
"""
function runthemodel_scaled(params::Para, ϕ0, h0, printit)
    @unpack ev, g, ρw, ρi, n, A, Σ, Γ, Λ, m, dx, dy, nx, ny, k, α, β, small,
            H, zb, ub, hr, lr, dt, ttot, tol, itMax, damp = params
    # Array allocation
    vo, vc, dϕ_dx, dϕ_dy, qx, qy, ix, iy,
    dϕ_dτ, dh_dτ, Res_ϕ, Res_h, Err_ϕ, d_eff, dτ = array_allocation(params)

    ϕ0[1, :]   = ρw .* g .* zb[1, :] # boundary condition, zero water pressure
    ϕ0[end, :] = ϕ0[end-1, :] # no flux boundary condition
    ϕ0[:, 1]   = ϕ0[:, 2]
    ϕ0[:, end] = ϕ0[:, end-1]
    ϕ_old   = copy(ϕ0)
    ϕ       = copy(ϕ0)
    h_old   = copy(h0)
    h       = copy(h0)

    t, it, ittot = 0.0, 0, 0

    # Physical time loop
    while t<ttot
        iter, err_ϕ, err_h = 0, 2*tol, 2*tol

        # Pseudo-transient iteration
        while max(err_ϕ, err_h)>tol && iter<itMax
            # save current ϕ for error calculation
            Err_ϕ   .= ϕ

            # quantities occurring in the equations
            dϕ_dx   .= LazyArrays.Diff(ϕ, dims=1) ./ dx                  # hydraulic gradient
            dϕ_dy   .= LazyArrays.Diff(ϕ, dims=2) ./ dy

            ix = upstream.(collect(1:nx-1), collect(1:ny)', nx, dϕ_dx; dims=1) # determine indexes of h that are upstream of dϕ/dx
            iy = upstream.(collect(1:nx), collect(1:ny-1)', nx, dϕ_dy; dims=2)
            qx      .= calc_q.(h[ix], dϕ_dx, k, α, β, small)
            qy      .= calc_q.(h[iy], dϕ_dy, k, α, β, small)

            # test whether boundary conditions are fulfilled
            if any(qx[end, :] .!= 0.0) || any(ϕ[1, :] .!= ρw .* g .* zb[1, :])
                @warn "qx b.c. not fulfilled"
            end
            if any(qy[:, 1] .!= 0.0) || any(qy[:, end] .!= 0.0)
                @warn "qy b.c. not fulfilled"
            end

            vo     .= calc_vo.(h, ub, hr, lr)                 # opening rate
            vc     .= calc_vc.(ϕ, h, ρi, ρw, g, H, zb, n, A)  # closure rate

            # calculate residuals
            Res_ϕ       .=      - ev/(ρw*g) * (ϕ[2:end-1, 2:end-1] .- ϕ_old[2:end-1, 2:end-1])/dt .-         # dhe/dt
                                (LazyArrays.Diff(qx, dims=1)[:, 2:end-1]/dx .+ LazyArrays.Diff(qy, dims=2)[2:end-1, :]/dy) .-      # div(q)
                                (Σ * vo[2:end-1, 2:end-1] .- Γ * vc[2:end-1, 2:end-1])            .+         # dh/dt
                                Λ * m[2:end-1, 2:end-1]                                                      # source term
            Res_h       .=      - (h .- h_old) / dt  .+
                                (Σ * vo .- Γ * vc)

            # determine pseudo-time step
            d_eff .= k * h.^α # effective diffusivity
            dτ    .= (1.0 ./ (min(dx, dy)^2 ./ d_eff / 4.1) .+ 1.0 / dt) .^(-1) # pseudo-time step; TODO: define how?

            # damped rate of change
            dϕ_dτ      .= Res_ϕ .+ damp .* dϕ_dτ
            dh_dτ      .= Res_h .+ damp .* dh_dτ

            # update fields
            ϕ[2:end-1, 2:end-1]  .= ϕ[2:end-1, 2:end-1] .+ dτ[2:end-1, 2:end-1] .* dϕ_dτ   # update ϕ
            h  .= h .+ dτ .* dh_dτ   # update h

            # boundary conditions (no flux)
            ϕ[end, :] .= ϕ[end-1, :]
            ϕ[:, 1]   .= ϕ[:, 2]
            ϕ[:, end] .= ϕ[:, end-1]

            Err_ϕ .= abs.(Err_ϕ .- ϕ)
            err_ϕ = norm(Err_ϕ) # this error is smaller than the error using Res_ϕ
            # err_ϕ = norm(Res_ϕ)/length(Res_ϕ) # but with this error it also converges
            err_h   = norm(Res_h)/length(Res_h)
            iter += 1

            if mod(iter, printit) == 0
                @printf("iterations = %d, error ϕ = %1.2e, error h = %1.2e \n", iter, err_ϕ, err_h)
            end
        end
        ittot += iter; it += 1; t += dt
        @printf("time = %d, number of iterations = %d \n", it, iter)

        ϕ_old .= ϕ
        h_old .= h
    end

    # give the effective pressure as output instead of the hydraulic potential
    N = calc_N.(ϕ, ρi, ρw, g, H, zb)

    return N, h
end

function plot_output(xc, yc, N, h)
    x_plt = [0; xc .+ (xc[2]-xc[1])]
    y_plt = [0; yc .+ (yc[2]-yc[1])]
    pygui(true)
    figure()
    subplot(1, 2, 1)
    pcolor(x_plt, y_plt, h', edgecolors="black")
    colorbar()
    title("h")
    subplot(1, 2, 2)
    pcolor(x_plt, y_plt, N', edgecolors="black")
    colorbar()
    title("N")
end


end # module
