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
    yc::LinRange{Float64} = LinRange(-ly/2, ly/2, ny) # vector of y-coordinates

    # Field parameters (defined on every grid point)
    calc_zs::Function
    calc_zb::Function
    calc_m_xyt::Function # fct(x, y, t)
    fct_pos::Function = x -> x > 0.0 ? x : 0.0 # turn all negative numbers into 0.0
    H::Matrix{Float64} = fct_pos.(calc_zs.(xc, yc') .- calc_zb.(xc, yc'))  # ice thickness, m
    zb::Matrix{Float64} = calc_zb.(xc, yc')                      # bed elevation, m
    calc_m_t::Function = t -> calc_m_xyt.(xc, yc', t)                # source term, m/s, fct(t)

    # Physical time stepping
    ttot        # total simulation time
    dt          # physical time step

    # misc fudge factors
    small = eps(Float32) # maybe a Float64 is needed here

    # Pseudo-time iteration
    tol    = 1e-6       # tolerance
    itMax  = 10^3       # max number of iterations
    γ_ϕ    = 1e-3        # damping parameter for ϕ update
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
Convert input parameters to non-dimensional quantities
"""
function scaling(p::Para, ϕ0, h0)
    @unpack g, ρw, k, A, lr, hr,
            H, zb, calc_m_t, ub,
            lx, ly, dx, dy, xc, yc,
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
    xy_ = max(lx, ly)

    H_ = mean(H)
    zb_ = H_
    m_ = mean(calc_m_t(0.0)) # for time-dependent input: temporal peak
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
        calc_m_t = t -> calc_m_t(t) ./ m_,

        # Numerical domain
        lx = lx ./ xy_,
        ly = ly ./ xy_,
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

    return scaled_params, ϕ0, h0, ϕ_, h_
end

"""
Return arrays of initial conditions for ϕ and h
"""
function initial_conditions(xc, yc, H; calc_ϕ = (x,y) -> 0.0, calc_h = (x,y) -> 0.01)
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
    m      = zeros(nx, ny)
    dϕ_dτ  = zeros(nx-2, ny-2)
    dh_dτ  = zeros(nx, ny)
    Res_ϕ  = zeros(nx-2, ny-2)
    Res_h  = zeros(nx, ny)
    Err_ϕ  = zeros(nx, ny)
    Err_h  = zeros(nx, ny)
    d_eff  = zeros(nx, ny)
    dτ_ϕ   = zeros(nx, ny)
    return vo, vc, dϕ_dx, dϕ_dy, qx, qy, ix, iy, m, dϕ_dτ, dh_dτ, Res_ϕ, Res_h, Err_ϕ, Err_h, d_eff, dτ_ϕ
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
Apply the boundary conditions to ϕ, at the moment pw(x=0) = 0 and no flux at the other boundaries
"""
function apply_bc(ϕ, h, H, ρw, g, zb) # TODO: shouldn't have any function with arrays as arguments
                                   # TODO: don't hard-wire, give bc as input parameters
    nx, ny = size(ϕ)
    ϕ[H .== 0.0] .= ρw .* g .* zb[H .== 0.0]  # zero water pressure outside of glacier domain
    h[H .== 0.0] .= 0.0                       # zero sheet thickness outside of glacier domain
    for i = 2:nx-1, j = 2:ny-1
        if H[i, j] > 0.0
            if H[i-1, j] == 0.0  || i==2 # x1 boundary
                ϕ[i-1, j] = ϕ[i, j] # no flux
            end
            if H[i+1, j] == 0.0 || i==nx-1 # xend boundary
                ϕ[i+1, j] = ϕ[i, j] # no flux
            end
            if H[i, j-1] == 0.0 || j==2 # y1 boundary
                ϕ[i, j-1] = ϕ[i, j] # no flux
            end
            if H[i, j+1] == 0.0 || j==ny-1 # yend boundary
                ϕ[i, j+1] = ϕ[i, j] # no flux
            end
            if H[i, j] > 0.0 && i==2 # x1 boundary
                ϕ[i-1, j] = ρw .* g .* zb[i-1, j] # zero water pressure
            end
        end
    end

    # corner points
    ϕ[nx, 1] = H[nx, 2] > 0.0 ? ϕ[nx, 2] : ϕ[nx, 1]
    ϕ[nx, 1] = H[nx-1, 1] > 0.0 ? ϕ[nx-1, 1] : ϕ[nx, 1]
    ϕ[nx, ny] = H[nx, ny-1] > 0.0 ? ϕ[nx, ny-1] : ϕ[nx, ny]
    ϕ[nx, ny] = H[nx-1, ny] > 0.0 ? ϕ[nx-1, ny] : ϕ[nx, ny]
    ϕ[1, 1]  = H[2, 1] > 0.0 ? ρw .* g .* zb[1, 1] : ϕ[1, 1]
    ϕ[1, ny] = (H[1, ny-1] > 0.0 || H[2, ny] > 0.0) ?  ρw .* g .* zb[1, ny] : ϕ[1, ny]

    return ϕ, h
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
                    printit=10^5,         # error is printed after `printit` iterations of pseudo-transient time
                    printtime=10^5)       # time step and number of PT iterations is printed after `printtime` number of physical time steps
    params, ϕ0, h0, ϕ_, h_ = scaling(input, ϕ0, h0)
    N, ϕ, h, nit = runthemodel_scaled(params::Para, ϕ0, h0, printit, printtime)
    ϕ .= ϕ .* ϕ_
    h .= h .* h_
    return N, ϕ, h, nit
end

"""
Run the model with scaled parameters
"""
function runthemodel_scaled(params::Para, ϕ0, h0, printit, printtime)
    @unpack ev, g, ρw, ρi, n, A, Σ, Γ, Λ, calc_m_t, dx, dy, nx, ny, k, α, β, small,
            H, zb, ub, hr, lr, dt, ttot, tol, itMax, γ_ϕ, γ_h, dτ_ϕ_, dτ_h_ = params
    # Array allocation
    vo, vc, dϕ_dx, dϕ_dy, qx, qy, ix, iy, m,
    dϕ_dτ, dh_dτ, Res_ϕ, Res_h, Err_ϕ, Err_h, d_eff, dτ_ϕ = array_allocation(params)

    # Apply boundary conditions
    ϕ0, h0 = apply_bc(ϕ0, h0, H, ρw, g, zb)

    ϕ_old   = copy(ϕ0)
    ϕ       = copy(ϕ0)
    h_old   = copy(h0)
    h       = copy(h0)

    # determine indices of glacier domain
    idx_ice = H .> 0.0

    # initiate time loop parameters
    t, it, ittot = 0.0, 0, 0

    # Physical time loop
    while t<ttot
        iter, err_ϕ, err_h = 0, 2*tol, 2*tol

        m .= calc_m_t(t)
        # Pseudo-transient iteration
        while !(max(err_ϕ, err_h) < tol) && iter<itMax # with the ! the loop also continues for NaN values of err
            # save current ϕ for error calculation
            Err_ϕ   .= ϕ
            Err_h   .= h

            # quantities occurring in the equations
            dϕ_dx   .= LazyArrays.Diff(ϕ, dims=1) ./ dx                  # hydraulic gradient
            dϕ_dy   .= LazyArrays.Diff(ϕ, dims=2) ./ dy

            ix = upstream.(collect(1:nx-1), collect(1:ny)', nx, dϕ_dx; dims=1) # determine indexes of h that are upstream of dϕ/dx
            iy = upstream.(collect(1:nx), collect(1:ny-1)', nx, dϕ_dy; dims=2)
            qx      .= calc_q.(h[ix], dϕ_dx, k, α, β, small)
            qy      .= calc_q.(h[iy], dϕ_dy, k, α, β, small)

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
            d_eff .= k * h.^α                                                                 # effective diffusivity, defined on each grid point
            dτ_ϕ  .= (1.0/dτ_ϕ_) .* (1.0 ./ (min(dx, dy)^2 ./ d_eff / 4.1) .+ 1.0 / dt) .^(-1) # pseudo-time step for ϕ, defined on each grid point
            dτ_h   = dt / dτ_h_                                                                # pseudo-time step for h, scalar

            # damped rate of change
            dϕ_dτ      .= Res_ϕ .+ γ_ϕ .* dϕ_dτ
            dh_dτ      .= Res_h .+ γ_h .* dh_dτ

            # update fields
            ϕ[2:end-1, 2:end-1]  .= ϕ[2:end-1, 2:end-1] .+ dτ_ϕ[2:end-1, 2:end-1] .* dϕ_dτ   # update ϕ (only interior points because fluxes only defined there)
            h                    .= h .+ dτ_h .* dh_dτ                                      # update h

            # apply boundary conditions
            ϕ, h = apply_bc(ϕ, h, H, ρw, g, zb)

            # determine the errors (only consider points where the ice thickness is > 0)
            Err_ϕ .= abs.(Err_ϕ .- ϕ)
            err_ϕ = norm(Err_ϕ[idx_ice]) # this error is smaller than the error using Res_ϕ
            # err_ϕ = norm(Res_ϕ)/length(Res_ϕ) # with this error it does not converge
            Err_h .= abs.(Err_h .- h)
            err_h = norm(Err_h[idx_ice])
            # err_h   = norm(Res_h)/length(Res_h)
            iter += 1

            if mod(iter, printit) == 0
                @printf("iterations = %d, error ϕ = %1.2e, error h = %1.2e \n", iter, err_ϕ, err_h)
            end
        end
        ittot += iter; it += 1; t += dt
        if mod(it, printtime) == 0
            @printf("time step = %d, number of iterations = %d \n", it, iter)
        end

        ϕ_old .= ϕ
        h_old .= h
    end

    # give the effective pressure as output instead of the hydraulic potential
    N = calc_N.(ϕ, ρi, ρw, g, H, zb)

    return N, ϕ, h, ittot
end

function plot_output(xc, yc, N, h)
    x_plt = [0; xc .+ (xc[2]-xc[1])]
    y_plt = [0; yc .+ (yc[2]-yc[1])]
    pygui(true)
    # pcolor of ϕ and h fields
    figure()
    subplot(2, 2, 1)
    pcolor(x_plt, y_plt, h', edgecolors="black")
    colorbar()
    title("h")
    subplot(2, 2, 2)
    pcolor(x_plt, y_plt, N', edgecolors="black")
    colorbar()
    title("N")
    # cross-sections of ϕ and h
    subplot(2, 2, 3)
    plot(xc, h[:, Int(length(yc)/2)])
    title(join(["h at y = ", string(round(yc[10], digits=1))]))
    subplot(2, 2, 4)
    plot(xc, N[:, Int(length(yc)/2)])
    title(join(["h at y = ", string(round(yc[10], digits=1))]))
end


end # module
