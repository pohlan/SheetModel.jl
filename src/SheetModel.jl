__precompile__(false)
module SheetModel
using Infiltrator

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
    xrange::Tuple{Float64, Float64}   # domain size
    yrange::Tuple{Float64, Float64}
    x1   = xrange[1]
    xend = xrange[2]
    y1   = yrange[1]
    yend = yrange[2]
    nx::Int64                        # number of grid points, including ghost points where ice thickness = 0
    ny::Int64
    dx = (xend-x1) / (nx-3)          # grid size
    dy = (yend-y1) / (ny-3)
    xc::LinRange{Float64} = LinRange(x1-dx, xend+dx, nx) # vector of x-coordinates
    yc::LinRange{Float64} = LinRange(y1-dy, yend+dy, ny) # vector of y-coordinates

    # Field parameters (defined on every grid point)
    calc_zs::Function
    calc_zb::Function
    calc_m_xyt::Function # fct(x, y, t)
    fct_pos::Function = x -> x > 0.0 ? x : 0.0 # turn all negative numbers into 0.0
    gp_x::Array{Float64, 1} = [0.0; ones(length(xc)-2); 0.0] # to achieve H = 0 at ghost points
    gp_y::Array{Float64, 1} = [0.0; ones(length(yc)-2); 0.0]
    H::Matrix{Float64} = (gp_x * gp_y') .* ( fct_pos.(calc_zs.(xc, yc') .- calc_zb.(xc, yc')) )  # ice thickness, m
    zb::Matrix{Float64} = calc_zb.(xc, yc')                      # bed elevation, m
    calc_m_t::Function = t -> calc_m_xyt.(xc, yc', t)                # source term, m/s, fct(t)

    # Physical time stepping
    ttot        # total simulation time
    dt          # physical time step

    # misc fudge factors
    small = eps(Float32) # maybe a Float64 is needed here

    # Pseudo-time iteration
    tol    = 1e-6       # tolerance
    itMax  = 5*10^3       # max number of iterations
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
    qx     = zeros(nx-1, ny-2)
    qy     = zeros(nx-2, ny-1)
    gp_ice  = zeros(Int, nx+2, ny+2)
    ix     = zeros(Int, nx-1, ny)
    iy     = zeros(Int, nx, ny-1)
    m      = zeros(nx, ny)
    div_q  = zeros(nx, ny)
    div_ϕ  = zeros(nx, ny)
    dϕ_dτ  = zeros(nx, ny)
    dh_dτ  = zeros(nx, ny)
    Res_ϕ  = zeros(nx, ny)
    Res_h  = zeros(nx, ny)
    Err_ϕ  = zeros(nx, ny)
    Err_h  = zeros(nx, ny)
    d_eff  = zeros(nx-2, ny-2)
    dτ_ϕ   = zeros(nx, ny)
    return vo, vc, dϕ_dx, dϕ_dy, qx, qy, gp_ice, ix, iy, m, div_q, div_ϕ, dϕ_dτ, dh_dτ, Res_ϕ, Res_h, Err_ϕ, Err_h, d_eff, dτ_ϕ
end

"""
Determine global index for upstream scheme
"""
function upstream(ix, iy, nx, dϕ_du; dims=1)
    di = dϕ_du >= 0 ? 0 : 1
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
    #ϕ[H .== 0.0] .= ρw .* g .* zb[H .== 0.0] # zero water pressure outside of glacier domain
    h[H .== 0.0] .= 0.0                       # zero sheet thickness outside of glacier domain
    ϕ[2, :] .= ρw .* g .* zb[2, :]
    #ϕ[end, :] .= ϕ[end-1, :]
    #ϕ[:, 1] .= ϕ[:, 2]
    #ϕ[:, end] .= ϕ[:, end-1]

    # ϕ[end-1,:] = ρw .* g .* zb[end-1,:]
    # ϕ[1, ny÷2+1] = ρw .* g .* zb[1, ny÷2+1]
    # if iseven(size(ϕ,2))
    #     ϕ[1, ny÷2] = ρw .* g .* zb[1, ny÷2]
    # end

    # for j = 2:ny-1, i = 2:nx-1
    #     if H[i, j] > 0.0
    #         if H[i-1, j] == 0.0 # x1 boundary
    #             ϕ[i, j] = ρw .* g .* zb[i-1, j] # zero water pressure
    #         end
    #         if H[i-1, j] > 0.0 && i == 2  # ice boundary = domain boundary
    #             ϕ[i-1, j] = ρw .* g .* zb[i-1, j] # zero water pressure
    #             #ϕ[i, j] = ρw .* g .* zb[i, j]
    #         end
    #     end
    # end

    # corner points
    # ϕ[1, 1]  = H[1, 1] > 0.0 ? ρw .* g .* zb[1, 1] : ϕ[1, 1]
    # ϕ[1, ny] = H[1, ny] > 0.0 ?  ρw .* g .* zb[1, ny] : ϕ[1, ny]

    return ϕ, h
end

"""
Calculate discharge
"""
function calc_q(h, dϕ_du1, dϕ_du2, k, α, β, small) # qx -> u1 = x, u2 = y; other way round for qy
    #@unpack k, α, β, small = p
    return - k * h^α * (sqrt(dϕ_du1^2 + dϕ_du2^2) + small)^(β-2) * dϕ_du1
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

@views av(A)    = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end]) # average
@views av_xa(A) = 0.5.*(A[1:end-1,:].+A[2:end,:]) # average x-dir
@views av_ya(A) = 0.5.*(A[:,1:end-1].+A[:,2:end]) # average y-dir

"""
Scale the parameters and call the model run function
"""
function runthemodel(input::Para, ϕ0, h0;
                    printit=10^5,         # error is printed after `printit` iterations of pseudo-transient time
                    printtime=10^5)       # time step and number of PT iterations is printed after `printtime` number of physical time steps
    params, ϕ0, h0, ϕ_, N_, h_, q_ = scaling(input, ϕ0, h0)
    N, ϕ, h, qx, qy, nit, err_ϕ, err_h, qx_ice, qy_ice = runthemodel_scaled(params::Para, ϕ0, h0, printit, printtime)
    N .= N .* N_ # scaling for N same as for ϕ
    ϕ .= ϕ .* ϕ_
    h .= h .* h_
    qx .= qx .* q_
    qy .= qy .* q_
    err_ϕ .= err_ϕ .* ϕ_
    err_h .= err_h .* h_
    return N, ϕ, h, qx, qy, nit, err_ϕ, err_h, qx_ice, qy_ice
end

"""
Run the model with scaled parameters
"""
function runthemodel_scaled(params::Para, ϕ0, h0, printit, printtime)
    @unpack ev, g, ρw, ρi, n, A, Σ, Γ, Λ, calc_m_t, dx, dy, nx, ny, k, α, β, small,
            H, zb, ub, hr, lr, dt, ttot, tol, itMax, γ_ϕ, γ_h, dτ_ϕ_, dτ_h_ = params
    # Array allocation
    vo, vc, dϕ_dx, dϕ_dy, qx, qy, gp_ice, ix, iy, m, div_q, div_ϕ,
    dϕ_dτ, dh_dτ, Res_ϕ, Res_h, Err_ϕ, Err_h, d_eff, dτ_ϕ = array_allocation(params)

    # determine indices of glacier domain
    idx_ice = H .> 0.0
    #qx_ice = zeros(Int, nx-1, ny)
    #qx_ice  .= idx_ice[1:end-1, :]
    #qx_ice .+= idx_ice[2:end, :]   # qx_ice = 1 for boundary, 2 for interior
    #qy_ice = zeros(Int, nx, ny-1)
    #qy_ice  .= idx_ice[:, 1:end-1]
    #qy_ice .+= idx_ice[:, 2:end]
    #qx_xlbound  = Array(LazyArrays.Diff(idx_ice, dims=1) .== 1) # on qx grid
    #qx_xubound  = Array(LazyArrays.Diff(idx_ice, dims=1) .== -1)
    qx_ice = zeros(Int, nx-1, ny-2)
    qx_ice  .= idx_ice[1:end-1, 2:end-1]
    qx_ice .+= idx_ice[2:end, 2:end-1]   # qx_ice = 1 for boundary, 2 for interior
    qy_ice = zeros(Int, nx-2, ny-1)
    qy_ice  .= idx_ice[2:end-1, 1:end-1]
    qy_ice .+= idx_ice[2:end-1, 2:end]
    qx_xlbound  = Array(LazyArrays.Diff(idx_ice[:, 2:end-1], dims=1) .== 1) # on qx grid
    qx_xubound  = Array(LazyArrays.Diff(idx_ice[:, 2:end-1], dims=1) .== -1)
    qy_ybound  = qy_ice .== 1 # on qy grid

    ϕ_outbound =    [LazyArrays.Diff(idx_ice, dims=1) .== 1; zeros(1, ny)] .+ # cells on ϕ grid that are just outside of H > 0 (idx_ice) area
                    [LazyArrays.Diff(idx_ice, dims=2) .== 1 zeros(nx, 1)] .+
                    [zeros(1, ny); LazyArrays.Diff(idx_ice, dims=1) .== -1] .+
                    [zeros(nx, 1) LazyArrays.Diff(idx_ice, dims=2) .== -1]

    # Apply boundary conditions
    ϕ0, h0 = apply_bc(ϕ0, h0, H, ρw, g, zb)

    ϕ_old   = copy(ϕ0)
    ϕ       = copy(ϕ0)
    h_old   = copy(h0)
    h       = copy(h0)

    # initiate time loop parameters
    t, it, ittot = 0.0, 0, 0

    # Physical time loop
    while t<ttot
        iter, err_ϕ, err_h = 0, 2*tol, 2*tol

        m .= calc_m_t(t+dt)
        # Pseudo-transient iteration
        while !(max(err_ϕ, err_h) < tol) && iter<itMax # with the ! the loop also continues for NaN values of err
        #while !(err_ϕ < tol) && iter<itMax # with the ! the loop also continues for NaN values of err
            # used indices:
            # - normal grid (i,j)
            #   e.g. ϕ[i,j]
            # - staggered grid (m,n)
            #   e.g. dϕ_dx[m,j]
            # - staggered grid with ghost points (p,q)
            #   e.g. qx[p,j], qy[i,q]

            # save current ϕ for error calculation
            Err_ϕ   .= ϕ
            Err_h   .= h

            # quantities occurring in the equations
            dϕ_dx   .= LazyArrays.Diff(ϕ, dims=1) ./ dx                  # hydraulic gradient
            dϕ_dy   .= LazyArrays.Diff(ϕ, dims=2) ./ dy

            #gradϕ = sqrt.(av_xa(dϕ_dx[:, 2:end-1]).^2 .+ av_ya(dϕ_dy[2:end-1, :]).^2) # on ϕ/h grid, size (nx-2, ny-2)
            #d_eff .= k*h[2:end-1, 2:end-1].^α .* (gradϕ .+ small).^(β-2) # on ϕ/h grid, size (nx-2, ny-2)

            # upstream
            #qx[2:end-1, 2:end-1] .= (dϕ_dx[2:end-1, 2:end-1] .>= 0.0) .* (.- d_eff[2:end, :]   .* dϕ_dx[2:end-1, 2:end-1]) .+
            #                        (dϕ_dx[2:end-1, 2:end-1] .< 0.0)  .* (.- d_eff[1:end-1, :] .* dϕ_dx[2:end-1, 2:end-1])
            #qy[2:end-1, 2:end-1] .= (dϕ_dy[2:end-1, 2:end-1] .>= 0.0) .* (.- d_eff[:, 2:end]   .* dϕ_dy[2:end-1, 2:end-1]) .+
            #                        (dϕ_dy[2:end-1, 2:end-1] .< 0.0)  .* (.- d_eff[:, 1:end-1] .* dϕ_dy[2:end-1, 2:end-1])
            # central differences
            #qx[2:end-1, 2:end-1] .= 0.5 .* (.- d_eff[2:end, :]   .* dϕ_dx[2:end-1, 2:end-1]) .+
            #                        0.5 .* (.- d_eff[1:end-1, :] .* dϕ_dx[2:end-1, 2:end-1])
            #qy[2:end-1, 2:end-1] .= 0.5 .* (.- d_eff[:, 2:end]   .* dϕ_dy[2:end-1, 2:end-1]) .+
            #                        0.5 .* (.- d_eff[:, 1:end-1] .* dϕ_dy[2:end-1, 2:end-1])

            # same as in Ludovic's ice flow code (gradϕ and d_eff defined on cells centered between ϕ/h points)
            # https://github.com/luraess/julia-parallel-course-EGU21/blob/b700a9ad0d1d14f26e0590ce41a8e63380a1df15/scripts/iceflow.jl#L66-L72
            gradϕ = sqrt.(av_ya(dϕ_dx).^2 .+ av_xa(dϕ_dy).^2) # grid centered between ϕ/h points, size (nx-1, ny-1)
            d_eff = k*av(h).^α .* (gradϕ .+ small).^(β-2) # same grid as gradϕ
            qx   .= .-av_ya(d_eff).*diff(ϕ[:,2:end-1], dims=1)/dx # qx grid, size (nx-1, ny-2)
            qy   .= .-av_xa(d_eff).*diff(ϕ[2:end-1,:], dims=2)/dy # qy grid, size (nx-2, nx-1)

            # # determine indexes of h that are upstream of dϕ/dx
            # ix = upstream.(1:nx-1, (1:ny)', nx, dϕ_dx; dims=1)
            # iy = upstream.(1:nx, (1:ny-1)', nx, dϕ_dy; dims=2)
            # qx[2:end-1, :]         .= calc_q.(h[ix], dϕ_dx, k, α, β, small)
            # qy[:, 2:end-1]         .= calc_q.(h[iy], dϕ_dy, k, α, β, small)

            # dϕ_dy_ = 0.0
            # dϕ_dx_ = 0.0
            # for j = 1:size(qx,2)
            #     for p = (1, size(qx,1))
            #         # outer boundary
            #         @assert qx[p,j] == 0
            #     end
            #     for p = 2:size(qx,1)-1
            #         if idx_ice[p,j] == 0 || idx_ice[p-1,j] == 0
            #             # BC: zero flux across boundary
            #             @assert qx[p,j] == 0 # (already set at initialization)
            #         else
            #             dϕ_dx_ = (ϕ[p,j] - ϕ[p-1,j])/dx # i.e. between [i,j] and [i-1,j]
            #             @assert (dϕ_dx_ == dϕ_dx[p-1,j]) (p,j)
            #             i = dϕ_dx_ >= 0 ?
            #                 p : # flux from cell [p,j] to [p-1,j]
            #                 p-1 # flux from cell [p-1,j] to [p,j]\
            #             if j-1 >= 1 && j <= size(dϕ_dy, 2)
            #                 dϕ_dy_ = 0.5 * (dϕ_dy[i, j-1] + dϕ_dy[i, j]) # take average of upstream dϕ_dy
            #             elseif j-1 < 1
            #                 dϕ_dy_ = dϕ_dy[i, j] # lower y-boundary, only one upstream value available
            #             elseif j > size(dϕ_dy, 2)
            #                 dϕ_dy_ = dϕ_dy[i, j-1] # upper y-boundary, ...
            #             end
            #             qx[p,j] = calc_q(h[i,j], dϕ_dx_, dϕ_dy_, k, α, β, small)                  # upstream scheme
            #             # qx[p,j] = calc_q((h[p ,j] + h[p-1 ,j])/2, dϕ_dx_, k, α, β, small) # central differences
            #         end
            #     end
            # end
            # for q=1:size(qy,2)
            #     if q == 1 || q == size(qy,2)
            #         # outer boundary
            #         @assert all(qy[:,q] .== 0)
            #         continue
            #     end
            #     for i=1:size(qy,1)
            #         if idx_ice[i,q] == 0 || idx_ice[i,q-1] == 0
            #             # BC: zero flux across boundary
            #             @assert qy[i,q] == 0 # (already set at initialization)
            #         else
            #             dϕ_dy_ = (ϕ[i,q] - ϕ[i,q-1])/dy # i.e. between [i,j] and [i,j-1]
            #             @assert (dϕ_dy_ == dϕ_dy[i,q-1]) (i,q)
            #             j = dϕ_dy_ >= 0 ?
            #                 q : # flux from cell [i,q] to [i,q-1]
            #                 q-1 # flux from cell [i,q-1] to [i,q]
            #             # @assert qy[i,q] == calc_q(h[i,j], dϕ_dy_, k, α, β, small) i,q
            #             if i-1 >= 1 && i <= size(dϕ_dx, 1)
            #                 dϕ_dx_ = 0.5 * (dϕ_dx[i-1, j] + dϕ_dx[i, j])
            #             elseif i-1 < 1
            #                 dϕ_dx_ = dϕ_dx[i, j]
            #             elseif i > size(dϕ_dx, 1)
            #                 dϕ_dx_ = dϕ_dx[i-1, j]
            #             end
            #             qy[i,q] = calc_q(h[i,j], dϕ_dy_, dϕ_dx_, k, α, β, small)                  # upstream scheme
            #             # qy[i,q] = calc_q((h[i,q] + h[i,q-1])/2, dϕ_dy_, k, α, β, small) # central differences
            #         end
            #     end
            # end
       #     @infiltrate iter==4000

            # set flux boundary conditions
            qx[qx_ice .== 0]   .= 0.0
            qy[qy_ice .== 0]   .= 0.0
            qx[qx_xubound] .= 0.0 # no flux boundary condition
            qx[qx_xlbound] .= 0.0
            qy[qy_ybound] .= 0.0

            vo     .= calc_vo.(h, ub, hr, lr)                 # opening rate
            vc     .= calc_vc.(ϕ, h, ρi, ρw, g, H, zb, n, A)  # closure rate
            #div_q[2:end-1, 2:end-1]  .= LazyArrays.Diff(qx, dims=1)[:, 2:end-1]/dx .+ LazyArrays.Diff(qy, dims=2)[2:end-1, :]/dy .+ small
            div_q[2:end-1, 2:end-1]  .= LazyArrays.Diff(qx, dims=1)/dx .+ LazyArrays.Diff(qy, dims=2)/dy .+ small

            # calculate residuals
            Res_ϕ   .=  idx_ice .* (
                            - ev/(ρw*g) * (ϕ .- ϕ_old)/dt .-         # dhe/dt
                            div_q .-                                 # div(q)
                            (Σ * vo .- Γ * vc)            .+         # dh/dt
                            Λ * m                                    # source term
                            )
            Res_h   .=  idx_ice .* (
                            - (h .- h_old) / dt  .+
                            (Σ * vo .- Γ * vc)
                            )

            # determine effecive diffusivity for pseudo-time step of ϕ: d_eff = divergence(q) / divergence(ϕ)
            #div_ϕ[2:end-1, 2:end-1] .= LazyArrays.Diff(dϕ_dx, dims=1)[:, 2:end-1] / dx .+ LazyArrays.Diff(dϕ_dy, dims=2)[2:end-1, :] / dy .+ small
            # dϕ_dy, dϕ_dx not defined at y-, x-boundary and corner points
            # could be calculated from mean between qy, qx fluxes but maybe ok like this
            # since mostly no flux boundaries (and for dirichlet boundaries it's not relevant)
            #div_ϕ[2:end-1, [1, end]] .= LazyArrays.Diff(dϕ_dx, dims=1)[:, [1, end]] / dx .+ small
            #div_ϕ[[1, end], 2:end-1] .= LazyArrays.Diff(dϕ_dy, dims=2)[[1, end], :] / dy .+ small
            # div_ϕ[[1, 1, end, end], [1, end, 1, end]] .= small # corner points

            #d_eff .=  abs.(div_q ./ div_ϕ)                          # effective diffusivity, spatially variable
            #d_eff[end, 1] = 0.5 * (d_eff[end-1, 1] + d_eff[end, 2]) # for corner points take average of neighbours, otherwise they are far too large
            #d_eff[end, end] = 0.5 * (d_eff[end-1, end] + d_eff[end, end-1])

            #dτ_ϕ[2:end-1, 2:end-1] .= (1.0/dτ_ϕ_) .* (1.0 ./ (min(dx, dy)^2 ./ d_eff / 4.1) .+ 1.0 / dt) .^(-1)
            dτ_ϕ[2:end-1, 2:end-1] .= (1.0/dτ_ϕ_) .* (1.0 ./ (min(dx, dy)^2 ./ av(d_eff) / 4.1) .+ 1.0 / dt) .^(-1)
            dτ_h   = dt / dτ_h_   # pseudo-time step for h, scalar

            # damped rate of change
            dϕ_dτ .= Res_ϕ .+ γ_ϕ .* dϕ_dτ
            dh_dτ .= Res_h .+ γ_h .* dh_dτ

            # update fields
            ϕ .= ϕ .+ dτ_ϕ .* dϕ_dτ   # update ϕ
            h .= h .+ dτ_h .* dh_dτ   # update h

            # apply boundary conditions
            ϕ, h = apply_bc(ϕ, h, H, ρw, g, zb)

            # determine the errors (only consider points where the ice thickness is > 0)
            Err_ϕ .= abs.(Err_ϕ .- ϕ)
            err_ϕ = norm(Err_ϕ) # this error is smaller than the error using Res_ϕ
            # err_ϕ = norm(Res_ϕ[idx_ice])/sum(idx_ice) # with this error it also converges but more slowly
            Err_h .= abs.(Err_h .- h)
            err_h = norm(Err_h)
            # err_h   = norm(Res_h[idx_ice])/sum(idx_ice)
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

    return N, ϕ, h, qx, qy, ittot, Err_ϕ, Err_h, qx_ice, qy_ice
end

function plot_output(xc, yc, H, N, h, qx, qy, qx_ice, qy_ice)
    x_plt = [xc[1]; xc .+ (xc[2]-xc[1])]
    y_plt = [yc[1]; yc .+ (yc[2]-yc[1])]
    N[H .== 0.0] .= NaN
    h[H .== 0.0] .= NaN
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
    title("ϕ")
    # cross-sections of ϕ and h
    subplot(2, 2, 3)
    ind = size(N,2)÷2
    plot(xc, h[:, ind])
    title(join(["h at y = ", string(round(yc[ind], digits=1))]))
    subplot(2, 2, 4)
    plot(xc, N[:, ind])
    title(join(["ϕ at y = ", string(round(yc[ind], digits=1))]))

    qx_plot = copy(qx)
    qy_plot = copy(qy)
    qx_plot[qx_ice .== 0] .= NaN
    qy_plot[qy_ice .== 0] .= NaN
    figure()
    subplot(1, 2, 1)
    pcolor(qx_plot')
    colorbar()
    title("qx (m/s)")
    subplot(1, 2, 2)
    pcolor(qy_plot')
    colorbar()
    title("qy (m/s)")
end


end # module
