using Base: Float64, Int64
using LinearAlgebra, Parameters, Statistics, PyPlot, LaTeXStrings, Infiltrator

# misc fudge factors to avoid dividing by zero
small = eps(Float64)

"""
Create a vector of lenght nx where boundary points are zeros and interior points ones.
Used to achieve zero ice thickness (H=0) at ghost points
"""
function ghost!(A)
    A[[1, end], :] .= 0.
    A[:, [1, end]] .= 0.
    return
end



function make_model_input(H, zb, Lx, Ly, dx, dy, ttot, dt, itMax, tol, γ_ϕ, γ_h, dτ_ϕ_, dτ_h_, ev, ev_num, ϕ_init, h_init, calc_m, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux)
    # make outermost rows and columns to ghost points
    ghost!(H)
    ghost!(zb)
    ghost!(ϕ_init)
    ghost!(h_init)
    ghost!(ice_mask)

    # struct of input parameters
    params_struct = model_input(;H, zb, Lx, Ly, dx, dy, ttot, dt, itMax, tol, γ_ϕ, γ_h, dτ_ϕ_, dτ_h_, ev, ev_num)

    return (;params_struct, ϕ_init, h_init, calc_m, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux)
end


"""
Create struct including all model parameters, physical and numerical
"""
@with_kw struct model_input{T} @deftype Float64
    # Physics
    g      = 9.81              # gravitational acceleration, m/s^2
    ρw     = 1000.0            # water density, kg/m^3
    ρi     = 910.0             # ice density, kg/m^3
    α      = 1.25              # first sheet flow exponent
    β      = 1.5               # second sheet flow exponent
    k      = 0.005             # sheet conductivity, m^(7/4)kg^(-1/2)
    n      = 3.0               # Glen's flow law exponent
    A      = 3.375e-24         # ice flow constant, Pa^(-n)s^(-1)
    ev     = 0              # englacial void ratio; SHMIP: 0 for ice-sheet, 1e-3 for valley glacier
    ev_num = 1e-1              # regularisation void ratio
    lr     = 2.0               # horizontal cavity spacing, m
    hr     = 0.1               # bedrock bump height, m
    ub     = 1e-6              # basal sliding speed, m/s
    H::T                      # ice thickness, m
    zb::T                     # bed elevation, m

    # Domain
    Lx
    Ly
    dx
    dy

    # Physical time stepping
    ttot                      # total simulation time
    dt                        # physical time step

    # Pseudo-time iteration
    tol   = 1e-6              # tolerance
    itMax = 10^6              # max number of iterations
    γ_ϕ   = 1e-3              # damping parameter for ϕ update
    γ_h   = 0.8               # damping parameter for h update
    dτ_ϕ_ = 1e6               # scaling factor for dτ_ϕ
    dτ_h_ = 50.0              # scaling factor for dτ_h

    # Dimensionless numbers
    Ψ   = NaN
    Σ   = NaN
    Γ   = NaN
    Λ   = NaN
end
Broadcast.broadcastable(p::model_input) = Ref(p)


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
function array_allocation(nx, ny)
    qx     = @zeros(nx-1, ny  )
    qy     = @zeros(nx,   ny-1)
    d_eff  = @zeros(nx,   ny  )
    Λ_m    = @zeros(nx,   ny  )
    N      = @zeros(nx,   ny  )
    dϕ_dτ  = @zeros(nx,   ny  )
    dh_dτ  = @zeros(nx,   ny  )
    Res_ϕ  = @zeros(nx-2, ny-2)
    Res_h  = @zeros(nx-2, ny-2)
    return qx, qy, d_eff, Λ_m, N, dϕ_dτ, dh_dτ, Res_ϕ, Res_h
end

"""
Convert input parameters to non-dimensional quantities.
Convert arrays of H and zb to correct types depending on whether CPU or GPU is used.
"""
function scaling(p::model_input, ϕ_init, h_init, calc_m, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux)
    @unpack g, ρw, ρi, k, A, lr, hr,
            H, zb, ub,
            Lx, Ly, dx, dy,
            ttot, dt,
            α, β, n, ev,
            Ψ, Σ, Γ, Λ = p

    # Scaling factors
    g_ = g
    ρw_ = ρw
    ρi_ = ρi
    k_ = k
    A_ = A
    lr_ = lr
    h_ = hr
    xy_ = max(Lx, Ly)

    H_ = maximum(H)
    zb_ = H_ * ρi_ / ρw_
    m_ = mean(calc_m(1:size(H, 1), (1:size(H, 2))', 0.0)) # for time-dependent input: temporal peak
    ub_ = ub

    ϕ_ = g_ * H_ * ρi_
    N_ = ϕ_
    q_ = k_ * h_^α * ( ϕ_  / xy_ )^(β-1)
    t_ = xy_ * h_ / q_
    vo_ = ub_ * h_ / lr_
    vc_ = 2 / n^n * A_ * h_ * N_ ^n

    if any(.!isnan.([Σ, Γ, Λ]))
        @warn "Ψ, Σ, Γ and Λ have already been assigned."
    end

    # Dimensionless parameters
    scaled_params = reconstruct(model_input{Data.Array}, p,
        # Scalars (one global value)
        g = g / g_,
        ρw = ρw / ρw_,
        ρi = ρi / ρi_,
        k = k / k_,
        A = A / A_,
        lr = lr / lr_,
        hr = hr / h_,
        ub = ub / ub_,

        # Field parameters (defined on every grid point)
        H = Data.Array(H ./ H_),
        zb = Data.Array(zb ./ zb_),

        # Numerical domain
        dx = dx ./ xy_,
        dy = dy ./ xy_,

        ttot = ttot / t_,
        dt = dt / t_,

        Ψ = ev * ϕ_ / (ρw_ * g_ * h_),
        Σ = vo_ * xy_ / q_,
        Γ = vc_ * xy_ / q_,
        Λ = m_ * xy_ / q_
        )::model_input{Data.Array}

    # initial conditions
    ϕ_init = Data.Array(ϕ_init ./ ϕ_)
    h_init = Data.Array(h_init ./ h_)

    # source term function
    calc_m_sc = (ix, iy, t) -> calc_m(ix, iy, t) / m_

    # masks
    ice_mask = Data.Array(ice_mask)
    bc_diric = Data.Array(bc_diric)
    bc_no_xflux = Data.Array(bc_no_xflux)
    bc_no_yflux = Data.Array(bc_no_yflux)

    return scaled_params, ϕ_init, h_init, calc_m_sc, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux, ϕ_, N_, h_, q_
end

mat(x) = Matrix{Float64}(x) # shorter notation for use in reformat() function

"""
Convert unitless output parameters back to dimensional quantities.
Convert all arrays to Matrix{Float64} type.
"""
function descale(output::model_output, N_, ϕ_, h_, q_)
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
