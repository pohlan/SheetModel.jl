module SheetModel

using Base: Float64
using Parameters

export Para

@with_kw struct Para{F<:Float64, I<:Int64, Arr<:Array{Float64}, Lin<:LinRange{Float64}}
    # Scalars (one global value)
    g::F     = 1.0              # gravitational acceleration
    ρw::F    = 1.0              # water density
    ρi::F    = 910.0 / 1000.0
    α::F     = 1.25
    β::F     = 1.5
    k::F     = 1.0
    n::F     = 3.0
    A::F     = 1.0
    ev::F    = 0.5
    lr::F    = 1.0              # horizontal cavity spacing, spatially uniform
    hr::F    = 1.0              # bedrock bump height, spatially uniform

    # Numerical domain
    lx::F = 12.0       # domain size
    ly::F = 10.0
    nx::I = 64         # number of grids
    ny::I = 32
    dx::F = lx/nx      # grid size
    dy::F = ly/ny
    xc::Lin = LinRange(dx/2, lx-dx/2, nx) # vector of x-coordinates
    yc::Lin = LinRange(dy/2, ly-dy/2, ny) # vector of y-coordinates

    # Field parameters (defined on every grid point)
    zs::Arr = 1/lx * xc * ones(ny)' # surface elevation
    zb::Arr = zeros(nx, ny) # bed elevation
    m::Arr  = zeros(nx-2, ny-2)  # source term
    ub::Arr = zeros(nx, ny) # basal sliding speed

    # Dimensionless numbers
    Σ::F   = 1e2
    Γ::F   = 1e5
    Λ::F   = 0.5
    r_ρ::F = 1000 / 910 # ρw / ρi

    # Physical time stepping
    ttot::F   = 2.0        # total simulation time
    dt::F     = 0.1        # physical time step

    # Pseudo-time iteration
    tol::F    = 1e-6       # tolerance
    itMax::I  = 10^3        # max number of iterations
    damp::F   = 1-41/nx    # damping (this is a tuning parameter, dependent on e.g. grid resolution) # TODO: define how?
    dτ::F     = (1.0/(dx^2/k/2.1) + 1.0/dt)^-1 # pseudo-time step; TODO: define how?
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





end # module
