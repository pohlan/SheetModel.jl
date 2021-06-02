using Pkg, PyPlot, Printf, LinearAlgebra
Pkg.activate(@__DIR__)
using SheetModel
S = SheetModel

# enable plotting by default
if !@isdefined do_visu; do_visu = true end
pygui(true)


@views function runthemodel(; do_visu=true)

    # Numerical domain
    lx, ly     = 10.0, 12.0       # domain size
    nx, ny     = 64, 32        # numerical grid resolution
    dx     = lx/nx      # grid size
    dy     = ly/ny
    xc     = LinRange(dx/2, lx-dx/2, nx)
    yc     = LinRange(dy/2, ly-dy/2, ny)

    # Dimensionless numbers
    sigma = 0
    gamma = 0
    lambda = 0
    r_rho = 1000 / 910 # rhow / rhoi

    # Physical parameters
    k = 1
    n = 3
    A = 1
    g = 1
    rhow = 1
    rhoi = rhow / r_rho
    ev = 0.5
    lr = 1 # horizontal cavity spacing, spatially uniform
    hr = 1 # bedrock bump height, spatially uniform

    zs = 1/lx * xc * ones(ny)' # surface elevation
    zb = zeros(nx, ny) # bed elevation
    m = zeros(nx-2, ny-2) # source term
    ub = zeros(nx, ny) # basal sliding speed

    # Time stepping
    ttot   = 2.0        # total simulation time
    dt     = 0.1        # physical time step

    # Pseudo-time iteration
    tol    = 1e-6       # tolerance
    itMax  = 1e3        # max number of iterations
    damp   = 1-41/nx    # damping (this is a tuning parameter, dependent on e.g. grid resolution)
    dtau   = (1.0/(dx^2/k/2.1) + 1.0/dt)^-1 # iterative "timestep"

    # Array allocation
    qx     = zeros(nx-1, ny-1)
    qy     = zeros(nx-1, ny-1)

    dphi_dtau = zeros(nx-2, ny-2)
    Res_phi   = zeros(nx-2, ny-2)

    dh_dtau = zeros(nx-2, ny-2)
    Res_h     = zeros(nx-2, ny-2)

    # Initial condition
    # phi0    = exp.(- 1e-2*(xc.-lx/2).^2) * exp.(-1e-2*(yc.-ly/2).^2)'
    # phi0 = 1/lx * xc * ones(ny)'
    phi0 = rand(nx, ny)
    phi_old   = copy(phi0)
    phi      = copy(phi0)

    h0 = 0.5 * ones(nx, ny)
    h_old = copy(h0)
    h = copy(h0)

    t = 0.0; it = 0; ittot = 0

    # Physical time loop
    while t<ttot
        iter = 0; err_phi = 2*tol; err_h = 2*tol;

        # Pseudo-transient iteration
        while max(err_phi, err_h) >tol && iter<itMax
            dphi_dx, dphi_dy = diff(phi, dims=1) ./ dx, diff(phi, dims=2) ./ dy   # hydraulic gradient
            qx, qy     = S.calc_q(h[1:end-1, :], dphi_dx, k), S.calc_q(h[:, 1:end-1], dphi_dy, k) # TODO: which h values to take?!

            pw     = phi .- rhow .* g .* zb                   # water pressure
            N      = rhoi * g * (zs .- zb) .- pw            # effective pressure
            vo     = ub .* (hr .- h) ./ lr                  # opening rate
            vc     = 2 / n^n * A * h .* abs.(N).^(n-1) .* N  # closure rate

            # calculate residuals
            Res_phi        =    - ev/(rhow*g) * (phi[2:end-1, 2:end-1] .- phi_old[2:end-1, 2:end-1])/dt .-    # dhe/dt
                                (diff(qx, dims=1)[:, 2:end-1]/dx .+ diff(qy, dims=2)[2:end-1, :]/dy)    .+    # div(q)
                                sigma * vo[2:end-1, 2:end-1] .- gamma * vc[2:end-1, 2:end-1]            .-    # dh/dt
                                lambda * m                                                                    # source term
            Res_h          =    (h[2:end-1, 2:end-1] .- h_old[2:end-1, 2:end-1]) / dt  .-
                                sigma * vo[2:end-1, 2:end-1] .+ gamma * vc[2:end-1, 2:end-1]

            # damped rate of change
            dphi_dtau      = Res_phi .+ damp * dphi_dtau
            dh_dtau        = Res_h .+ damp * dh_dtau

            # update fields
            phi[2:end-1, 2:end-1]  = phi[2:end-1, 2:end-1] + dtau * dphi_dtau   # update rule, sets the BC as phi[1]=phi[end]=0
            h[2:end-1, 2:end-1]    = h[2:end-1, 2:end-1] + dtau * dh_dtau

            err_phi = norm(Res_phi)/length(Res_phi)
            err_h   = norm(Res_h)/length(Res_h)
            iter += 1
        end
        ittot += iter; it += 1; t += dt

        phi_old = phi
        h_old = h
    end

    # Visualize
    if do_visu
        figure()
        subplot(1, 2, 1)
        pcolor(xc, yc, phi0')
        subplot(1, 2, 2)
        pcolor(xc, yc, phi')
    end

end

runthemodel(; do_visu=do_visu);