using Pkg
Pkg.activate(@__DIR__)
using SheetModel
S = SheetModel


# Define numerical domain and input parameters
Lx = 12.0
Ly = 10.0
Nx = 64
Ny = 32
dx = Lx/Nx      # grid size
dy = Ly/Ny
xc = LinRange(dx/2, Lx-dx/2, Nx) # vector of x-coordinates
yc = LinRange(dy/2, Ly-dy/2, Ny) # vector of y-coordinates

input_params = S.Para(
        lx = Lx,
        ly = Ly,
        nx = Nx,
        ny = Ny,
        ttot = 2.0,
        dt = 0.1,

        H = 1/Lx .* xc .* ones(Ny)', # surface elevation
        zb = zeros(Nx, Ny), # bed elevation
        m = 0.1* ones(Nx-2, Ny-2),  # source term
        ub = 1e-5 * ones(Nx, Ny) # basal sliding speed
    )

# Initial condition
ϕ0    = exp.(- 1e-2*(xc.-Lx/2).^2) * exp.(-1e-2*(yc.-Ly/2).^2)'
#ϕ0 = 1/Lx * xc * ones(Ny)'
#ϕ0 = rand(Nx, Ny)
h0 = 0.5 * ones(Nx, Ny)

xc, yc, ϕ0, ϕ = S.runthemodel(input_params, ϕ0, h0);
S.plot_output(xc, yc, ϕ0, ϕ)
