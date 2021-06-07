using Pkg
Pkg.activate(@__DIR__)
using SheetModel, Parameters
S = SheetModel

# SHMIP A3

# Define numerical domain and input parameters

# TODO: shouldn't be necessary in the end
Lx = 100e3
Ly = 20e3
Nx = 64
Ny = 32
dx = Lx/Nx      # grid size
dy = Ly/Ny
xc = LinRange(0, Lx, Nx) # vector of x-coordinates
yc = LinRange(0, Ly, Ny) # vector of y-coordinates

input_params = S.Para(
    lx = Lx,  # domain length in x-direction, m
    ly = Ly,  # domain length in y-direction, m
    nx = Nx,
    ny = Ny,
    dx = dx,
    ttot = 1000.0,
    dt = 100.0,

    H = (6 .*( sqrt.(xc.+5e3) .- sqrt(5e3) ) .+ 1 ) .* ones(Ny)', # ice thickness, m
    zb = zeros(Nx, Ny), # bed elevation, m
    m = 5.79e-09* ones(Nx-2, Ny-2),  # source term
)

# Initial condition
@unpack nx, ny, H = input_params
#ϕ0    = exp.(- 1e-2*(xc.-Lx/2).^2) * exp.(-1e-2*(yc.-Ly/2).^2)'
ϕ0 = 1e6/Lx * xc * ones(Ny)'
#ϕ0 = 1e6*rand(Nx, Ny)
#ϕ0 = 5e6 * rand(nx, ny)
h0 = 0.05/Lx * xc * ones(Ny)'

xc, yc, ϕ, h = S.runthemodel(input_params, ϕ0, h0);
S.plot_output(xc, yc, ϕ, h)
