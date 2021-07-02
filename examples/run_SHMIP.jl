include("SHMIP_cases.jl")

nit, ϕ, H, qx, qy, inputs = run_SHMIP("A1", Nx=64, Ny=32,
                   γ_ϕ= 0.2, γ_h=0.2, dτ_ϕ_=1.0, dτ_h_= 10.0,
                   printtime=50, printit=1000, make_plot=true);

# nit, ϕ, H, qx, qy, inputs = run_SHMIP("E1", Nx=64, Ny=32,
#                     γ_ϕ= 0.05, γ_h=0.5, dτ_ϕ_= 100.0, dτ_h_= 10.0, #  does not converge
#                     printtime=1, printit=1000, make_plot=true);
