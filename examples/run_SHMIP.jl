include("SHMIP_cases.jl")

# nit, ϕ, H, qx, qy = run_SHMIP("A1", Nx=64, Ny=32,
#                    γ_ϕ= 0.98, γ_h=0.92, dτ_ϕ_=1.04, dτ_h_= 1.06, #  -> nit = 1369
#                    # γ_ϕ= 0.9, γ_h=0.9, dτ_ϕ_=1.0, dτ_h_= 1.0 -> nit = 8876
#                    printtime=50, printit=1000, make_plot=true);

nit, ϕ, H, qx, qy, inputs = run_SHMIP("E1", Nx=64, Ny=32,
                    γ_ϕ= 0.2, γ_h=0.9, dτ_ϕ_=10.0, dτ_h_= 1.0, #  does not converge
                    printtime=1, printit=1000, make_plot=true);
