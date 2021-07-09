include("SHMIP_cases.jl")

 #nit, ϕ, h, H, qx, qy, inputs = run_SHMIP("A1", Nx=64, Ny=32,
 #                   γ_ϕ= 0.8, γ_h=0.9, dτ_ϕ_=1.0, dτ_h_= 7e-6,
 #                   printtime=50, printit=1000, make_plot=true);

 nit, ϕ, h, H, qx, qy, inputs = run_SHMIP("E1", Nx=64, Ny=32,
                     γ_ϕ= 0.7, γ_h=0.5, dτ_ϕ_= 1.0, dτ_h_= 1e-2,
                     printtime=1, printit=1000, make_plot=true);

                     #γ_ϕ= 0.9, γ_h=0.5, dτ_ϕ_= 0.1, dτ_h_= 5e-3,
