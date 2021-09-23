include("SHMIP_cases.jl")


 inputs, outputs = run_SHMIP("A1", Nx=128, Ny=64, dt = 5e7, tsteps=1, # dt in seconds; 5e7 s to reach +/- steady-state
                   γ_ϕ= 0.45, γ_h=0.45, dτ_ϕ_=1.0, dτ_h_= 7e-6,
                   printtime=1, make_plot=false);

# inputs, outputs = run_SHMIP("E1", Nx=64, Ny=32, dt = 5e7, tsteps=1,
#                     γ_ϕ= 0.7, γ_h=0.5, dτ_ϕ_= 1.0, dτ_h_= 1e-2,
#                     printtime=1, make_plot=true);
