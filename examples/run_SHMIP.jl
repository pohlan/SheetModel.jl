include("SHMIP_cases.jl")


# inputs, outputs = run_SHMIP(test_case="A1", nx=1024, ny=1024, dt = 5e7, tsteps=1, # dt in seconds; 5e7 s to reach +/- steady-state
#                   γ_ϕ= 0.2, γ_h=0.2, dτ_ϕ_=1.0, dτ_h_= 7e-6,
#                   printtime=1, make_plot=false);

 inputs, outputs = run_SHMIP(test_case="E1", nx=256, ny=256, dt = 5e7, tsteps=1,
                     γ_ϕ= 0.7, γ_h=0.5, dτ_ϕ_= 1.0, dτ_h_= 1e-2, itMax=5*10^4,
                     printtime=1, make_plot=true);
