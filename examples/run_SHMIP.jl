include("SHMIP_cases.jl")


 inputs, outputs = run_SHMIP(test_case="A1", nx=64, ny=32, dt = 1e9, # dt in seconds; 5e7 s to reach +/- steady-state
                             γ_ϕ= 0.91, γ_h=0.91, dτ_ϕ_=1.0, dτ_h_=1e-5, itMax=3*10^4,
                             update_h_only=false, make_plot=true);

# inputs, outputs = run_SHMIP(test_case="E1", nx=256, ny=256, dt = 1e9,
#                     γ_ϕ= 0.7, γ_h=0.5, dτ_ϕ_= 1.0, dτ_h_= 7e-4,
#                     make_plot=true);
