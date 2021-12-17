include("SHMIP_cases.jl")


 inputs, outputs = run_SHMIP(nx=64, ny=32,
                             γ_ϕ= 0.8, γ_h=0.8, dτ_h_=7e-6,
                             warmup=0, itMax=10^6, make_plot=true);

# inputs, outputs = run_SHMIP(test_case="E1", nx=256, ny=256, dt = 1e9,
#                     γ_ϕ= 0.7, γ_h=0.5, dτ_ϕ_= 1.0, dτ_h_= 7e-4,
#                     make_plot=true);
