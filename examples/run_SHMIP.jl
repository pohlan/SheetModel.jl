include("SHMIP_cases.jl")


 inputs, outputs = run_SHMIP(nx=128, ny=64,
                             γ_ϕ= 0.9, γ_h=0.8, dτ_h_=7e-6,
                             warmup=0, itMax=2*10^4, make_plot=true);

# inputs, outputs = run_SHMIP(test_case="E1", nx=256, ny=256, dt = 1e9,
#                     γ_ϕ= 0.7, γ_h=0.5, dτ_ϕ_= 1.0, dτ_h_= 7e-4,
#                     make_plot=true);
