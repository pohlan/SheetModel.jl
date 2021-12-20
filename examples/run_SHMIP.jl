include("SHMIP_cases.jl")

 inputs, outputs = run_SHMIP(nx=512, ny=512, test_case="A1",
                             γ_ϕ= 0.3, γ_h=0.3, dτ_h_=1e-8,
                             warmup=0, itMax=10^6, make_plot=true);

#  inputs, outputs = run_SHMIP(test_case="E4", nx=64, ny=64,
#                      γ_ϕ= 0.7, γ_h=0.5, dτ_ϕ_= 1.0, dτ_h_= 7e-4,
#                      make_plot=true);

# test cases E1 to E4 work with m = 5e-10
