# Reproducing the figures from the Master thesis

# plotting functions
#include(joinpath(@__DIR__, "../test/plotting_fcts.jl")
#include(joinpath(@__DIR__, "../examples/SHMIP_cases.jl"))


# figure
# for figure xxx
# inout1 = run_SHMIP(nx=64,  ny=32,  dτ_h_=1e-5)
# inout2 = run_SHMIP(nx=1024, ny=1024, dτ_h_=5e-6, itMax=4*10^5)
# save("test/error_count.jld2", "inout", (inout1, inout2))
plot_error()

# figure
#plot_benchmarks()

