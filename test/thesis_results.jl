# Reproducing the figures from the Master thesis

## produce the results for Fig. 3.2 and 3.5 (takes > 1h)
#include(joinpath(@__DIR__, "../test/benchmarking.jl"))

## produce the results for Fig. 3.3
# include(joinpath(@__DIR__, "../examples/SHMIP_cases.jl"))
# inout1 = run_SHMIP(nx=64,  ny=32,  dτ_h=1e-5)
# inout2 = run_SHMIP(nx=1024, ny=1024, dτ_h=5e-6, itMax=4*10^5)
# save("test/error_count.jld2", "inout", (inout1, inout2))

# plotting functions
include(joinpath(@__DIR__, "../test/plotting_fcts.jl"))

# Figure 3.2
############
plot_fieldresults("jsb")

# Figure 3.3
############
plot_error()

# Figure 3.4
############
plot_benchmarks()
