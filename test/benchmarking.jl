# ------------------------------------------------------------------------------------------------------------------#
# This script runs the model over run_SHMIP.jl and saves the parameters                                             #
# in a nested dictionary, which is then saved in the "benchmarks.jld2" file.                                        #
# Structure:                                                                                                        #
#                                                                                                                   #
# :Achtzack01_GPU => :commit1 => :SHMIP_case   = ["A1", "A1", "A2", ...]                                            #
#                                :run_time     = [28.1, 34.5, 38.1, ...]  # absolute time, in s                     #
#                                :T_eff        = [20.3, 19.1, 18.7, ...]  # effective memory throughput, GB/s       #
#                                :nx           = [1024, 4096, 4096, ...]                                            #
#                                :ny           = [512, 1024, 4096, ...]                                             #
#                                :iterations   = [10^3, 28200, 10^3, ...]                                           #
#                    :commit2  => ...                                                                               #
#                                                                                                                   #
# :Achtzack01_CPU-1threads      => :commit1  => ...                                                                 #
# :Achtzack01_CPU-16threads     => ...                                                                              #
# :Achtzack01_CPU-32threads     => ...                                                                              #
# :annegret-laptop_CPU-1threads => ...                                                                              #
# :Octopus_GPU                  => ...                                                                              #
# ------------------------------------------------------------------------------------------------------------------#

using Base: AbstractVecOrTuple
using DrWatson, JLD2, Printf, CUDA
using SheetModel

# get the name of the server where the model was run
hostname = read(`hostname`, String)[1:end-1] # the last character is always "\n"

# get the type of processing unit and merge strings with hostname
if USE_GPU
    unitname = join([hostname "_GPU"])
else
    unitname = join([hostname, "_CPU-", string(Threads.nthreads()), "threads"])
end

# load file with previous benchmarks
file = joinpath(@__DIR__, "benchmarks.jld2")
if isfile(file)
    benchmarks = load(file)
else
    benchmarks = Dict()
end

# get current commit hash and check that repo is clean
gitcommit = gitdescribe()
@assert !isdirty() "Make sure the repo is clean before benchmarking."

# check if unitname and gitcommit entries in the dictionaries already exist
if !haskey(benchmarks, unitname)
    benchmarks[unitname] = Dict()
end
if haskey(benchmarks[unitname], gitcommit)
    @printf("This commit has already been benchmarked by %s. Press enter to continue or ^C to abort.", unitname)
    readline() # will stop the execution and only continue when pressing enter
else
    benchmarks[unitname][gitcommit] = Dict()
end

# define the test sets
include(joinpath(@__DIR__, "../examples/SHMIP_cases.jl"))
test_sets = [# 10^3 iterations without reaching steady state (and without calculating errors)
            (test_case="A1", nx=128,   ny=128,  itMax=10^3),
            (test_case="A1", nx=256,   ny=256,  itMax=10^3),
            (test_case="A1", nx=512,   ny=512,  itMax=10^3),
            (test_case="A1", nx=1024,  ny=1024, itMax=10^3),
            (test_case="A1", nx=2048,  ny=2048, itMax=10^3),
            (test_case="A1", nx=4096,  ny=4096, itMax=10^3),
            #(test_case="A1", nx=8192,  ny=8192, itMax=10^3),
            #(test_case="A1", nx=16384, ny=8192, itMax=10^3),

             # going into steady state
            (test_case="A1", nx=128, ny=128, itMax=10^6, γ_ϕ= 0.9, γ_h=0.8, dτ_ϕ_=1.0, dτ_h_=6e-6),
            (test_case="A1", nx=256, ny=256, itMax=10^6, γ_ϕ= 0.9, γ_h=0.8, dτ_ϕ_=1.0, dτ_h_=6e-6),
            (test_case="A1", nx=512, ny=512, itMax=10^6, γ_ϕ= 0.9, γ_h=0.8, dτ_ϕ_=1.0, dτ_h_=6e-6),
            #(test_case="A1", nx=1024, ny=1024, itMax=10^6, γ_ϕ= 0.9, γ_h=0.8, dτ_ϕ_=1.0, dτ_h_=8e-7),   # takes more than 30 min and the result is still a tiny bit off
            #(test_case="A1", nx=4096, ny=2048, itMax=10^6),
            #(test_case="A3", nx=1024, ny=512, itMax=10^5),
            #(test_case="A3", nx=4096, ny=1024, itMax=10^5),
            #(test_case="F1", nx=1024, ny=512, itMax=10^5),
            #(test_case="F1", nx=4096, ny=2048, itMax=10^5),
             ]

## run the benchmarking for each test set
for test_set in test_sets
    # run the model
    inputs, outputs = try
        run_SHMIP(;test_set...);
    catch e
        if isa(e, CUDA.OutOfGPUMemoryError)
            @printf("nx =%d, ny=%d could not be run, not enough memory. \n", test_set.nx, test_set.ny)
            continue
        end
    end

    # save output in dictionary
    d = Dict(
        :SHMIP_case => inputs.SHMIP_case,
        :run_time => outputs.time_tot,
        :T_eff => outputs.T_eff,
        :nx => inputs.input_params.nx,
        :ny => inputs.input_params.ny,
        :iterations => outputs.ittot,
    )

    # push to each entry of the dictionary
    if !haskey(benchmarks[unitname][gitcommit], :SHMIP_case)
        [benchmarks[unitname][gitcommit][k] = typeof(d[k])[] for k in keys(d)]
    end
    for key in filter(x->x!=:gitpatch, keys(benchmarks[unitname][gitcommit]))
        push!(benchmarks[unitname][gitcommit][key], d[key])
    end
end

# write updated dictionary back to file
save(file, benchmarks)

include("plot_benchmarks.jl")
