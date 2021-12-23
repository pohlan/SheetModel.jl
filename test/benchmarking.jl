# ---------------------------------------------------------#
# Running benchmarking test cases and saving to jld2 file  #
# ---------------------------------------------------------#

using Base: AbstractVecOrTuple
using JLD2, Printf, CUDA, DataFrames
using SheetModel

# get the name of the server where the model was run
hostname = read(`hostname`, String)[1:end-1] # the last character is always "\n"

# get the type of processing unit and merge strings with hostname
if USE_GPU
    unitname = join([hostname "_GPU"])
else
    unitname = join([hostname, "_CPU_", string(Threads.nthreads()), "threads"])
end

# define the test sets
test_sets = Dict(# test case A1, 900 iterations not reaching steady-state, no error calculation; for T_eff vs. dof plot
                "Teff"      => [(nx=128,   ny=128,  itMax=900, warmup=10),
                                (nx=256,   ny=256,  itMax=900, warmup=10),
                                (nx=1024,  ny=1024, itMax=900, warmup=10),
                                (nx=2048,  ny=2048, itMax=900, warmup=10),
                                (nx=4096,  ny=4096, itMax=900, warmup=10),
                                (nx=8192,  ny=4096, itMax=900, warmup=10),
                                ],

                # test case A1, going into steady state; for wall time/#it vs. dof plot
                "std-state" => [(nx=32,  ny=32,  dτ_h=1e-5),
                                (nx=64,  ny=32,  dτ_h=1e-5),
                                (nx=64,  ny=64,  dτ_h=1e-5),
                                (nx=128, ny=64,  dτ_h=1e-5),
                                (nx=128, ny=128, dτ_h=7e-6),
                                (nx=256, ny=128, dτ_h=7e-6),
                                (nx=512, ny=256, dτ_h=6e-6)
                                ],

                # run different shmip cases to compare with reference from other models
                "A-suite"   => [(test_case="A1", dτ_h=1e-5),
                                (test_case="A2", dτ_h=1e-5),
                                (test_case="A3", dτ_h=1e-5),
                                (test_case="A4", dτ_h=1e-5),
                                (test_case="A5", dτ_h=1e-5)
                                ]
)

# get run_SHMIP function
include(joinpath(@__DIR__, "../examples/SHMIP_cases.jl"))

## run the benchmarking for each test set
dic = Dict()
for set in keys(test_sets)
    nt = (  wall_time = [],   # wall time in s
            T_eff     = [],   # effective memory throughput in GB/s
            nit       = [],   # number of iterations
            dof       = []    # degrees of freedeom
    )
    if set != "Teff"
        nt = merge(nt, (test_case = [],
                        ϕ         = [],
                        h         = [],
                        iters     = [],
                        errs_ϕ    = [],
                        errs_h    = [])
        )
    end

    for kwargs in test_sets[set]
        # for CPU only do the dof_Teff test runs
        if !USE_GPU && set != "Teff"
            continue
        end

        # run the model
        inputs, outputs = try
            run_SHMIP(;kwargs...);
        catch e
            if isa(e, CUDA.OutOfGPUMemoryError)
                @printf("nx =%d, ny=%d could not be run, not enough memory. \n", test_set.nx, test_set.ny)
                continue
            end
        end

        # push results to arrays
        push!(nt.wall_time, outputs.time_tot)
        push!(nt.T_eff,     outputs.T_eff)
        push!(nt.nit,       outputs.ittot)
        push!(nt.dof,       length(outputs.h))

        if set != "Teff"
            push!(nt.test_case, inputs.test_case)
            push!(nt.ϕ,         outputs.ϕ[2:end-1,end÷2])
            push!(nt.h,         outputs.h[2:end-1,end÷2])
            push!(nt.iters,     outputs.iters)
            push!(nt.errs_ϕ,    outputs.errs_ϕ)
            push!(nt.errs_h,    outputs.errs_h)
        end

    end
    dic[set] = DataFrame(;nt...)
end

# save to file
jldopen("test/bm_results.jld2", "a+") do file
    file[unitname] = dic
end
