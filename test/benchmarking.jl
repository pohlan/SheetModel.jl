using DrWatson, JLD2
using SheetModel

# get current commit hash
gitcommit = gitdescribe()

# load the correct file with previous benchmarks; there is one file for the GPU and one for the CPU
if USE_GPU
    file = joinpath(@__DIR__, "benchmarks_GPU.jld2")
else
    file = joinpath(@__DIR__, "benchmarks_CPU.jld2")
end
if isfile(file)
    benchmarks = load(file)
else
    benchmarks = Dict()
end

# only continue if the current version hasn't been benchmarked yet
@assert !haskey(benchmarks, gitcommit) "The current version has already been benchmarked for the chosen processing unit."

# run the model
scriptname = "examples/run_SHMIP.jl"
include(joinpath(@__DIR__, "../" * scriptname))

# get the name of the device where the model was run
hostname = read(`hostname`, String)

# save output in dictionary
d = Dict(
    :hostname => hostname,
    :which => scriptname,
    :total_time => outputs.time_tot,            # in s
    :T_eff => outputs.T_eff,                 # effective memory throuhput in GB/s
    :nx => inputs.input_params.nx,
    :ny => inputs.input_params.ny,
    :iterations => outputs.ittot,
)

# make a new entry to the loaded dictionary
benchmarks[gitcommit] = d

# write updated dictionary back to file
save(file, benchmarks)
