using Base: AbstractVecOrTuple
using DrWatson, JLD2
using SheetModel

# get current commit hash
gitcommit = gitdescribe()
if gitcommit[end-4:end] == "dirty"
    print("Press enter to continue with dirty repo or ^C to abort.")
    readline() # will stop the execution and only continue when pressing enter
    print("_dirty suffix was removed.")
    gitcommit = split(gitcommit, "_")[1]
end

# get the name of the server where the model was run
hostname = read(`hostname`, String)[1:end-1] # the last character is always "\n"

# get the type of processing unit and merge strings with hostname
if USE_GPU
    unitname = join([hostname "_GPU"])
else
    unitname = join([hostname, "_CPU-", string(Threads.nthreads()), "threads"])
end

# run the model
scriptname = "examples/run_SHMIP.jl"
include(joinpath(@__DIR__, "../" * scriptname))

# save output in dictionary
d = Dict(
    :SHMIP_case => inputs.SHMIP_case,
    :total_time => outputs.time_tot,         # in s
    :T_eff => outputs.T_eff,                 # effective memory throuhput in GB/s
    :nx => inputs.input_params.nx,
    :ny => inputs.input_params.ny,
    :iterations => outputs.ittot,
)

# load file with previous benchmarks
file = joinpath(@__DIR__, "benchmarks.jld2")
if isfile(file)
    benchmarks = load(file)
else
    benchmarks = Dict()
end

# check if unitname and gitcommit entries in the dictionaries already exist
if !haskey(benchmarks, unitname)
    benchmarks[unitname] = Dict()
end
if !haskey(benchmarks[unitname], gitcommit)
    benchmarks[unitname][gitcommit] = Dict(k=>[] for k in keys(d) )
end

# push to each entry of the dictionary
for key in keys(benchmarks[unitname][gitcommit])
    push!(benchmarks[unitname][gitcommit][key], d[key])
end

# write updated dictionary back to file
save(file, benchmarks)
