__precompile__(false)
module SheetModel

using ParallelStencil
export Para, USE_GPU

# Take command line argument to decide which processing unit should be used.
const USE_GPU = any(ARGS .== "-gpu") # defaults to CPU

# Initiate ParallelStencil
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end

include("helpers.jl")      # defining struct types, (de-)scaling functions, plotting functions etc.
include("modelonly.jl")    # the actual model

end # module
