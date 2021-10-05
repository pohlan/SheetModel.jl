__precompile__(false)
module SheetModel

# ParallelStencil, decide whether to use CPU or GPU
const USE_GPU = true
using ParallelStencil #, CUDA
#CUDA.allowscalar(false) # don't allow scalar indexing which is not GPU compatible
#using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end


export Para, USE_GPU

include("helpers.jl")
include("modelonly.jl")

end # module
