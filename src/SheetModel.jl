__precompile__(false)
module SheetModel

# ParallelStencil, decide whether to use CPU or GPU
const USE_GPU = false
using ParallelStencil
#using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end


export Para

include("helpers.jl")
include("modelonly.jl")

end # module
