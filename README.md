# SheetModel.jl

### Run the model
Execute the script `examples/run_SHMIP.jl`, in the bash with:
```
../SheetModel.jl$ julia examples/run_SHMIP.jl -gpu       # for running it on GPU
../SheetModel.jl$ julia -t16 examples/run_SHMIP.jl       # running it on CPU with 16 threads
```
To see the plots, one must run the script from the Julia REPL:
```
../SheetModel.jl$ julia -t16
julia> include("examples/run_SHMIP.jl")
```
(Note: I don't know yet how this works with -gpu as this is an input to the script and not to julia in general.)

To get optimal performance, one can replace `julia` by `julia -O3 --check-bounds=no`.

### Benchmark the current version
Execute the script `test/benchmarking.jl`, in the same manner as above.
This script will automatically call `test/plot_benchmarks.jl`, which produces Unicode plots of all the results from the current commit. Since these plots are Unicode they can also be shown when the script is called from the shell.
If the benchmarking has already been done and one just wants to plot the results, `test/plot_benchmarks.jl` can be executed directly.

### Where to find what
All the physics and numerics of the model is described in `src/modelonly.jl`.

**Documentation** can be found at https://github.com/pohlan/MA-notes. In particular, the `spatial_discretisation.tex` explains which grids are used, how derivatives, gradients etc. are discretised and how boundary conditions are imposed.