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
(Note: I don't know yet how this works with -gpu as this is an input to the script and not to julia in general. One can manually set const USE_GPU = true in `src/SheetModel.jl`.)

To get optimal performance, one can replace `julia` by `julia -O3 --check-bounds=no`.

### Where to find what
**src/**
- All the physics and numerics of the model is described in `src/modelonly.jl`.
- Helper functions defining struct types, performing the scaling etc. are in `src/helpers.jl`

**examples/**
- `examples/run_SHMIP.jl` calls `examples/SHMIP_cases.jl`, where all the SHMIP cases are defined.
- Note that `examples/Antarctica_case.jl` and `examples/minimal_example` are not working currently.

**test/**
- `test/runtests.jl` tests whether the current implementation agrees with a reference solution that has been saved using `test/test_references.jl` in an earlier commit
- `test/benchmarking.jl` produces the results for Fig. 3.2 and 3.4 in the thesis
- `test/compare_params.jl` performs a systematic search of gamma and dtau_h parameter space
- `test/thesis_results.jl` calls the plotting functions defined in `test/plotting_fcts.jl` to produce Figures 3.2, 3.3 and 3.4 from the thesis

### Documentation
The thesis can be found here https://github.com/pohlan/MA-notes/tree/main/main_document