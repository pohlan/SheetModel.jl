# SheetModel.jl

- **Run the model** by executing the script `examples/run_SHMIP.jl`.

- **Benchmark** the current version by running `test/benchmarking.jl`

- **Plot benchmarking results** with `test/plot_benchmarks.jl`

The most relevant stuff is in `src/modelonly.jl`

**Documentation** can be found at https://github.com/pohlan/MA-notes. In particular, the `spatial_discretisation.tex` explains which grids are used, how derivatives, gradients etc. are discretised and how boundary conditions are imposed.