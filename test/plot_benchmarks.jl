# ---------------------------------------------------#
# This script plots the output of "benchmarking.jl"  #
# ---------------------------------------------------#

using JLD2, UnicodePlots, DrWatson

# check if file exists and load
file = joinpath(@__DIR__, "benchmarks.jld2")
if !@isdefined benchmarks
    @assert isfile(file) "No benchmarks.jld2 file available to read."
    benchmarks = load(file)
end

# get current commit hash
gitcommit = gitdescribe()

# draw plots
colors = Dict("GPU"           => :green,
              "CPU-1threads"  => :white,
              "CPU-4threads"  => :red,
              "CPU-16threads" => :cyan,
              "CPU-32threads" => :magenta)
symbols = Dict("achtzack01"      => :xcross,
               "octopus"         => :star5,
               "annegret-laptop" => :hexagon)
for (n, unit) in enumerate(keys(benchmarks))
    (host, PU) = split(unit, "_")
    res = benchmarks[unit][gitcommit][:nx] .* benchmarks[unit][gitcommit][:ny]
    t_tot = Float64.(benchmarks[unit][gitcommit][:run_time])
    T_eff = Float64.(benchmarks[unit][gitcommit][:T_eff])
    if n == 1
        global plt_ttot = scatterplot(res, t_tot, marker=symbols[host], color=colors[PU], name=unit, xlabel="Grid size (nx*ny)", ylabel="Run time (s)")
        global plt_Teff = scatterplot(res, T_eff, marker=symbols[host], color=colors[PU], name=unit, xlabel="Grid size (nx*ny)", ylabel="T_eff (GB/s)")
    else
        scatterplot!(plt_ttot, res, t_tot, marker=symbols[host], color=colors[PU], name=unit)
        scatterplot!(plt_Teff, res, T_eff, marker=symbols[host], color=colors[PU], name=unit)
    end
end

# make sure everything displays properly in the shell
print("Benchmarking results of the current commit: \n \n")
display(plt_ttot)
print("\n \n \n")
display(plt_Teff)
print("\n")
