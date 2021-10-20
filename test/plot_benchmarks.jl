# ---------------------------------------------------#
# This script plots the output of "benchmarking.jl"  #
# ---------------------------------------------------#
using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using JLD2, UnicodePlots, DrWatson

# check if file exists and load
file = joinpath(@__DIR__, "benchmarks.jld2")
if !@isdefined benchmarks
    @assert isfile(file) "No benchmarks.jld2 file available to read."
    benchmarks = load(file)
end

# get current commit hash
if !@isdefined gitcommit
    gitcommit = gitdescribe()
end
if isdirty()
    gitcommit = split(gitcommit, "_")[1]
end

# draw plots
colors = [:green, :white, :red, :cyan, :magenta, :blue, :yellow]
symbols = Dict("A1" => :xcross,
               "E1" => :star5)
for (n, unit) in enumerate(keys(benchmarks))
    if haskey(benchmarks[unit], gitcommit)
        (host, PU) = split(unit, "_")
        case = benchmarks[unit][gitcommit][:SHMIP_case]
        steadyst = benchmarks[unit][gitcommit][:iterations] .> 10^3
        res = benchmarks[unit][gitcommit][:nx] .* benchmarks[unit][gitcommit][:ny]
        t_tot = benchmarks[unit][gitcommit][:run_time]
        T_eff = benchmarks[unit][gitcommit][:T_eff]
        if n == 1
            if any(steadyst)
                global plt_ttot = scatterplot(res[steadyst], t_tot[steadyst], marker=symbols[host], color=colors[PU], name=unit, xlabel="Grid size (nx*ny)", ylabel="Run time (s)")
            end
            global plt_Teff = scatterplot(res[.!steadyst], T_eff[.!steadyst], xscale=:log10, marker=symbols[case[n]], color=colors[n], name=unit, xlabel="Grid size (nx*ny)", ylabel="T_eff (GB/s)")
        else
            if any(steadyst)
                scatterplot!(plt_ttot, res[steadyst], t_tot[steadyst], marker=symbols[host], color=colors[PU], name=unit)
            end
            scatterplot!(plt_Teff, res[.!steadyst], T_eff[.!steadyst], xscale=:log10, marker=symbols[case[n]], color=colors[n], name=unit)
        end
    end
end

# make sure everything displays properly in the shell
print("Benchmarking results of the current commit: \n \n")
display(plt_ttot)
print("\n \n \n")
display(plt_Teff)
print("\n")
