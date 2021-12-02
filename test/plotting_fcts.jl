using JLD2, PyPlot, LaTeXStrings, Parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 24

"""
Plot a map of wall time vs damping (γ) and pseudo-time step (dτ_h)
"""
function plot_para_comp()
    fig, axs = subplots(1, 3, figsize=(23,10))
    fig.subplots_adjust(hspace=0.3, wspace=0.5)
    p = (nothing for i=1:100)
    for (n, dof) in enumerate(["nx64_ny32", "nx128_ny64", "nx256_ny128"])
        γ, dτ_h, wt = load("test/params_comp_" * dof * ".jld2","γs", "dτ_hs", "wall_time")

        dγ       = γ[2]-γ[1]
        dτ       = dτ_h[2] - dτ_h[1]
        γ_plt    = [γ[1] - 0.5dγ; γ .+ 0.5dγ]
        dτ_h_plt = [dτ_h[1] - 0.5dτ; dτ_h .+ 0.5dτ]

        p = axs[n].pcolormesh(γ_plt, dτ_h_plt, wt', cmap="viridis_r"; norm=PyPlot.matplotlib.colors.LogNorm())
        axs[n].set_xlabel(L"Damping parameter $γ$ = $γ_h$ = $γ_ϕ$", labelpad=20)
        axs[n].set_ylabel(L"Pseudo-time step $dτ_h$", labelpad=20)

        axs[n].set_title("nx=" * split(dof, "_")[1][3:end] * ", ny=" * split(dof, "_")[2][3:end])
    end
    cb = fig.colorbar(p, ax=fig.get_axes())
    cb.set_ticks([cb.vmin, cb.vmax])
    cb.set_ticklabels(["slow", "fast"])
    cb.set_label("Wall time, normalised")
end

#plot_para_comp()

function plot_error(errs_ϕ, errs_h, iters; nx, ny)
    figure(figsize=(13, 10))
    semilogy(iters, errs_ϕ, label=L"ϕ", color="gold", lw=2)
    semilogy(iters, errs_h, label=L"h", color="deepskyblue", lw=2)
    title("nx = " * string(nx) * ", ny = " * string(ny))

    xlabel("Iteration")
    ylabel(L"RMS($f_h$)")
    legend()
end

# after running run_SHMIP
#plot_error(outputs.errs_ϕ, outputs.errs_h, outputs.iters; nx=size(outputs.h, 1), ny=size(outputs.h, 2))

function plot_fieldresults(;inputs, outputs)
    @unpack Lx, Ly, dx, dy, H = inputs.params_struct
    @unpack ϕ, h = outputs
    nx, ny = size(ϕ)
    xc                = LinRange(-dx, Lx+dx, nx)         # including ghost points
    yc                = LinRange(-dy, Ly+dy, ny)
    x_plt = [xc[1]; xc .+ (xc[2]-xc[1])]
    y_plt = [yc[1]; yc .+ (yc[2]-yc[1])]
    ϕ[H .== 0.0] .= NaN
    h[H .== 0.0] .= NaN

    include("shmip_results.jl")
    xs     = LinRange(0, 100e3, nx-2)
    ϕ_ref = get_ϕ_ref(xs)
    h_ref = get_h_ref(xs)

    fig, axs = subplots(1, 2, figsize=(17,10))
    fig.subplots_adjust(hspace=0.3, wspace=0.5)
    lw = 2.5

    ind = size(ϕ ,2)÷2

    axs[1].plot(xc.*1e-3, ϕ[:, ind], label="GPU model"; lw)
    axs[1].plot(xs.*1e-3, ϕ_ref, label="GlaDS"; lw)
    axs[1].set_xlabel("x-coordinate (km)")
    axs[1].set_ylabel("ϕ (Pa)")

    axs[2].plot(xc.*1e-3, h[:, ind], label="GPU model"; lw)
    axs[2].plot(xs.*1e-3, h_ref, label="GlaDS"; lw)
    axs[2].legend()
    axs[2].set_xlabel("x-coordinate (km)")
    axs[2].set_ylabel("h (m)")

    fig.suptitle("cross-sections")
end

#plot_fieldresults(;inputs, outputs)


function plot_benchmarks()
    file = joinpath(@__DIR__, "benchmarks.jld2")
    benchmarks = load(file)
    glads = Dict("ode15s" => (res   = [2.64e3, 31.42e3, 59.9e3, 116.58e3, 314.5e3, 627.86e3],
                              t_tot = [0.087, 0.701, 1.359, 2.65, 6.586, 12.9780]),
                "ode113"  => (res   = [586, 7665, 11956, 19368, 48391, 90473],
                              t_tot = [0.087, 0.701, 1.359, 2.65, 6.586, 12.978])
    )
    commits = ["aba611ca46563ed1340675245febae8c714520ba"]#, "cceaa1b7d6913f66b04baaefaf492f91fd54bf04"]

    colors = Dict("achtzack01_GPU" => "deepskyblue",
                  "achtzack01_CPU-32threads" => "gold",
                  "achtzack01_CPU-16threads" => "purple",
                  "achtzack01_CPU-4threads" => "green",
                  "node35.octopoda_GPU" => "black")

    fig, axs = subplots(1, 2, figsize=(17,10))
    fig.subplots_adjust(hspace=0.3, wspace=0.5)

    for gitcommit in commits
        for unit in keys(benchmarks)
            if haskey(benchmarks[unit], gitcommit)
                (host, PU) = split(unit, "_")
                case = benchmarks[unit][gitcommit][:SHMIP_case]
                steadyst = benchmarks[unit][gitcommit][:iterations] .> 10^3
                res = benchmarks[unit][gitcommit][:nx] .* benchmarks[unit][gitcommit][:ny]
                t_tot = benchmarks[unit][gitcommit][:run_time]
                T_eff = benchmarks[unit][gitcommit][:T_eff]
                #nx = benchmarks[unit][gitcommit][:nx][steadyst]
                its = benchmarks[unit][gitcommit][:iterations][steadyst]

                if any(steadyst)
                    axs[1].plot(res[steadyst], t_tot[steadyst], "o", ms=5, color=colors[unit], label=unit) #marker="x", markersize=5,
                    axs[1].set_ylim(0, 400)
                    axs[1].set_xscale("log")
                    axs[1].set_xlabel("Grid size (nx*ny)")
                    axs[1].set_ylabel("Run time (s)")
                end
                axs[2].plot(res[.!steadyst], T_eff[.!steadyst], "-o", color=colors[unit], label=unit)
                axs[2].set_ylim(1, 500)
                axs[2].set_xscale("log")
                axs[2].set_yscale("log")
                axs[2].set_xlabel("Grid size (nx*ny)")
                axs[2].set_ylabel("T_eff (GB/s)")
                axs[2].legend()
            end
        end
    end
    for k in keys(glads)
        axs[1].plot(glads[k].res, glads[k].t_tot, "o", ms=5, color="cyan", label="GlaDS")
        axs[1].legend()
    end
end

#plot_benchmarks()