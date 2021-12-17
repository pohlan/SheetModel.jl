using JLD2, PyPlot, LaTeXStrings, Parameters, DataFrames

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 26

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
    file = "test/bm_results.jld2"
    bm = load(file)

    glads = Dict("ode15s" => (wt  = [2.64, 31.42, 59.9, 116.58, 314.5, 627.86],
                              dof = [0.087, 0.701, 1.359, 2.65, 6.586, 12.9780].*1e3),
                #"ode113"  => (t_tot = [586, 7665, 11956, 19368, 48391, 90473],
                #              res   = [0.087, 0.701, 1.359, 2.65, 6.586, 12.978].*1e3)
    )

    colors = Dict("achtzack01_GPU" => "deepskyblue",
                  "achtzack01_CPU_32threads" => "darkorange",
                  "achtzack01_CPU_16threads" => "darkorange",#"purple",
                  "achtzack01_CPU_8threads" => "darkorange",
                  "achtzack01_CPU_4threads" => "darkorange",#"green",
                  "node35.octopoda_GPU" => "black")

    function get_ls(kw)
        if endswith(kw, "32threads")
            style = "--"
        elseif endswith(kw, "16threads") || kw == "ode113"
            style = ":"
        elseif endswith(kw, "4threads")
            style = "-."
        elseif endswith(kw, "8threads")
            style = (0, (5, 10))
        else
            style = "-"
        end
        return style
    end

    function get_lab(unit)
        if startswith(unit, "node")
            lab = "GPU (Tesla V100)"
        elseif unit == "achtzack01_GPU"
            lab = "GPU (RTX 2070)"
        else
            lab = "CPU " * split(split(unit, "_")[3], "thr")[1] * " threads"
        end
    end

    lw = 3 # linewidth
    ms = 7 # markersize

    # going into steady state #
    # ------------------------#
    fig, axs = subplots(1, 3, figsize=(23,10), gridspec_kw=Dict(:width_ratios => [1,1,0.4], :wspace=>0.45, :hspace=>0.3))

    unit = "achtzack01_GPU"
    key  = "std-state"
    dof   = bm[unit][key].dof
    wt    = bm[unit][key].wall_time
    T_eff = bm[unit][key].T_eff
    nit   = bm[unit][key].nit

    axs[1].plot(dof, wt, "o", ls=get_ls(unit), ms=5, color=colors[unit], label=get_lab(unit); lw) #marker="x", markersize=5,
    #axs[1].set_ylim(0, 400)
    axs[1].set_xscale("log")
    axs[1].set_yscale("log")
    axs[1].set_xlabel(L"Degrees of freedom (nx $\cdot$ ny)", labelpad=10)
    axs[1].set_ylabel("Wall time (s)", labelpad=10)
    for k in keys(glads)
        axs[1].plot(glads[k].dof, glads[k].wt, marker="o", ls=get_ls(k), color="purple", label="GlaDS " * k; lw, ms)
        axs[1].legend()
    end

    # axs[2].plot(dof, nit, marker="o", ls=get_ls(unit), color=colors[unit], label=get_lab(unit); lw, ms)
    # #axs[2].set_ylim(1, 500)
    # axs[2].set_xscale("log")
    # axs[2].set_yscale("log")
    # axs[2].set_xlabel(L"Degrees of freedom (nx $\cdot$ ny)", labelpad=10)
    # axs[2].set_ylabel("Number of iterations", labelpad=10)

    #            T_eff        #
    # ------------------------#

    key = "Teff"
    for unit in ["achtzack01_GPU", "achtzack01_CPU_32threads", "achtzack01_CPU_16threads", "achtzack01_CPU_8threads", "achtzack01_CPU_4threads"] #, "node35.octopoda_GPU"]
        dof   = bm[unit][key].dof
        wt    = bm[unit][key].wall_time
        T_eff = bm[unit][key].T_eff
        nit   = bm[unit][key].nit

        axs[2].plot(dof, T_eff, "o", ls=get_ls(unit), color=colors[unit], label=get_lab(unit); lw, ms=10)
        axs[2].set_xscale("log")
        axs[2].set_yscale("log")
        axs[2].set_xlabel(L"Degrees of freedom (nx $\cdot$ ny)", labelpad=10)
        axs[2].set_ylabel(L"$T_{eff}$ (GB/s)", labelpad=10)
    end
    axs[2].hlines(900, dof[1], dof[end], color=colors["node35.octopoda_GPU"], ls=":", label="Tesla V100 bandwidth"; lw)
    axs[2].hlines(448, dof[1], dof[end], color=colors["achtzack01_GPU"], ls=":", label="RTX2070 bandwidth"; lw)
    axs[2].legend(handlelength=4, fontsize="small",bbox_to_anchor=(1.1,1), loc="upper left")

    axs[3].axis("off")
end

plot_benchmarks()