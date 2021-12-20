using JLD2, PyPlot, LaTeXStrings, Parameters, DataFrames, NetCDF, Statistics, SmoothingSplines

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 27

L2D = PyPlot.matplotlib.lines.Line2D

"""
Plot a map of wall time vs damping (γ) and pseudo-time step (dτ_h)
"""
function plot_para_comp()
    fig, axs = subplots(1, 3, figsize=(23,10))
    fig.subplots_adjust(hspace=0.3, wspace=0.5)
    p = (nothing for i=1:100)
    for (n, dof) in enumerate(["nx64_ny32", "nx128_ny64"]) #, "nx256_ny128"])
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

# function plot_error()
#     file = "test/bm_results.jld2"
#     bm = load(file)
#     unit = "achtzack01_GPU"
#     key  = "std-state"

#     errs_ϕ = bm[unit][key].errs_ϕ
#     errs_h = bm[unit][key].errs_h
#     iters  = bm[unit][key].iters
#     dof    = bm[unit][key].dof

#     lw = 2
#     fig, axs = subplots(1, 2, figsize=(17,10), gridspec_kw=Dict(:wspace=>0.45, :hspace=>0.3))

#     for n in 1:2:length(bm[unit][key].dof)
#         axs[1].plot(iters[n], errs_ϕ[n], label=string(n); lw) # color="gold"
#         axs[2].plot(iters[n], errs_h[n], label=string(n); lw) # , color="deepskyblue"
#     end

#     axs[1].set_xscale("log")
#     axs[1].set_yscale("log")
#     axs[1].set_xlabel("Iteration count")
#     axs[1].set_ylabel(L"RMS($f_\phi$)")
#     axs[1].legend()

#     axs[2].set_xscale("log")
#     axs[2].set_yscale("log")
#     axs[2].set_xlabel("Iteration count")
#     axs[2].set_ylabel(L"RMS($f_h$)")

# end

# after running run_SHMIP
#plot_error()


function shmip_results(case, model, x)
    path = "../shmip_results/" * model * "/"
    # h, ϕ and N
    h_shmip = ncread(path * case * "_" * model * ".nc", "h")
    H_shmip = ncread(path * case * "_" * model * ".nc", "H")
    N_shmip = ncread(path * case * "_" * model * ".nc", "N")
    ϕ_shmip = 910 * 9.81 .* H_shmip .- N_shmip
    coords1 = ncread(path * case * "_" * model * ".nc", "coords1")
    x1 = coords1[:, 1]
    y1 = coords1[:, 2]
    i_cross = abs.(mean(y1).-y1) .< 30

    h_cross = h_shmip[i_cross]
    ϕ_cross = ϕ_shmip[i_cross]

    # get an interpolated cross-section averaged over all y coordinates
    spl_h = fit(SmoothingSpline, x1[i_cross], vec(h_cross), 250.0) # λ=250.0
    spl_ϕ = fit(SmoothingSpline, x1[i_cross], vec(ϕ_cross), 250.0)
    ϕ_ref = predict(spl_ϕ, x)
    h_ref = predict(spl_h, x)
    return ϕ_ref, h_ref
end

function plot_error()
    inout = load("test/error_count.jld2", "inout")
    lw = 2.5
    figure(figsize=(13,10))
    for (n, ((inputs, outputs), color)) in enumerate(zip(inout, ["blue", "purple"]))
        @unpack ϕ, iters, errs_ϕ, errs_h = outputs
        nx, ny = size(ϕ)
        dof = nx * ny

        plot(iters, errs_ϕ, label=L"$f_\phi$, "*string(nx)*" x "*string(ny)*" grid", ls="-"; lw, color)
        plot(iters, errs_h, label=L"$f_h$, "*string(nx)*" x "*string(ny)*" grid", ls=":"; lw, color)
        #xscale("log")
        xticks([10^4, 10^5, 2*10^5, 3*10^5, 4*10^5], [L"$1\times10^4$", L"$1\times10^5$", L"$2\times10^5$", L"$3\times10^5$", L"$4\times10^5$"])
        yscale("log")
        xlabel("Iteration count", labelpad=10)
        ylabel(L"RMS($f$)", labelpad=10)
        legend()
    end
end

function plot_fieldresults()
    bm = load("test/bm_results.jld2")

    key = "A-suite"
    unit = "achtzack01_GPU"
    ϕs  = bm[unit][key].ϕ
    hs  = bm[unit][key].h
    tcs = bm[unit][key].test_case
    xs     = LinRange(0, 100e3, length(ϕs[1]))

    fig, axs = subplots(1, 2, figsize=(20,12), gridspec_kw=Dict(:wspace=>0.4, :hspace=>0.3))
    lw = 1.5
    lw_ref = 2

    custom_legend = [L2D([0], [0], color="black", lw=lw_ref)]

    for (n, (ϕ, h, tc)) in enumerate(zip(ϕs, hs, tcs))
        ϕ_ref, h_ref = shmip_results(tc, "jd", xs)

        axs[1].plot(xs.*1e-3, ϕ, zorder=2; lw)
        axs[1].plot(xs.*1e-3, ϕ_ref, color="black", zorder=1; lw=lw_ref)
        axs[1].set_xlabel("x-coordinate (km)", labelpad=10)
        axs[1].set_ylabel("ϕ (Pa)", labelpad=10)
        axs[1].legend(custom_legend, ["Reference ('jd')"], fontsize="small")
        axs[1].text(-0.2, 1., L"\bf{a}", transform=axs[1].transAxes, ha="right", va="top", fontsize=27)

        axs[2].plot(xs.*1e-3, h, label=tc, zorder=2; lw)
        axs[2].plot(xs.*1e-3, h_ref, color="black", zorder=1; lw=lw_ref)
        axs[2].set_xlabel("x-coordinate (km)", labelpad=10)
        axs[2].set_ylabel("h (m)", labelpad=10)
        axs[2].legend(title="GPU sheet model", fontsize="small")
        axs[2].text(-0.2, 1., L"\bf{b}", transform=axs[2].transAxes, ha="right", va="top", fontsize=27)
    end
end


function plot_benchmarks()
    bm = load("test/bm_results.jld2")

    glads = Dict("ode15s" => (wt  = [2.64, 31.42, 59.9, 116.58, 314.5, 627.86],
                              dof = [0.087, 0.701, 1.359, 2.65, 6.586, 12.9780].*1e3),
                #"ode113"  => (wt  = [586, 7665, 11956, 19368, 48391, 90473],
                #              dof = [0.087, 0.701, 1.359, 2.65, 6.586, 12.978].*1e3)
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
    fig, axs = subplots(1, 3, figsize=(26,10), gridspec_kw=Dict(:width_ratios => [1,1,0.4], :wspace=>0.5, :hspace=>0.3))

    key  = "std-state"

    for unit in ["node35.octopoda_GPU", "achtzack01_GPU"]
        dof   = bm[unit][key].dof
        wt    = bm[unit][key].wall_time
        T_eff = bm[unit][key].T_eff
        nit   = bm[unit][key].nit

        axs[1].plot(dof, wt, "o", ls=get_ls(unit), color=colors[unit], label=get_lab(unit); lw, ms) #marker="x", markersize=5,
        #axs[1].set_ylim(0, 400)
        axs[1].set_xscale("log")
        axs[1].set_yscale("log")
        axs[1].set_xlabel(L"Degrees of freedom (nx $\cdot$ ny)", labelpad=10)
        axs[1].set_ylabel("Wall time (s)", labelpad=10)
    end
    for k in keys(glads)
        axs[1].plot(glads[k].dof, glads[k].wt, marker="o", ls=get_ls(k), color="purple", label="GlaDS " * k; lw, ms)
        axs[1].legend()
    end
    axs[1].text(-0.2, 1., L"\bf{a}", transform=axs[1].transAxes, ha="right", va="top", fontsize=27)

    #            T_eff        #
    # ------------------------#

    key = "Teff"
    for unit in ["node35.octopoda_GPU", "achtzack01_GPU", "achtzack01_CPU_32threads", "achtzack01_CPU_16threads", "achtzack01_CPU_8threads", "achtzack01_CPU_4threads"] #, ]
        dof   = bm[unit][key].dof
        wt    = bm[unit][key].wall_time
        T_eff = bm[unit][key].T_eff
        nit   = bm[unit][key].nit

        axs[2].plot(dof, T_eff, "o", ls=get_ls(unit), color=colors[unit], label=get_lab(unit); lw, ms)
        axs[2].set_xscale("log")
        axs[2].set_yscale("log")
        axs[2].set_xlabel(L"Degrees of freedom (nx $\cdot$ ny)", labelpad=10)
        axs[2].set_ylabel(L"$T_{eff}$ (GB/s)", labelpad=10)
    end
    dof = bm["achtzack01_GPU"][key].dof
    axs[2].hlines(900, dof[1], dof[end], color=colors["node35.octopoda_GPU"], ls=":", label="Tesla V100 bandwidth"; lw)
    axs[2].hlines(448, dof[1], dof[end], color=colors["achtzack01_GPU"], ls=":", label="RTX2070 bandwidth"; lw)
    axs[2].legend(handlelength=4, fontsize="small",bbox_to_anchor=(1.1,1), loc="upper left")
    axs[2].text(-0.2, 1., L"\bf{b}", transform=axs[2].transAxes, ha="right", va="top", fontsize=27)

    axs[3].axis("off") # trick to make legend fit in the figure
end
