using JLD2, PyPlot, LaTeXStrings

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 24

"""
Plot a map of wall time vs damping (γ) and pseudo-time step (dτ_h)
"""
function plot_para_comp()
    fig, axs = subplots(1, 2, figsize=(23,10))
    fig.subplots_adjust(hspace=0.3, wspace=0.5)
    for (n, dof) in enumerate(["nx64_ny32", "nx128_ny64"])
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


plot_para_comp()
