using JLD2, PyPlot, LaTeXStrings

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 24

"""
Plot a map of wall time vs damping (γ) and pseudo-time step (dτ_h)
"""
function plot_para_comp()

    γ, dτ_h, wt = load("test/params_comp_nx64_ny32.jld2","γs", "dτ_hs", "wall_time")
    wt_norm = wt ./ maximum(wt[.!isnan.(wt)])
    wt_min  = minimum(wt_norm[.!isnan.(wt_norm)])

    dγ       = γ[2]-γ[1]
    dτ       = dτ_h[2] - dτ_h[1]
    γ_plt    = [γ[1] - 0.5dγ; γ .+ 0.5dγ]
    dτ_h_plt = [dτ_h[1] - 0.5dτ; dτ_h .+ 0.5dτ]
    figure(figsize=(15, 10))
    pcolormesh(γ_plt, dτ_h_plt, wt_norm', cmap="viridis_r")
    xlabel(L"Damping parameter $γ$ = $γ_h$ = $γ_ϕ$", labelpad=20)
    ylabel(L"Pseudo-time step $dτ_h$ (s)", labelpad=20)
    cb = colorbar(cmap="viridis")
    cb.set_ticks([wt_min, 1])
    cb.set_ticklabels(["fast", "slow"])
    cb.set_label("Wall time, normalised")

    title("nx=64, ny=32")
end


plot_para_comp()

