using NetCDF, PyPlot, Statistics, SmoothingSplines

case = "A1"
#cd(@__DIR__)
# ncinfo(case * "_mw.nc") # to see which variables are available

# h, ϕ and N
print(pwd())
h_shmip = ncread("test/" * case * "_mw.nc", "h")
H_shmip = ncread("test/" * case * "_mw.nc", "H")
N_shmip = ncread("test/" * case * "_mw.nc", "N")
ϕ_shmip = 910 * 9.81 .* H_shmip .- N_shmip
coords1 = ncread("test/" * case * "_mw.nc", "coords1")
x1 = coords1[:, 1]
y1 = coords1[:, 2]

# get an interpolated cross-section averaged over all y coordinates
spl_h = fit(SmoothingSpline, x1, vec(h_shmip), 250.0) # λ=250.0
spl_ϕ = fit(SmoothingSpline, x1, vec(ϕ_shmip), 250.0)
get_ϕ_ref(x) = predict(spl_ϕ, x)
get_h_ref(x) = predict(spl_h, x)


# plot reference fields
function plot_GlaDs_result(x1, y1, h_shmip, spl_h, ϕ_shmip, spl_ϕ)
    x_crsec = 0:1e3:maximum(x1)
    h_crsec = predict(spl_h, x_crsec)
    ϕ_crsec = predict(spl_ϕ, x_crsec)

    figure()
    subplot(2, 2, 1)
    scatter(x1, y1, 20, h_shmip)
    colorbar()
    title("h")
    subplot(2, 2, 3)
    plot(x_crsec, h_crsec)
    # plot(x1[isapprox.(y1, mean(y1), atol=500)], h_shmip[isapprox.(y1, mean(y1), atol=500)])  # to not average over the entire y dimension
    title("h at y ≈ " * string(round(Int, mean(y1))) * " m")

    subplot(2, 2, 2)
    scatter(x1, y1, 20, ϕ_shmip)
    colorbar()
    title("ϕ")
    subplot(2, 2, 4)
    plot(x_crsec, ϕ_crsec)
    title("ϕ at y ≈ " * string(round(Int, mean(y1))) * " m")
end

# plot_GlaDs_result(x1, y1, h_shmip, h_crsec, ϕ_shmip, ϕ_crsec)
