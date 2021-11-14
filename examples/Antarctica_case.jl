using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters, GeoData, Infiltrator      # GeoData contains function NCDstack()
import Plots; Plt = Plots
# using PyPlot
const S = SheetModel

s_per_day  = 24 * 60 * 60
s_per_year = 365 * s_per_day

# all of continental Antarctica plus some 100km border:
Box(x::Tuple, y::Tuple) = (X(Between(x[1],x[2])), Y(Between(y[1],y[2])))
const box_antarctica = Box((-2750000, 2780000+1), (-2200000, 2300000+1))

datadir = joinpath(@__DIR__, "data_Antarctica")
!isdir(datadir) && mkdir(datadir)

# Note: To download the data the ISG shares need to be mounted (as described here https://vawiki.ethz.ch/vaw/informatics:samba_for_linux?s[]=samba)
filepath = homedir() * "/itet-stor/glazioarch/GlacioData/BedMachine_Antarctica/168596330/BedMachineAntarctica_2019-11-05_v01.nc"
url = "file://" * filepath
filename = if !isfile(joinpath(datadir, splitdir(url)[2]))
    print("Downloading data ...")
    try
        download(url, joinpath(datadir, splitdir(url)[2]))
    catch e
        if !isfile(filepath)
            print("Cannot find the file to be downloaded. Make sure ISG shares are mounted.")
        else
            print(e)
        end
    end
    println("done.")
else
    joinpath(datadir, splitdir(url)[2])
end

"""
Convert `missing` or `==missingval` to NaN and return a `GeoArray{Float32,2}`.
"""
function missing2nan(ar::A, T=Float32) where A<:AbstractGeoArray{<:Number,2}
    # convert to T first
    data = convert(Matrix{T}, ar.data)
    data[data.==ar.missingval] .= NaN
    return GeoArray(data; ar.dims, ar.name, ar.refdims, ar.metadata)
end
function missing2nan(ar::AbstractGeoArray, T=Float32) # Union{Missing,...}
    data = convert(Matrix{Union{T,Missing}}, copy(ar.data))
    data[ismissing.(data)] .= NaN
    data = convert(Matrix{T}, data)
    return GeoArray(data; ar.dims, ar.name, ar.refdims, ar.metadata)
end

## Read data; function taken from https://github.com/mauro3/FourDAntarcticaSubglacialRouting.jl
function read_bedmachine(thin=1)
    nc = NCDstack(datadir * "/BedMachineAntarctica_2019-11-05_v01.nc") # this doesn't do anything: , childkwargs=(crs=crs,))
    # read arrays into memory and thin them, if desired
    gas = []
    for k in keys(nc)
        if k in [:mapping, :surface,  :source, :geoid, :mask]
            continue # drop these fields
        end
        ga = (reverse(nc[k][1:thin:end, 1:thin:end], dims=2))[box_antarctica...] # this also loads it into memory
        ga = if k==:bed || k==:errbed || k==:surface || k==:firn
            # remove "missing" for bed, errbed and surface
            missing2nan(ga)
        else
            ga
        end
        push!(gas, ga)
    end
    mask = reverse(nc[:mask][1:thin:end, 1:thin:end], dims=2)

    # https://nsidc.org/data/nsidc-0756
    # mask==
    #       0 -> ocean, 1 -> ice-free land, 2 -> grounded ice, 3 -> floating ice, 4 -> Lake Vostok
    seamask  = (mask.==0) .| (mask.==3)  # treat floating ice as ocean
    icemask  = (mask.==2) .| (mask.==4)
    landmask = (mask.==1)
    push!(gas, GeoArray(seamask, name=:seamask)[box_antarctica...])
    push!(gas, GeoArray(icemask, name=:icemask)[box_antarctica...])
    push!(gas, GeoArray(landmask, name=:landmask)[box_antarctica...])


    # make dims a range:
    x, y = dims(gas[1])

    dx = x[2]-x[1]; @assert y[2]-y[1]==dx
    D = (X(x[1]:dx:x[end], mode=mode(x), metadata=metadata(x)),
         Y(y[1]:dx:y[end], mode=mode(y), metadata=metadata(y)))
    @assert val(D[1])==val(x) && val(D[2])==val(y)
    # the double-GeoStack is necessary as the first cannot use the dims-kw:
    return GeoStack(GeoStack(gas..., metadata=nc.metadata), dims=D), nc
#    return GeoStack(gas..., metadata=nc.metadata), nc
end
# topo, nc = read_bedmachine();
# to plot
# import Plots; Plt = Plots
# topo[:surface] |> Plt.plot

function run_Antarctica(;itMax=2*10^4, make_plot=false,
    dt=1e9, tsteps=1, γ_ϕ= 0.9, γ_h=0.8, dτ_ϕ_=1.0, dτ_h_=6e-6)      # parameters for pseudo-transient time stepping

    # read data
    topo, nc = read_bedmachine(50);

    # input parameters
    x, y              = Vector{Float64}.(Array.(dims(topo)))
    dx                = x[2] - x[1]
    dy                = y[2] - y[1]
    Lx                = x[end] - x[1]
    Ly                = y[end] - y[1]
    zb                = Matrix{Float64}(topo[:bed])
    H                 = Matrix{Float64}(topo[:thickness])
    calc_m(ix, iy, t) = 1e-6
    ttot = tsteps * dt

    ice_mask = Matrix{Bool}(topo[:icemask])
    bc_diric = Matrix{Bool}(ice_mask[2:end-1, 2:end-1] .& (topo[:seamask][1:end-2, 2:end-1] .| topo[:seamask][3:end, 2:end-1] .| topo[:seamask][2:end-1, 1:end-2] .| topo[:seamask][2:end-1, 3:end]))
    bc_no_xflux = Matrix{Bool}(abs.(diff(ice_mask .- topo[:landmask], dims=1)) .== 2)
    bc_no_yflux = Matrix{Bool}(abs.(diff(ice_mask .- topo[:landmask], dims=2)) .== 2)

    # to plot
    # import Plots; Plt = Plots
    # coords = findall(bc_diric .== 1)
    # xs     = [coords[i][1] for i in 1:length(coords)]
    # ys     = [coords[i][2] for i in 1:length(coords)]
    # Plots.scatter(xs, ys)

    # initial conditions
    ϕ_init = 1e6 * ones(size(H))
    h_init = 0.05 * ones(size(H))

    # call the SheetModel
    input = make_model_input(H, zb, Lx, Ly, dx, dy, ttot, dt, itMax, γ_ϕ, γ_h, dτ_ϕ_, dτ_h_, ϕ_init, h_init, calc_m, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux)
    output = runthemodel(;input...);

    # plotting
    @unpack N, ϕ, h, qx, qy,
            ittot, iters, Res_ϕ, Res_h, errs_ϕ, errs_h = output

    return input, output
end

input, output = run_Antarctica(γ_ϕ=0.8, γ_h=0.9, dτ_ϕ_=1e-8, dτ_h_=1e-14, itMax=10^6);
ice = input.ice_mask
ϕ = output.ϕ; ϕ[ice .== 0] .= NaN
h = output.h; h[ice .== 0] .= NaN
Plt.plot(Plt.heatmap(output.ϕ', aspect_ratio=:equal, title="ϕ"),
         Plt.heatmap(output.h', aspect_ratio=:equal, title="h"))