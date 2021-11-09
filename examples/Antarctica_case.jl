using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters, GeoData, Infiltrator      # GeoData contains function NCDstack()
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
topo, nc = read_bedmachine();
# to plot
# import Plots; Plt = Plots
# topo[:surface] |> Plt.plot

function run_Antarctica(;test_case, nx, ny, itMax=10^6, make_plot=false, printtime=10^5,
    dt=1e9, tsteps=1, γ_ϕ= 0.9, γ_h=0.8, dτ_ϕ_=1.0, dτ_h_=6e-6)      # parameters for pseudo-transient time stepping

    topo, nc = read_bedmachine();
    x, y = Array.(dims(topo))

    # physical domain (without ghost points)
    Lx                = x[end] - x[1]
    Ly                = y[end] - y[1]

    zb                = topo[:bed]
    H                 = topo[:thickness]
    calc_m(ix, iy, t) = 1e-6

    ttot = tsteps * dt

    # Initial condition
    ϕ_init, h_init = initial_conditions(
        xc,
        yc,
        calc_ϕ = (x, y) -> 100.0,
        #calc_ϕ = (x, y) -> 1e6/lx * x,
        #calc_ϕ = (x, y) -> exp(- 1e-2*(x-Lx/2)^2) * exp(-1e-2*(yc-Ly/2)^2),
        #calc_ϕ = (x, y) -> rand(),
        calc_h = (x, y) -> 0.04
    )
    input = make_model_input(H, zb, Lx, Ly, ttot, dt, itMax, γ_ϕ, γ_h, dτ_ϕ_, dτ_h_, ϕ_init, h_init, calc_m)
    output = runthemodel(;input...);
    @unpack N, ϕ, h, qx, qy,
            ittot, iters, Res_ϕ, Res_h, errs_ϕ, errs_h = output

    # plotting

    return (;input, SHMIP_case=test_case), output
end

# todo: masks for b.c.

