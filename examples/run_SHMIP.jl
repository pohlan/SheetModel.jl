include("SHMIP_cases.jl")

nit, ϕ = run_SHMIP("A1", Nx=64, Ny=32,
                   γ_ϕ= 0.108907, γ_h=0.636720, dτ_ϕ_=799042.498250, dτ_h_= 25.472691, printtime=50, make_plot=true);

# run_SHMIP("F1", Nx=64, Ny=32,
#           γ_ϕ=0.017, γ_h=0.74, τ_ϕ_=8.6e5, τ_h_=31.0, printtime=50, make_plot=true)
