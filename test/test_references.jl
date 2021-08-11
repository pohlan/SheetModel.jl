# This file writes or reads the reference test cases.
using JLD2

ref = ["write", "read"][2]

if ref == "write"
    include(joinpath(@__DIR__, "../examples/SHMIP_cases.jl"))
    ϕ_ref = Dict()
    h_ref = Dict()
    for test_case in ["A1", "A2", "A3", "A4", "A5", "A6",
                      "D1", "D2", "D3", "D4", "D5"]
        inputs, outputs = run_SHMIP(test_case, Nx=64, Ny=32, dt = 5e7, tsteps=1,
                                   γ_ϕ= 0.8, γ_h=0.9, dτ_ϕ_=1.0, dτ_h_= 7e-6,
                                   printtime=1, printit=1000, make_plot=false);
        ϕ_ref[test_case] =  outputs.ϕ[1:20:end, 1:10:end]
        h_ref[test_case] =  outputs.h[1:20:end, 1:10:end]
    end

    # save the reference ϕ and h dictionaries in a file
    save(joinpath(@__DIR__, "phi_reference.jld2"), ϕ_ref)
    save(joinpath(@__DIR__, "h_reference.jld2"), h_ref)

elseif ref == "read"
    using JLD2, FileIO
    ϕ_ref = load(joinpath(@__DIR__, "phi_reference.jld2"))
    h_ref = load(joinpath(@__DIR__, "h_reference.jld2"))
end
