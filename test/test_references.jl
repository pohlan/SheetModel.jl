# This file writes or reads the reference test cases.
using JLD2

ref = ["write", "read"][2]

if ref == "write"
    include(joinpath(@__DIR__, "../examples/SHMIP_cases.jl"))
    ϕ_ref = Dict()
    h_ref = Dict()
    # A1 test case: sqrt geometry
    inputs, outputs = run_SHMIP(test_case="A1", nx=64, ny=32, dt = 1e9,
                                γ_ϕ= 0.8, γ_h=0.8, dτ_ϕ_=1.0, dτ_h_= 7e-6);
    ϕ_ref["A1"] =  outputs.ϕ[2:20:end-1, 2:10:end-1]
    h_ref["A1"] =  outputs.h[2:20:end-1, 2:10:end-1]

    # E1 test case: valley geometry
    inputs, outputs = run_SHMIP(test_case="E1", nx=256, ny=256, dt = 1e9,
                                γ_ϕ= 0.7, γ_h=0.5, dτ_ϕ_= 1.0, dτ_h_= 7e-4);
    ϕ_ref["E1"] =  outputs.ϕ[2:20:end-1, 2:10:end-1]
    h_ref["E1"] =  outputs.h[2:20:end-1, 2:10:end-1]

    # save the reference ϕ and h dictionaries in a file
    save(joinpath(@__DIR__, "phi_reference.jld2"), ϕ_ref)
    save(joinpath(@__DIR__, "h_reference.jld2"), h_ref)

elseif ref == "read"
    using JLD2, FileIO
    ϕ_ref = load(joinpath(@__DIR__, "phi_reference.jld2"))
    h_ref = load(joinpath(@__DIR__, "h_reference.jld2"))
end
