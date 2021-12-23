# note: let it run first with other parameters such that the values

using Printf, JLD2

include("../examples/SHMIP_cases.jl")


nx, ny = [(64, 32), (256, 128), (1024, 512)][2]
tol = Dict((64, 32)  => 5e-3,
           (128, 64) => 5e-2,
           (256, 128) => 1e-3,
           (1024, 512) => 5e-3)
γs = 0.5:0.05:0.7
dτ_hs = 3e-5:2e-6:5.4e-5

wall_time = zeros(length(γs), length(dτ_hs))
run_SHMIP(; nx, ny, γ_ϕ= 0.9, γ_h=0.9, dτ_h=7e-6, tol=1e-2, do_print=false); # "warming up" to make also the first run representative

for (col, γ) in enumerate(γs)
    for (row, dτ_h) in enumerate(dτ_hs)
        @printf("Running for γ = %1.2f, dτ_h = %1.2e, set %d out of %d. \n", γ, dτ_h, size(wall_time, 2)*(col-1) + row, length(wall_time))
        inputs, outputs = run_SHMIP(; nx, ny, γ_ϕ= γ, γ_h=γ, dτ_h, tol=tol[(nx,ny)], itMax=4*10^5, do_print=false, warmup=0);
        if any(isnan.([outputs.errs_h outputs.errs_ϕ]))
            time = NaN
            println("Unstable, only NaNs.")
        else
            time = outputs.time_tot
            @printf("Wall time = %1.2f seconds. \n", outputs.time_tot)
        end
        wall_time[col, row] = time
    end
end

jldsave("test/params_comp_nx" * string(nx) * "_ny" * string(ny) * ".jld2"; γs, dτ_hs, wall_time)
