# note: let it run first with other parameters such that the values

using Printf, JLD2

include("../examples/SHMIP_cases.jl")

nx, ny = 64, 32
γs = 0.8:0.01:0.9
dτ_hs = 5e-6:1e-6:2.5e-5

wall_time = zeros(length(γs), length(dτ_hs))
run_SHMIP(; nx, ny, γ_ϕ= 0.9, γ_h=0.9, dτ_h_=7e-6, tol=1e-3, do_print=false); # "warming up" to make also the first run representative

for (col, γ) in enumerate(γs)
    for (row, dτ_h_) in enumerate(dτ_hs)
        @printf("Running for γ = %1.2f, dτ_h_ = %1.2e, set %d out of %d. \n", γ, dτ_h_, size(wall_time, 2)*(col-1) + row, length(wall_time))
        inputs, outputs = run_SHMIP(; nx, ny, γ_ϕ= γ, γ_h=γ, dτ_h_, tol=1e-3, itMax=25*10^4, do_print=false, warmup=0);
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

save("params_comp_nx" * string(nx) * "_ny" * string(ny) * ".jld2", "wall_time", wall_time)
