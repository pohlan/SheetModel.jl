# -------------------------------------------------------------------------#
# This script runs the model directly with dimensionless numbers as input  #
# Without scaling and descaling as in SHMIP_cases.jl                       #
# -------------------------------------------------------------------------#

using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters, ParallelStencil # for ProfileView, gtk needs to be installed
const S = SheetModel

function run_example(;dt,
                      nx, ny,
                      tsteps)         # number of timesteps

    input_params = S.model_input(
            g     = 1.0,              # gravitational acceleration
            ρw    = 1.0,              # water density
            ρi    = 0.91,             # ice density
            α     = 1.25,             # first sheet flow exponent
            β     = 1.5,              # second sheet flow exponent
            k     = 1.0,              # sheet conductivity
            n     = 3.0,              # Glen's flow law exponent
            A     = 1.0,              # ice flow constant
            ev    = 0.0,              # englacial void ratio; SHMIP: 0 for ice-sheet, 1e-3 for valley glacier
            lr    = 1.0,              # horizontal cavity spacing
            hr    = 1.0,              # bedrock bump height
            ub    = 1.0,              # basal sliding speed

            xrange = (0.0, 1.0),  # domain length in x-direction
            yrange = (0.0, 1.5),  # domain length in y-direction
            nx = nx,
            ny = ny,

            calc_zs =  (x, y) -> 0.5 * sqrt(x+0.05),    # surface elevation (SHMIP)
            calc_zb = (x, y) -> 0.0,                    # bed elevation
            calc_m_phys  = (x, y, t) -> 1.0,             # source term

            dt   = dt,                  # time step
            ttot = dt * tsteps,         # total time

            itMax = 1,
            γ_ϕ  = 0.6,     # damping parameter for ϕ
            γ_h  = 0.8,     # damping parameter for h
            dτ_ϕ_ = 1.0,    # scaling factor for dτ_ϕ
            dτ_h = 3e-3,   # scaling factor for dτ_h

            # Dimensionless numbers
            Ψ   = 1.0,
            Σ   = 2.0,
            Γ   = 6e3,
            Λ   = 3e-3
        )

    calc_Λ_m! = @parallel_indices (ix,iy)   function calc_Λ_m!(Λ_m, Λ, t)
                                                if (ix <= size(Λ_m, 1) && iy <= size(Λ_m, 2))
                                                    Λ_m[ix, iy] = Λ * 1
                                                end
                                                return
                                            end

    ϕ_init = 0.5 * ones(nx, ny)
    h_init = 0.5 * ones(nx, ny)

    ice_mask = ϕ_init .> 0.
    bc_diric    = falses(nx, ny); bc_diric[2, :] .= true
    bc_no_xflux = diff(ice_mask, dims=1) .!= 0
    bc_no_yflux = diff(ice_mask, dims=2) .!= 0

    #ProfileView.@profview S.runthemodel_scaled(input_params, ϕ_init, h_init, S.CuParams(nx=Nx, ny=Ny), 1000) # for profiling
    # In VSCode, ProfileView.xx is necessary (https://github.com/timholy/ProfileView.jl/pull/172/commits/5ea809fe6409a41b96cfba0800b78d708f1ad604)

    inputs, outputs = S.runthemodel_scaled(scaled_params, ϕ_init, h_init, calc_Λ_m!, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux, do_print=true, warmup=0)

    # plot output
    #@unpack xc, yc, H = params_struct
    #@unpack N, ϕ, h, qx, qy,
    #        ittot, iters, Err_ϕ, Err_h, Res_ϕ, Res_h,
    #        errs_ϕ, errs_h, errs_ϕ_rel, errs_h_rel,
    #        errs_ϕ_res, errs_h_res, errs_ϕ_resrel, errs_h_resrel = output
    #if plotting
    #    S.plot_output(xc, yc, H, N, h, qx, qy, Err_ϕ, Err_h,
    #                  iters, errs_h, errs_ϕ, errs_ϕ_rel, errs_h_rel,
    #                  errs_ϕ_res, errs_h_res, errs_ϕ_resrel, errs_h_resrel)
    #end
end

run_example(dt=1e8, tsteps=1, nx=1024, ny=1024)
#run_example(dt=2e7, tsteps=5)
