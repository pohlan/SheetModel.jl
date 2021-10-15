# -------------------------------------------------------------------------#
# This script runs the model directly with dimensionless numbers as input  #
# Without scaling and descaling as in SHMIP_cases.jl                       #
# -------------------------------------------------------------------------#

using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters, ProfileView, ParallelStencil # for ProfileView, gtk needs to be installed
const S = SheetModel

function run_example(;dt,
                      Nx, Ny,
                      tsteps,               # number of timesteps
                      plotting=true)     # whether to produce plots from model output
    input_params = S.Para(
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
            nx = Nx,
            ny = Ny,

            calc_zs =  (x, y) -> 0.5 * sqrt(x+0.05),    # surface elevation (SHMIP)
            calc_zb = (x, y) -> 0.0,                    # bed elevation
            calc_m_xyt  = (x, y, t) -> 1.0,             # source term

            dt   = dt,                  # time step
            ttot = dt * tsteps,         # total time

            itMax = 1,
            γ_ϕ  = 0.6,     # damping parameter for ϕ
            γ_h  = 0.8,     # damping parameter for h
            dτ_ϕ_ = 1.0,    # scaling factor for dτ_ϕ
            dτ_h_ = 3e-3,   # scaling factor for dτ_h

            # Dimensionless numbers
            Σ   = 2.0,
            Γ   = 6e3,
            Λ   = 3e-3
        )

    ϕ_init = 0.5 * ones(Nx, Ny)
    h_init = 0.5 * ones(Nx, Ny)

    #ProfileView.@profview S.runthemodel_scaled(input_params, ϕ_init, h_init, S.CuParams(nx=Nx, ny=Ny), 1000) # for profiling
    # In VSCode, ProfileView.xx is necessary (https://github.com/timholy/ProfileView.jl/pull/172/commits/5ea809fe6409a41b96cfba0800b78d708f1ad604)

    output = S.runthemodel_scaled(input_params, ϕ_init, h_init, 1000, 1)

    # plot output
    #@unpack xc, yc, H = input_params
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

run_example(dt=1e8, tsteps=1, Nx=64, Ny=64, plotting=false)
#run_example(dt=2e7, tsteps=5, plotting=true)
