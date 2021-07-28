# -------------------------------------------------------------------------#
# This script runs the model directly with dimensionless numbers as input  #
# Without scaling and descaling as in SHMIP_cases.jl                       #
# -------------------------------------------------------------------------#

using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters
const S = SheetModel

function run_example(;dt,
                      tsteps,               # number of timesteps
                      plot_output=true)     # whether to produce plots from model output
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
            nx = 64,
            ny = 64,

            calc_zs =  (x, y) -> 0.5 * sqrt(x+0.05),    # surface elevation (SHMIP)
            calc_zb = (x, y) -> 0.0,                    # bed elevation
            calc_m_xyt  = (x, y, t) -> 1.0,             # source term

            dt   = dt,                  # time step
            ttot = dt * tsteps,         # total time

            itMax = 5*10^4,
            γ_ϕ  = 0.6,     # damping parameter for ϕ
            γ_h  = 0.7,     # damping parameter for h
            dτ_ϕ_ = 1.0,    # scaling factor for dτ_ϕ
            dτ_h_ = 1e-7,   # scaling factor for dτ_h

            # Dimensionless numbers
            Σ   = 2.0,
            Γ   = 6e3,
            Λ   = 3e-3
        )

        @unpack nx, ny = input_params
        ϕ0 = 0.5 * ones(nx, ny)
        h0 = 0.5 * ones(nx, ny)

    output = S.runthemodel_scaled(input_params, ϕ0, h0, 100, 1) # output is a struct (see model_output in SheetModel.jl) containing all the variable and error fields as well as number of iterations

    # plot output
    @unpack xc, yc, H = input_params
    @unpack N, ϕ, h, qx, qy, ittot, Err_ϕ, Err_h, errs_ϕ, errs_h, errs_ϕ_res, errs_h_res, errs_ϕ_rel, errs_h_rel, Res_ϕ, Res_h, qx_ice, qy_ice = output
    if plot_output
        S.plot_output(xc, yc, H, N, h, qx, qy, qx_ice, qy_ice, Err_ϕ, Err_h, errs_h, errs_ϕ, errs_ϕ_res, errs_h_res, errs_ϕ_rel, errs_h_rel)
    end
end

run_example(dt=100, tsteps=1, plot_output=true)