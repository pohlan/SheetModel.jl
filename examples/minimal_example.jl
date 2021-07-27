using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using SheetModel, Parameters
const S = SheetModel

function run_example(;dt,
                      tsteps,               # number of timesteps
                      plot_output=true)     # whether to produce plots from model output
    input_params = S.Para(
            g     = 1.0,              # gravitational acceleration, m/s^2
            ρw    = 1.0,              # water density, kg/m^3
            ρi    = 0.91,             # ice density, kg/m^3
            α     = 1.25,             # first sheet flow exponent
            β     = 1.5,              # second sheet flow exponent
            k     = 1.0,              # sheet conductivity, m^(7/4)kg^(-1/2)
            n     = 3.0,              # Glen's flow law exponent
            A     = 1.0,              # ice flow constant, Pa^(-n)s^(-1)
            ev    = 0.0,              # englacial void ratio; SHMIP: 0 for ice-sheet, 1e-3 for valley glacier
            lr    = 1.0,              # horizontal cavity spacing, m
            hr    = 1.0,              # bedrock bump height, m
            ub    = 1.0,              # basal sliding speed, m/s

            xrange = (0.0, 1.0),  # domain length in x-direction, m
            yrange = (0.0, 1.5),  # domain length in y-direction, m
            nx = 64,
            ny = 64,

            calc_zs =  (x, y) -> 0.5 * sqrt(x+0.05),    # surface elevation, m (SHMIP)
            calc_zb = (x, y) -> 0.0,                    # bed elevation, m
            calc_m_xyt  = (x, y, t) -> 1.0,             # source term, m/s

            dt   = dt, #  TODO: Adaptive time stepping, in the end it shouldn't be specified as input
            ttot = dt * tsteps,

            itMax = 5*10^4,
            γ_ϕ  = 0.8,     # damping parameter for ϕ 0.6
            γ_h  = 0.7,     # damping parameter for h 0.8
            dτ_ϕ_ = 1.0,    # scaling factor for dτ_ϕ 1.0
            dτ_h_ = 1e-4,   # scaling factor for dτ_h 1e-7

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