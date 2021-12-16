using Printf, Infiltrator

# used grids:
# - normal grid (i,j), size (nx, ny)
#   ϕ, h, Res_ϕ, Res_h, d_eff, ...
# - qx grid (m,n) staggered in x-dir.: (nx-1, ny)
# - qy grid (m,n) staggered in y-dir.: (nx, ny-1)


### MACROS ###

"""
Calculate water pressure; input coordinates on ϕ/h grid.
Needs access to ϕ, zb
"""
macro pw(ix, iy) esc(:(ϕ[$ix, $iy] - zb[$ix, $iy])) end

"""
Calculate effective pressure; input coordinates on ϕ/h grid.
Needs access to ϕ, H, zb
"""
macro N(ix, iy) esc(:(H[$ix, $iy] - @pw($ix, $iy))) end

"""
Calculate closure rate; input coordinates on ϕ/h grid.
Needs access to ϕ, H, zb, n, h, Γ, A
"""
macro vc(ix, iy) esc(:(A * h[$ix, $iy] * abs(@N($ix, $iy))^(n-1) * @N($ix, $iy))) end

"""
Calculate opening rate; input coordinates on ϕ/h grid.
Needs access to h, hr, Σ, ub, lr
"""
macro vo(ix, iy) esc(:(h[$ix, $iy] < hr ? ub / lr  * (hr - h[$ix, $iy]) : 0.0)) end

"""
Calculate hydraulic gradient in x-direction; input coordinates on qx grid.
Implicitly sets zero-flux boundary conditions.
Needs access to H, ϕ, dx_
"""
macro dϕ_dx(ix, iy) esc(:( bc_no_xflux[$ix, $iy] == 1 ? 0. : (ϕ[$ix+1, $iy] - ϕ[$ix, $iy]) * dx_ )) end

"""
Calculate hydraulic gradient in y-direction; input coordinates on qy grid.
Implicitly sets zero-flux boundary conditions.
Needs access to H, ϕ, dy_
"""
macro dϕ_dy(ix, iy) esc(:( bc_no_yflux[$ix, $iy] == 1 ? 0. : (ϕ[$ix, $iy+1] - ϕ[$ix, $iy]) * dy_ )) end

"""
Calculate absolute hydraulic gradient, |∇ϕ|;
input coordinates on ϕ/h grid (nx, ny), but only possible to calculate on inner points.
Needs access to ϕ, dx_, dy_, H
"""
macro gradϕ(ix, iy) esc(:( sqrt(
                                  (0.5 * (@dϕ_dx($ix, $iy) + @dϕ_dx($ix-1, $iy)))^2
                                + (0.5 * (@dϕ_dy($ix, $iy) + @dϕ_dy($ix, $iy-1)))^2
                                  ))) end
"""
Calculate effective diffusivity;
input coordinates on ϕ/h grid (nx, ny), but only possible to calculate on inner points.
Needs access to k, h, α, β, small, ϕ, dx_, dy_, H
"""
macro d_eff(ix, iy) esc(:( k * h[$ix, $iy]^α * (@gradϕ($ix, $iy) + small)^(β-2) )) end

"""
Calculate pseudo-time step of ϕ,
input coordinates on ϕ/h grid (nx, ny), but only possible to calculate on inner points.
Needs access to dτ_ϕ_, dx_, dy_, min_dxy2, dt, k, h, α, β, small, ϕ, H
"""
macro dτ_ϕ(ix, iy) esc(:( dτ_ϕ_ * min(min_dxy2 / d_eff[$ix, $iy] / 4.1, dt))) end
#macro dτ_ϕ(ix, iy) esc(:( dτ_ϕ_ *  (min_dxy2_ / d_eff[$ix, $iy] / 4.1) .+ 1.0 * dt_) .^(-1))) end # other definitions...

"""
Calculate flux in x-direction using an upstream scheme; input coordinates on qx grid.
Needs access to k, h, α, β, small, ϕ, dx_, dy_, H
"""
macro qx(ix, iy) esc(:(
    - d_eff[$ix+1, $iy] * max(@dϕ_dx($ix, $iy), 0.) +   # flux in negative x-direction
    - d_eff[$ix,   $iy] * min(@dϕ_dx($ix, $iy), 0.)    # flux in positive x-direction
    )) end
"""
Calculate flux in y-direction using an upstream scheme; input coordinates on qy grid.
Needs access to k, h, α, β, small, ϕ, dx_, dy_, H
"""
macro qy(ix, iy) esc(:(
    - d_eff[$ix, $iy+1] * max(@dϕ_dy($ix, $iy), 0.) +   # flux in negative y-direction
    - d_eff[$ix, $iy  ] * min(@dϕ_dy($ix, $iy), 0.)    # flux in positive y-direction
    )) end


"""
Calculate residual of ϕ; input coordinates on ϕ/h grid (but only defined on inner grid points 2:nx-1, 2:ny-1, due to d_eff).
Needs access to ϕ, ϕ_old, h, H, k, α, β, small, dx_, dy_, min_dxy2, dt, dt_, Λ_m, hr, zb, n, Ψ, Σ, ub, lr
"""
macro Res_ϕ(ix, iy) esc(:(  ice_mask[$ix, $iy] == 1 ? (                                                      # only calculate at points with non-zero ice thickness
                                                # - Ψ * (ϕ[$ix, $iy] - ϕ_old[$ix, $iy]) * dt_                                                                # dϕ/dt
                                                 - ( (@qx($ix, $iy) - @qx($ix-1, $iy)) * dx_ + (@qy($ix, $iy) - @qy($ix, $iy-1)) * dy_ )    # divergence
                                                 - (Σ * @vo($ix, $iy) - Γ * @vc($ix, $iy))                                                                          # dh/dt
                                                 + Λ_m[$ix, $iy]                                                                                            # source term Λ * m
                                                ) : 0.0
)) end

"""
Calculate residual of h; input coordinates on ϕ/h grid.
Needs access to ϕ, h, h_old, H, dt_, ub, lr, hr, zb, n, Σ, Γ
"""
macro Res_h(ix, iy) esc(:(  ice_mask[$ix, $iy] == 1 ? (
                                                # - (h[$ix, $iy] - h_old[$ix, $iy]) * dt_ +
                                                 (Σ * @vo($ix, $iy) - Γ * @vc($ix, $iy))
                                                ) : 0.0
)) end

### KERNEL FUNCTIONS ###
@parallel_indices (ix, iy) function update_deff!(d_eff, ϕ, h, dx_, dy_, k, α, β, dt, Σ, Γ, ub, lr, hr, n, H, zb, small, bc_no_xflux, bc_no_yflux,)
    nx, ny = size(ϕ)
    if (1 < ix < nx && 1 < iy < ny)
        d_eff[ix, iy] = @d_eff(ix, iy)
    end
    return
end

"""
Calculate residuals of ϕ and h and store them in arrays.
Used for error calculation and only to be carried out every xx iterations, e.g. every thousand.
"""
@parallel_indices (ix,iy) function residuals!(ϕ, ϕ_old, h, h_old, Res_ϕ, Res_h, Λ_m, d_eff, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux,
                                              dx_, dy_, min_dxy2, k, α, β, dt, dt_, ub, lr, hr, Ψ, Σ, Γ, A, n, H, zb, small)
    nx, ny = size(ϕ)
    if (1 < ix < nx && 1 < iy < ny)
        # residual of ϕ
        if  bc_diric[ix, iy] == 1    # position where dirichlet b.c. are imposed
            Res_ϕ[ix-1, iy-1] = 0.
        else
            Res_ϕ[ix-1, iy-1] = @Res_ϕ(ix, iy)
        end

        # residual of h
        Res_h[ix-1, iy-1] = @Res_h(ix, iy)
    end
    return
end

@parallel_indices (ix, iy) function update_h_only!(ϕ, ϕ2, ϕ_old, h, h2, h_old, Λ_m, d_eff, iter, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux,
                                                   dx_, dy_, min_dxy2, k, α, β, dt, dt_, ub, lr, hr, Ψ, Σ, Γ, A, n, H, zb, small,
                                                   dϕ_dτ, dh_dτ, γ_ϕ, γ_h, dτ_h, dτ_ϕ_)
    nx, ny = size(ϕ)
    if (1 < ix < nx && 1 < iy < ny)
        dh_dτ[ix, iy] = @Res_h(ix, iy) + γ_h * dh_dτ[ix, iy]
        h2[ix, iy] = h[ix, iy] + dτ_h * dh_dτ[ix, iy]
    end
    return
end

"""
Update the fields of ϕ and h using the pseudo-transient method with damping.
"""
@parallel_indices (ix,iy) function update_fields!(ϕ, ϕ2, ϕ_old, h, h2, h_old, Λ_m, d_eff, iter, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux,
                                                  dx_, dy_, min_dxy2, k, α, β, dt, dt_, ub, lr, hr, Ψ, Σ, Γ, A, n, H, zb, small,
                                                  dϕ_dτ, dh_dτ, γ_ϕ, γ_h, dτ_h, dτ_ϕ_)
    nx, ny = size(ϕ)
    if (1 < ix < nx && 1 < iy < ny)
        # update ϕ
        if bc_diric[ix, iy] == 1
            ϕ2[ix, iy] = 0.             # dirichlet boundary conditions
        else
            dϕ_dτ[ix, iy] = @Res_ϕ(ix, iy) + γ_ϕ * dϕ_dτ[ix, iy]
            ϕ2[ix, iy] = ϕ[ix, iy] + @dτ_ϕ(ix, iy) * dϕ_dτ[ix, iy]
        end

        # update h
        dh_dτ[ix, iy] = @Res_h(ix, iy) + γ_h * dh_dτ[ix, iy]
        h2[ix, iy] = max(h[ix, iy] + dτ_h * dh_dτ[ix, iy], 0.0)
    end
    return
end

"""
Apply Dirichlet boundary conditions to ϕ=0.
Neumann boundary conditions are applied when calculating the hydraulic gradient (@dϕ_dx and @dϕ_dy macros).
"""
@parallel_indices (ix,iy) function apply_bc!(ϕ, bc_diric)
    nx, ny = size(ϕ)
    if (ix <= nx && iy <= ny) && bc_diric[ix, iy] == 1
        ϕ[ix, iy] = 0.
    end
    return
end

"""
Calculate the difference between previous and updated ϕ and h, used for error calculation.
"""
@parallel_indices (ix,iy) function update_difference!(Δϕ, ϕ, ϕ2, Δh, h, h2)
    if (ix <= size(ϕ, 1) && iy <= size(ϕ, 2))
        Δϕ[ix, iy] = abs(ϕ[ix, iy] - ϕ2[ix, iy])
        Δh[ix, iy] = abs(h[ix, iy] - h2[ix, iy])
    end
    return
end

"""
Assign updated ϕ and h values to ϕ_old and h_old, respectively.
Carried out after each physical time step.
"""
@parallel_indices (ix,iy) function old2new!(ϕ, ϕ_old, h, h_old)
    if (ix <= size(ϕ, 1) && iy <= size(ϕ, 2))
        ϕ_old[ix, iy] = ϕ[ix, iy]
        h_old[ix, iy] = h[ix, iy]
    end
    return
end

"""
Calculate effective pressure N and fluxes (qx, qy) at the end of model run, for plotting.
"""
@parallel_indices (ix,iy) function output_params!(N, qx, qy, ϕ, h, d_eff, H, zb, k, α, β, dx_, dy_, small, bc_no_xflux, bc_no_yflux)
    nx, ny = size(ϕ)
    if (ix <= nx && iy <= ny)
        N[ix, iy]  = @N(ix, iy)
        if (1 < ix < nx && 1 < iy < ny)
            qx[ix, iy] = @qx(ix, iy)
            qy[ix, iy] = @qy(ix, iy)
        end
    end
    return
end

"""
Run the model with scaled parameters.
"""
@views function runthemodel_scaled(params, ϕ_init, h_init, calc_Λ_m!, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux, update_h_only, do_print, warmup)
    @unpack ev, n, A, Ψ, Σ, Γ, Λ, dx, dy, k, α, β,
            H, zb, ub, hr, lr, dt, ttot, tol, itMax, γ_ϕ, γ_h, dτ_ϕ_, dτ_h_ = params

    # Pre-calculate reciprocals for better performance
    min_dxy2 = min(dx, dy)^2
    dx_, dy_, dt_, min_dxy2_ = (dx, dy, dt, min_dxy2).^(-1)

    # Array allocation
    nx, ny = size(ϕ_init)
    qx, qy, d_eff, Λ_m, N, dϕ_dτ, dh_dτ, Res_ϕ, Res_h = array_allocation(nx, ny)

    # Apply boundary conditions
    @parallel apply_bc!(ϕ_init, bc_diric)

    ϕ_old   = copy(ϕ_init)
    ϕ       = copy(ϕ_init)
    ϕ2      = copy(ϕ_init)
    h_old   = copy(h_init)
    h       = copy(h_init)
    h2      = copy(h_init)

    # for iterations vs. error plot
    iters      = Int64[]
    errs_ϕ     = Float64[]
    errs_h     = Float64[]
    err_ϕ      = 1e10
    err_h      = 1e10

    # initiate time loop parameters
    t = 0.0; tstep=0; ittot = 0; t_tic = 0.

    # Physical time loop
    while t<ttot
        iter = 0

        @parallel calc_Λ_m!(Λ_m, Λ, t)

        # Pseudo-transient iteration
        while !(max(err_ϕ, err_h) < tol) && iter<itMax && !any(isnan.([err_ϕ, err_h]))

            if  err_h > 1e-3 && update_h_only # once update_h_only = false it cannot go back
                dτ_h = 1e-3
            else
                dτ_h = dτ_h_
                update_h_only = false
            end

            # don't consider first ten iterations for performance measure
            if (iter == warmup) t_tic = Base.time() end

            if update_h_only
                @parallel update_h_only!(ϕ, ϕ2, ϕ_old, h, h2, h_old, Λ_m, d_eff, iter, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux,
                                         dx_, dy_, min_dxy2, k, α, β, dt, dt_, ub, lr, hr, Ψ, Σ, Γ, A, n, H, zb, small,
                                         dϕ_dτ, dh_dτ, γ_ϕ, γ_h, dτ_h, dτ_ϕ_)
            else
                # update ϕ and h
                @parallel update_deff!(d_eff, ϕ, h, dx_, dy_, k, α, β, dt, Σ, Γ, ub, lr, hr, n, H, zb, small, bc_no_xflux, bc_no_yflux)
                @parallel update_fields!(ϕ, ϕ2, ϕ_old, h, h2, h_old, Λ_m, d_eff, iter, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux,
                                         dx_, dy_, min_dxy2, k, α, β, dt, dt_, ub, lr, hr, Ψ, Σ, Γ, A, n, H, zb, small,
                                         dϕ_dτ, dh_dτ, γ_ϕ, γ_h, dτ_h, dτ_ϕ_)
            end

            # pointer swap
            ϕ, ϕ2 = ϕ2, ϕ
            h, h2 = h2, h

            iter += 1

            # determine the errors
            if (update_h_only || iter % 1000 == 0) && itMax > 10^4     # only calculate errors every 1000 time steps
                                                                       # and if itMax is high, i.e. if convergence is the goal
                # update the residual arrays
                @parallel residuals!(ϕ, ϕ_old, h, h_old, Res_ϕ, Res_h, Λ_m, d_eff, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux,
                                     dx_, dy_, min_dxy2, k, α, β, dt, dt_, ub, lr, hr, Ψ, Σ, Γ, A, n, H, zb, small)

                # residual error
                err_ϕ = norm(Res_ϕ) / sqrt(length(Res_ϕ))
                err_h = norm(Res_h) / sqrt(length(Res_h))

                # save error evolution in vector
                append!(iters, iter)
                append!(errs_ϕ, err_ϕ)
                append!(errs_h, err_h)

                if do_print @printf("iterations = %d, error ϕ = %1.2e, error h = %1.2e \n", iter, err_ϕ, err_h) end

            end


        end
        ittot += iter; tstep += 1; t += dt

        @parallel old2new!(ϕ, ϕ_old, h, h_old)
    end

    # Perfomance measures
    t_toc = Base.time() - t_tic                # execution time, s
    A_eff = (5+10)/1e9*nx*ny*sizeof(Float64)   # effective main memory access per iteration [GB];
                                               # 5 write arrays (dϕ_dτ, ϕ2, dh_dτ, h2, d_eff)
                                               # 10 read arrays (ϕ, ϕ_old, dϕ_dτ, h, h_old, dh_dτ, H, zb, m, d_eff)
    t_it  = t_toc/(ittot-10)                   # execution time per iteration, s
    T_eff = A_eff/t_it                         # effective memory throughput, GB/s
    if do_print @printf("Time = %1.3f sec, T_eff = %1.1f GB/s, iterations total = %d, (nx, ny) = (%d, %d)\n", t_toc, round(T_eff, sigdigits=3), ittot, nx, ny) end

    # calculate N, qx and qy as output parameters
    @parallel output_params!(N, qx, qy, ϕ, h, d_eff, H, zb, k, α, β, dx_, dy_, small, bc_no_xflux, bc_no_yflux,)

    return model_output{Data.Array}(;   N, ϕ, h, qx, qy,
                            Res_ϕ, Res_h, errs_ϕ, errs_h, ittot, iters,
                            time_tot=t_toc, T_eff)
end

"""
Scale the parameters and call the model run function.
"""
@views function runthemodel(;params_struct::model_input, ϕ_init, h_init, calc_m, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux, update_h_only, do_print=true, warmup=10)
    scaled_params, ϕ_init, h_init, calc_m, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux, ϕ_, N_, h_, q_ = scaling(params_struct, ϕ_init, h_init, calc_m, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux)
    calc_Λ_m! = @parallel_indices (ix,iy)   function calc_Λ_m!(Λ_m, Λ, t)
                                                if (ix <= size(Λ_m, 1) && iy <= size(Λ_m, 2))
                                                    Λ_m[ix, iy] = Λ * calc_m(ix, iy, t)
                                                end
                                                return
                                            end
    output = runthemodel_scaled(scaled_params, ϕ_init, h_init, calc_Λ_m!, ice_mask, bc_diric, bc_no_xflux, bc_no_yflux, update_h_only, do_print, warmup)
    output_descaled = descale(output, N_, ϕ_, h_, q_)
    return output_descaled
end
