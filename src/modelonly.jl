using LazyArrays: Diff
using Printf, Infiltrator

# used grids:
# - normal grid (i,j), size (nx, ny)
#   ϕ, h, Res_ϕ, Res_h, ..
# - d_eff grid (p, q), size (nx-2, ny-2)
#   only inner points of ϕ/h grid, e.g. p=1 corresponds to i=2
# - qx grid (m,n) staggered in x-dir.: (nx-1, ny)
# - qy grid (m,n) staggered in y-dir.: (nx, ny-1)


### MACROS ###

"""
Calculate water pressure; input coordinates on ϕ/h grid.
Needs access to ϕ, ρw, g, zb.
"""
macro pw(ix, iy) esc(:(ϕ[$ix, $iy] - ρw * g * zb[$ix, $iy])) end

"""
Calculate effective pressure; input coordinates on ϕ/h grid.
Needs access to ϕ, H, ρw, ρi, g, zb
"""
macro N(ix, iy) esc(:(ρi * g * H[$ix, $iy] - @pw($ix, $iy))) end

"""
Calculate closure rate; input coordinates on ϕ/h grid.
Needs access to ϕ, H, ρw, ρi, g, zb, n, A, h
"""
macro vc(ix, iy) esc(:(2 / n^n * A * h[$ix, $iy] * abs(@N($ix, $iy))^(n-1) * @N($ix, $iy))) end

"""
Calculate opening rate; input coordinates on ϕ/h grid.
Needs access to h, hr, ub, lr
"""
macro vo(ix, iy) esc(:(h[$ix, $iy] < hr ? ub * (hr - h[$ix, $iy]) / lr : 0.0)) end

"""
Calculate hydraulic gradient in x-direction; input coordinates on qx grid.
Implicitly sets zero-flux boundary conditions.
Needs access to H, ϕ, dx
"""
macro dϕ_dx(ix, iy) esc(:( (H[$ix, $iy] > 0. && H[$ix+1, $iy] > 0.) # only consider ice interior; gradients at boundary and outside are zero
                           * (ϕ[$ix+1, $iy] - ϕ[$ix, $iy]) * dx_
                        )) end

"""
Calculate hydraulic gradient in y-direction; input coordinates on qy grid.
Implicitly sets zero-flux boundary conditions.
Needs access to H, ϕ, dy
"""
macro dϕ_dy(ix, iy) esc(:( (H[$ix, $iy] > 0. && H[$ix, $iy+1] > 0.) # only consider ice interior; gradients at boundary and outside are zero
                            * (ϕ[$ix, $iy+1] - ϕ[$ix, $iy]) * dy_
                         )) end

"""
Calculate absolute hydraulic gradient, |∇ϕ|;
input coordinates on ϕ/h grid (nx, ny), but only possible to calculate on inner points.
Needs access to ϕ, dx, dy, H
"""
macro gradϕ(ix, iy) esc(:( sqrt(
                                  (0.5 * (@dϕ_dx($ix, $iy) + @dϕ_dx($ix-1, $iy)))^2
                                + (0.5 * (@dϕ_dy($ix, $iy) + @dϕ_dy($ix, $iy-1)))^2
                                  ))) end
"""
Calculate effective diffusivity;
input coordinates on ϕ/h grid (nx, ny), but only possible to calculate on inner points.
Needs access to k, h, α, β, small, ϕ, dx, dy, H
"""
macro d_eff(ix, iy) esc(:( k * h[$ix, $iy]^α * (@gradϕ($ix, $iy) + small)^(β-2) )) end

"""
Calculate pseudo-time step of ϕ,
input coordinates on ϕ/h grid (nx, ny), but only possible to calculate on inner points.
Needs access to dτ_ϕ_, dx, dy, dt, k, h, α, β, small, ϕ, H
"""
macro dτ_ϕ(ix, iy) esc(:( dτ_ϕ_ * min(min_dxy2 / d_eff[$ix, $iy] / 4.1, dt))) end
#macro dτ_ϕ(ix, iy) esc(:( dτ_ϕ_ *  (min_dxy2_ / d_eff[$ix, $iy] / 4.1) .+ 1.0 * dt_) .^(-1))) end # other definitions...

"""
Calculate flux in x-direction using an upstream scheme; input coordinates on qx grid.
Needs access to k, h, α, β, small, ϕ, dx, dy, H
"""
macro flux_x(ix, iy) esc(:(
    - d_eff[$ix+1, $iy] * max(@dϕ_dx($ix, $iy), 0) +   # flux in negative x-direction
    - d_eff[$ix,   $iy] * min(@dϕ_dx($ix, $iy), 0.)    # flux in positive x-direction
    )) end
"""
Calculate flux in y-direction using an upstream scheme; input coordinates on qy grid.
Needs access to k, h, α, β, small, ϕ, dx, dy, H
"""
macro flux_y(ix, iy) esc(:(
    - d_eff[$ix, $iy+1] * max(@dϕ_dy($ix, $iy), 0.) +   # flux in negative y-direction
    - d_eff[$ix, $iy  ] * min(@dϕ_dy($ix, $iy), 0.)    # flux in positive y-direction
    )) end


"""
Calculate residual of ϕ; input coordinates on ϕ/h grid (but only defined on inner grid points 2:nx-1, 2:ny-1, due to d_eff).
Needs access to ϕ, ϕ_old, h, H, ev, ρw, ρi, g, k, α, β, small, dx, dy, dt, Σ, Γ, Λ, m, hr, ub, lr, zb, n, A
"""
macro Res_ϕ(ix, iy) esc(:(( H[$ix, $iy] > 0.) * (                                                      # only calculate at points with non-zero ice thickness
                                                  - ev/(ρw*g) * (ϕ[$ix, $iy] - ϕ_old[$ix, $iy]) * dt_                                                   # dhe/dt
                                                  - ( (@flux_x($ix, $iy) - @flux_x($ix-1, $iy)) * dx_ + (@flux_y($ix, $iy) - @flux_y($ix, $iy-1)) * dy_ )    # divergence
                                                  - (Σ * @vo($ix, $iy) - Γ * @vc($ix, $iy))                                                            # dh/dt
                                                  + Λ * m[$ix, $iy]                                                                                  # source term
                                                 )
)) end

"""
Calculate residual of h; input coordinates on ϕ/h grid.
Needs access to ϕ, h, h_old, H, dt, Σ, Γ, hr, ub, lr, zb, n, A, ρw, ρi, g
"""
macro Res_h(ix, iy) esc(:(( H[$ix, $iy] > 0.) * (
                                                  - (h[$ix, $iy] - h_old[$ix, $iy]) * dt_
                                                  + (Σ * @vo($ix, $iy) - Γ * @vc($ix, $iy))
                                                 )
)) end

### KERNEL FUNCTIONS ###

@parallel_indices (ix, iy) function update_deff!(d_eff, ϕ, h, dx_, dy_, k, α, β, dt, ev, hr, lr, ub, g, ρw, ρi, A, n, H, zb, small)
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
@parallel_indices (ix,iy) function residuals!(ϕ, ϕ_old, h, h_old, Res_ϕ, Res_h, qx, qy, m, d_eff,
                                              dx_, dy_, min_dxy2, k, α, β, dt, dt_, ev, hr, lr, ub, g, ρw, ρi, A, n, H, zb, Σ, Γ, Λ, small)
    nx, ny = size(ϕ)
    if (ix <= nx && iy <= ny)
        # residual of ϕ
        if ix == 2 # ϕ: without boundary points (divergence of q not defined there)
            Res_ϕ[ix, iy] = 0. # position where dirichlet b.c. are imposed
        elseif (1 < ix < nx && 1 < iy < ny)
            Res_ϕ[ix, iy] = @Res_ϕ(ix, iy)
        end

        # residual of h
        Res_h[ix, iy] = @Res_h(ix, iy)
    end
    return
end

"""
Update the fields of ϕ and h using the pseudo-transient method with damping.
"""
@parallel_indices (ix,iy) function update_fields!(ϕ, ϕ2, ϕ_old, h, h2, h_old, qx, qy, m, d_eff, iter,
                                                  dx_, dy_, min_dxy2, k, α, β, dt, dt_, ev, hr, lr, ub, g, ρw, ρi, A, n, H, zb, Σ, Γ, Λ, small,
                                                  dϕ_dτ, dh_dτ, γ_ϕ, γ_h, dτ_h_, dτ_ϕ_)
    nx, ny = size(ϕ)
    if (1 < ix < nx && 1 < iy < ny)
        # update ϕ
        dϕ_dτ[ix, iy] = @Res_ϕ(ix, iy) + γ_ϕ * dϕ_dτ[ix, iy]
        ϕ2[ix, iy] = ϕ[ix, iy] + @dτ_ϕ(ix, iy) * dϕ_dτ[ix, iy]
        # dirichlet boundary conditions to pw = 0
        if ix == 2
            ϕ2[ix, iy] = ρw * g * zb[ix, iy]
        end

        # update h
        dh_dτ[ix, iy] = @Res_h(ix, iy) + γ_h * dh_dτ[ix, iy]
        h2[ix, iy] = h[ix, iy] + dτ_h_ * dh_dτ[ix, iy]
    end
    return
end

"""
Apply Dirichlet boundary conditions to ϕ, at the moment pw(x=0) = 0.
Neumann boundary conditions are applied when calculating the hydraulic gradient (@dϕ_dx and @dϕ_dy macros).
"""
@parallel_indices (ix,iy) function apply_bc!(ϕ, h, H, ρw, g, zb) # TODO: don't hard-wire, give bc as input parameters
    nx, ny = size(ϕ)
    if (ix <= nx && iy <= ny)
        #if H[ix, iy] == 0.             # zero sheet thickness outside of glacier domain; necessary ??
        #    h[ix, iy] = 0.
        #end
        if ix == 2
            ϕ[ix, iy] = ρw * g * zb[ix, iy]
        end
        #if (ix == 2) && (iy == ny÷2+1)
        #    ϕ[2, ny÷2+1] = ρw .* g .* zb[2, ny÷2+1]
        #end
        #if (ix == 2) && (iy == ny÷2) && iseven(ny)
        #    ϕ[2, ny÷2] = ρw .* g .* zb[2, ny÷2]
        #end
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
@parallel_indices (ix,iy) function output_params!(N, qx, qy, ϕ, h, d_eff,
                                                  ρi, ρw, g, H, zb, k, α, β, dx_, dy_, small)
    nx, ny = size(ϕ)
    if (ix <= nx && iy <= ny)
        N[ix, iy]  = @N(ix, iy)
        if (1 < ix < nx && 1 < iy < ny)
            qx[ix, iy] = @flux_x(ix, iy)
            qy[ix, iy] = @flux_y(ix, iy)
        end
    end
    return
end

"""
Run the model with scaled parameters.
"""
@views function runthemodel_scaled(params::Para, ϕ0, h0, printtime)
    @unpack ev, g, ρw, ρi, n, A, Σ, Γ, Λ, calc_m_t, dx, dy, nx, ny, k, α, β,
            H, zb, ub, hr, lr, dt, ttot, tol, itMax, γ_ϕ, γ_h, dτ_ϕ_, dτ_h_ = params

    # Pre-calculate reciprocals for better performance
    min_dxy2 = min(dx, dy)^2
    dx_, dy_, dt_, min_dxy2_ = (dx, dy, dt, min_dxy2).^(-1)

    # Array allocation
    Δϕ, Δh, qx, qy, d_eff, m, N,
    dϕ_dτ, dh_dτ, Res_ϕ, Res_h = array_allocation(params)

    # Apply boundary conditions
    @parallel apply_bc!(ϕ0, h0, H, ρw, g, zb)

    ϕ_old   = copy(ϕ0)
    ϕ       = copy(ϕ0)
    ϕ2      = copy(ϕ0)
    h_old   = copy(h0)
    h       = copy(h0)
    h2      = copy(h0)

    # for iterations vs. error plot
    iters      = Int64[]
    errs_ϕ     = Float64[]
    errs_h     = Float64[]
    errs_ϕ_rel = Float64[]
    errs_h_rel = Float64[]
    errs_ϕ_res = Float64[]
    errs_h_res = Float64[]
    errs_ϕ_resrel = Float64[]
    errs_h_resrel = Float64[]
    err_ϕ = 0.
    err_h = 0.
    err_ϕ_rel = 0.
    err_h_rel = 0.
    err_ϕ_res = 0.
    err_h_res = 0.
    err_ϕ_resrel = 0.
    err_h_resrel = 0.
    err_ϕ_ini = 0.0
    err_h_ini = 0.0

    # initiate time loop parameters
    t = 0.0; tstep=0; ittot = 0; t_tic = 0.
    steady_state = false;

    # Physical time loop
    while t<ttot
        iter = 0
        err_ϕ_tol, err_h_tol = 2*tol, 2*tol

        m .= Data.Array(calc_m_t(t+dt))
        # Pseudo-transient iteration
        while !(max(err_ϕ_tol, err_h_tol) < tol) && iter<itMax # with the ! the loop also continues for NaN values of err

            # don't consider first ten iterations for performance measure
            if (iter == 10) t_tic = Base.time() end

            # update d_eff (the bottleneck in performance, only do it every 10 iterations)
            # but it has to be done in a seperate kernel to ensure every grid point is accessing the same version of ϕ
            #if iter % 1000 == 0.
                @parallel update_deff!(d_eff, ϕ, h, dx_, dy_, k, α, β, dt, ev, hr, lr, ub, g, ρw, ρi, A, n, H, zb, small)
            #end
            # update ϕ and h
            @parallel update_fields!(ϕ, ϕ2, ϕ_old, h, h2, h_old, qx, qy, m, d_eff, iter,
                                     dx_, dy_, min_dxy2, k, α, β, dt, dt_, ev, hr, lr, ub, g, ρw, ρi, A, n, H, zb, Σ, Γ, Λ, small,
                                     dϕ_dτ, dh_dτ, γ_ϕ, γ_h, dτ_h_, dτ_ϕ_)

            # pointer swap
            ϕ, ϕ2 = ϕ2, ϕ
            h, h2 = h2, h

            iter += 1

            # determine the errors (only consider points where the ice thickness is > 0)

            if iter % 1000 == 0 && itMax > 10^4     # only calculate errors every 1000 time steps
                                                    # and if itMax is high, i.e. if convergence is the goal
                # update the residual arrays
                @parallel residuals!(ϕ, ϕ_old, h, h_old, Res_ϕ, Res_h, qx, qy, m, d_eff,
                                     dx_, dy_, min_dxy2, k, α, β, dt, dt_, ev, hr, lr, ub, g, ρw, ρi, A, n, H, zb, Σ, Γ, Λ, small)

                # residual error
                err_ϕ_res = norm(Res_ϕ) / length(Res_ϕ) # or length(Res_ϕ) instead of sum(H .> 0.) ??
                err_h_res = norm(Res_h) / norm(h0)
                if (iter==0)
                    err_ϕ_ini = err_ϕ_res
                    err_h_ini = err_h_res
                end
                err_ϕ_resrel = err_ϕ_res / err_ϕ_ini
                err_h_resrel = err_h_res / err_h_ini

                # update error
                @parallel update_difference!(Δϕ, ϕ, ϕ2, Δh, h, h2)
                err_ϕ = norm(Δϕ) / length(Δϕ)
                err_h = norm(Δh) / norm(h0)
                if (iter==0)
                    err_ϕ_ini = err_ϕ
                    err_h_ini = err_h
                end
                err_ϕ_rel = err_ϕ / err_ϕ_ini
                err_h_rel = err_h / err_h_ini

                # decide which errors should be below the tolerance and be printed out
                err_ϕ_tol, err_h_tol = err_ϕ, err_h

                # save error evolution in vector
                append!(iters, iter)
                append!(errs_ϕ, err_ϕ)
                append!(errs_h, err_h)
                append!(errs_ϕ_rel, err_ϕ_rel)
                append!(errs_h_rel, err_h_rel)
                append!(errs_ϕ_res, err_ϕ_res)
                append!(errs_h_res, err_h_res)
                append!(errs_ϕ_resrel, err_ϕ_resrel)
                append!(errs_h_resrel, err_h_resrel)

                @printf("iterations = %d, error ϕ = %1.2e, error h = %1.2e \n", iter, err_ϕ_tol, err_h_tol)

            end


        end
        ittot += iter; tstep += 1; t += dt

        #if mod(tstep, printtime) == 0
        #    @printf("time step = %d, number of iterations = %d \n", tstep, iter)
        #end

        @parallel old2new!(ϕ, ϕ_old, h, h_old)
        if max(err_ϕ_tol, err_h_tol) < tol && t > 1e7
            steady_state = true
        end
    end

    # Perfomance measures
    t_toc = Base.time() - t_tic                # execution time, s
    A_eff = (5+10)/1e9*nx*ny*sizeof(Float64)  # effective main memory access per iteration [GB];
                                               # 5 write arrays (dϕ_dτ, ϕ2, dh_dτ, h2, d_eff)
                                               # 2 read arrays (ϕ, ϕ_old, dϕ_dτ, h, h_old, dh_dτ, H, zb, m, d_eff)
    t_it  = t_toc/(ittot-10)                   # execution time per iteration, s
    T_eff = A_eff/t_it                         # effective memory throughput, GB/s
    @printf("Time = %1.3f sec, T_eff = %1.1f GB/s, iterations total = %d, (nx, ny) = (%d, %d)\n", t_toc, round(T_eff, sigdigits=3), ittot, nx, ny)

    # calculate N, qx and qy as output parameters
    @parallel output_params!(N, qx, qy, ϕ, h, d_eff, ρi, ρw, g, H, zb, k, α, β, dx_, dy_, small)

    return model_output{Data.Array}(;   N, ϕ, h, qx, qy,
                            Err_ϕ=Δϕ, Err_h=Δh, Res_ϕ, Res_h,
                            ittot, iters,
                            errs_ϕ, errs_h, errs_ϕ_rel, errs_h_rel,
                            errs_ϕ_res, errs_h_res, errs_ϕ_resrel, errs_h_resrel,
                            time_tot=t_toc, T_eff, steady_state)
end

"""
Scale the parameters and call the model run function.
"""
@views function runthemodel(input::Para, ϕ0, h0;
                    printtime=10^5)       # time step is printed after `printtime` number of physical time steps
    params, ϕ0, h0, ϕ_, N_, h_, q_ = scaling(input, ϕ0, h0)
    output = runthemodel_scaled(params::Para, ϕ0, h0, printtime)
    output_descaled = descaling(output, N_, ϕ_, h_, q_)
    return output_descaled
end
