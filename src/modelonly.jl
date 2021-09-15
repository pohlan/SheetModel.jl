using LazyArrays: Diff
using Printf, Infiltrator

# used indices:
# - normal grid (i,j), size (nx, ny)
#   ϕ, h, Res_ϕ, Res_h, ..
# - qx grid (m,n) staggered in x-dir.: (nx-1, ny)
# - qy grid (m,n) staggered in y-dir.: (nx, ny-1)

"""
Calculate discharge; qx = calc_q(h, dϕ_dx, dϕ_dy, ...), qy = calc_q(h, dϕ_dy, dϕ_dx, ...)
"""
calc_q(h, dϕ_du1, dϕ_du2, k, α, β, small) = - k * h^α * (sqrt(dϕ_du1^2 + dϕ_du2^2) + small)^(β-2) * dϕ_du1

"""
Calculate water pressure
"""
calc_pw(ϕ, ρw, g, zb) = ϕ - ρw * g * zb

"""
Calculate effective pressure
"""
calc_N(ϕ, ρi, ρw, g, H, zb) = ρi * g * H - calc_pw(ϕ, ρw, g, zb)


"""
Calculate closure rate
"""
function calc_vc(ϕ, h, ρi, ρw, g, H, zb, n, A)
    N = calc_N(ϕ, ρi, ρw, g, H, zb)
    return 2 / n^n * A * h * abs(N)^(n-1) * N
end

"""
Calculate opening rate
"""
calc_vo(h, ub, hr, lr) = h < hr ? ub * (hr - h) / lr : 0.0


macro pw(ix, iy) esc(:(ϕ[$ix, $iy] - ρw * g * zb[$ix, $iy])) end

macro N(ix, iy) esc(:(ρi * g * H[$ix, $iy] - @pw($ix, $iy))) end

macro vc(ix, iy) esc(:(2 / n^n * A * h[$ix, $iy] * abs(@N($ix, $iy))^(n-1) * @N($ix, $iy))) end # scaled version

macro vo(ix, iy) esc(:(h[$ix, $iy] < hr ? ub * (hr - h[$ix, $iy]) / lr : 0.0)) end # scaled version

macro dϕ_dx(ix, iy) esc(:( (qx_ice[$ix, $iy] == 2) #* !qx_xubound[$ix, $iy] * !qx_xlbound[$ix, $iy] # leave qx_ice away?
                           * (ϕ[$ix+1, $iy] - ϕ[$ix, $iy]) / dx
                        )) end
macro dϕ_dy(ix, iy) esc(:( (qy_ice[$ix, $iy] == 2) #* !qy_yubound[$ix, $iy] * !qy_ylbound[$ix, $iy] # leave qy_ice away?
                            * (ϕ[$ix, $iy+1] - ϕ[$ix, $iy]) / dy
                         )) end

macro gradϕ(ix, iy) esc(:( sqrt(                                                      # gradϕ only defined on interior points
                                  (0.5 * (@dϕ_dx($ix+1, $iy+1) + @dϕ_dx($ix, $iy+1)))^2
                                + (0.5 * (@dϕ_dy($ix+1, $iy+1) + @dϕ_dy($ix+1, $iy)))^2
                                  ))) end
macro d_eff(ix, iy) esc(:( k * h[$ix+1, $iy+1]^α * (@gradϕ($ix, $iy) + small)^(β-2) )) end # d_eff only defined on interior points

#macro dτ_ϕ(ix, iy) esc(:( dτ_ϕ_ *  (1.0 ./ (min(dx, dy)^2 ./ @d_eff($ix, $iy) / 4.1) .+ 1.0 / dt) .^(-1))) end # other definitions...
macro dτ_ϕ(ix, iy) esc(:( dτ_ϕ_ * min(min(dx, dy)^2 / @d_eff($ix, $iy) / 4.1, dt))) end


### KERNEL functions ###

function output_params!(N, ϕ, ρi, ρw, g, H, zb, qx, qy, dx, dy, k, h, α, β, small, qx_ice, qy_ice, qx_xlbound, qx_xubound, qy_ylbound, qy_yubound)
    for iy = 1:size(ϕ, 2)
        for ix = 1:size(ϕ, 1)
            N[ix, iy] = @N(ix, iy)
        end
    end
    return
end

function update_difference!(Δϕ, ϕ, ϕ2, Δh, h, h2)
    #Threads.@threads for iy=1:size(ϕ, 2)
    for iy = 1:size(ϕ, 2)
        for ix = 1:size(ϕ, 1)
            Δϕ[ix, iy] = abs(ϕ[ix, iy] - ϕ2[ix, iy])
            Δh[ix, iy] = abs(h[ix, iy] - h2[ix, iy])
        end
    end
    return
end

"""
Calculate fluxes in x-direction using upstream scheme
"""
function flux_x!(qx, qy, ϕ, dx, dy, k, h, α, β, small, qx_ice, qy_ice, qx_xlbound, qx_xubound, qy_ylbound, qy_yubound)
    nx, ny = size(ϕ)
    #Threads.@threads for iy=1:ny
    for iy = 2:ny-1
        for ix = 2:nx-2
            qx[ix, iy] = - @d_eff(ix, iy-1)   * @dϕ_dx(ix, iy) * (@dϕ_dx(ix, iy) >= 0) +   # flux in negative x-direction
                         - @d_eff(ix-1, iy-1) * @dϕ_dx(ix, iy) * (@dϕ_dx(ix, iy) <  0)     # flux in positive x-direction
        end
    end
    return
end

"""
Calculate fluxes in y-direction using upstream scheme
"""
function flux_y!(qx, qy, ϕ, dx, dy, k, h, α, β, small, qx_ice, qy_ice, qx_xlbound, qx_xubound, qy_ylbound, qy_yubound)
    nx, ny = size(ϕ)
    #Threads.@threads for iy=1:ny-1
    for iy = 2:ny-2
        for ix = 2:nx-1
            qy[ix, iy] = - @d_eff(ix-1, iy)   * @dϕ_dy(ix, iy) * (@dϕ_dy(ix, iy) >= 0) +   # flux in negative y-direction
                         - @d_eff(ix-1, iy-1) * @dϕ_dy(ix, iy) * (@dϕ_dy(ix, iy) <  0)     # flux in positive y-direction
        end
    end
    return
end

"""
Calculate residuals of ϕ & h and update the fields
"""
function update_fields!(Res_ϕ, Res_h, dϕ_dτ, dh_dτ, idx_ice, qx_ice, qy_ice, qx_xlbound, qx_xubound, qy_ylbound, qy_yubound, dx, dy,
                            k, α, β, small, qx, qy,
                            ϕ, ϕ2, ϕ_old, h, h2, h_old, dt, ev, m, hr, lr, ub, g, ρw, ρi, A, n, H, zb, Σ, Γ, Λ,
                            γ_ϕ, γ_h, dτ_h_, dτ_ϕ_)
    nx, ny = size(ϕ)
    #Threads.@threads for iy=1:ny
    for iy = 1:ny
        for ix = 1:nx
            if (1 < ix < nx) && (1 < iy < ny)
                # residual of ϕ
                if ix == 2 # ϕ: without boundary points (divergence of q not defined there)
                    Res_ϕ[ix, iy] = 0. # position where dirichlet b.c. are imposed
                else
                    Res_ϕ[ix, iy] = idx_ice[ix, iy] * (
                                    - ev/(ρw*g) * (ϕ[ix, iy] - ϕ_old[ix, iy]) / dt                               # dhe/dt
                                    - ( (qx[ix, iy] - qx[ix-1, iy]) / dx + (qy[ix, iy] - qy[ix, iy-1]) / dy )    # divergence
                                    - (Σ * @vo(ix, iy) - Γ * @vc(ix, iy))                                        # dh/dt
                                    + Λ * m[ix, iy]                                                              # source term
                                    )
                end

                # update ϕ
                dϕ_dτ[ix, iy] = Res_ϕ[ix, iy] + γ_ϕ * dϕ_dτ[ix, iy]
                ϕ2[ix, iy] = ϕ[ix, iy] + @dτ_ϕ(ix-1, iy-1) * dϕ_dτ[ix, iy]
            end

            # residual of h
            Res_h[ix, iy] = idx_ice[ix, iy] * (
                            - (h[ix, iy] - h_old[ix, iy]) / dt
                            + (Σ * @vo(ix, iy) - Γ * @vc(ix, iy))
                            )

            # update h
            dh_dτ[ix, iy] = Res_h[ix, iy] + γ_h * dh_dτ[ix, iy]
            h2[ix, iy] = h[ix, iy] + dτ_h_ * dh_dτ[ix, iy]
        end
    end
    return
end

"""
Apply Dirichlet boundary conditions to ϕ, at the moment pw(x=0) = 0
"""
function apply_bc!(ϕ, h, H, ρw, g, zb) # TODO: don't hard-wire, give bc as input parameters
    nx, ny = size(ϕ)
    #Threads.@threads for iy=1:ny
    for iy = 1:ny
        for ix = 1:nx
            if H[ix, iy] == 0.             # zero sheet thickness outside of glacier domain; necessary ??
                h[ix, iy] = 0.
            end
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
    end
    return
end


"""
Run the model with scaled parameters
"""
function runthemodel_scaled(params::Para, ϕ0, h0, printit, printtime)
    @unpack ev, g, ρw, ρi, n, A, Σ, Γ, Λ, calc_m_t, dx, dy, nx, ny, k, α, β, small,
            H, zb, ub, hr, lr, dt, ttot, tol, itMax, γ_ϕ, γ_h, dτ_ϕ_, dτ_h_ = params

    # Array allocation
    Δϕ, Δh, qx, qy,  m, N,
    dϕ_dτ, dh_dτ, Res_ϕ, Res_h = array_allocation(params)

    # determine indices of glacier domain
    idx_ice = H .> 0.0
    qx_ice = zeros(Int, nx-1, ny)
    qx_ice  .= idx_ice[1:end-1, :]
    qx_ice .+= idx_ice[2:end, :]   # qx_ice = 1 for boundary, 2 for interior
    qy_ice = zeros(Int, nx, ny-1)
    qy_ice  .= idx_ice[:, 1:end-1]
    qy_ice .+= idx_ice[:, 2:end]
    qx_xlbound  = Diff(idx_ice, dims=1) .== 1 # on qx grid
    qx_xubound  = Diff(idx_ice, dims=1) .== -1
    qy_ylbound  = Diff(idx_ice, dims=2) .== 1 # on qy grid
    qy_yubound  = Diff(idx_ice, dims=2) .== -1

    # Apply boundary conditions
    apply_bc!(ϕ0, h0, H, ρw, g, zb)

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
    t = 0.0; tstep=0; ittot = 0; t_tic = Base.time()

    # Physical time loop
    while t<ttot
        iter = 0
        err_ϕ_tol, err_h_tol = 2*tol, 2*tol

        m .= calc_m_t::Float64(t+dt)
        # Pseudo-transient iteration
        while !(max(err_ϕ_tol, err_h_tol) < tol) && iter<itMax # with the ! the loop also continues for NaN values of err

            flux_x!(qx, qy, ϕ, dx, dy, k, h, α, β, small, qx_ice, qy_ice, qx_xlbound, qx_xubound, qy_ylbound, qy_yubound)
            flux_y!(qx, qy, ϕ, dx, dy, k, h, α, β, small, qx_ice, qy_ice, qx_xlbound, qx_xubound, qy_ylbound, qy_yubound)
            update_fields!(Res_ϕ, Res_h, dϕ_dτ, dh_dτ, idx_ice, qx_ice, qy_ice, qx_xlbound, qx_xubound, qy_ylbound, qy_yubound, dx, dy,
                            k, α, β, small, qx, qy,
                            ϕ, ϕ2, ϕ_old, h, h2, h_old, dt, ev, m, hr, lr, ub, g, ρw, ρi, A, n, H, zb, Σ, Γ, Λ,
                            γ_ϕ, γ_h, dτ_h_, dτ_ϕ_)

            # apply boundary conditions
            apply_bc!(ϕ2, h2, H, ρw, g, zb)

            # switch pointer
            ϕ, ϕ2 = ϕ2, ϕ
            h, h2 = h2, h

            # determine the errors (only consider points where the ice thickness is > 0)

            if iter % 100 == 0
                # residual error
                err_ϕ_res = norm(Res_ϕ[idx_ice]) / sum(idx_ice) # or length(Res_ϕ) instead of sum(idx_ice) ??
                err_h_res = norm(Res_h[idx_ice]) / norm(h0)
                if (iter==0)
                    err_ϕ_ini = err_ϕ_res
                    err_h_ini = err_h_res
                end
                err_ϕ_resrel = err_ϕ_res / err_ϕ_ini
                err_h_resrel = err_h_res / err_h_ini

                # update error
                update_difference!(Δϕ, ϕ, ϕ2, Δh, h, h2)
                err_ϕ = norm(Δϕ[idx_ice]) / sum(idx_ice)
                err_h = norm(Δh[idx_ice]) / norm(h0)
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

            end

            if iter % printit == 0
                @printf("iterations = %d, error ϕ = %1.2e, error h = %1.2e \n", iter, err_ϕ_tol, err_h_tol)
            end

            iter += 1

        end
        ittot += iter-1; tstep += 1; t += dt

        #if mod(tstep, printtime) == 0
        #    @printf("time step = %d, number of iterations = %d \n", tstep, iter)
        #end

        ϕ_old .= ϕ
        h_old .= h
    end

    # Perfomance measures
    t_toc = Base.time() - t_tic                # execution time, s
    A_eff = (2*5+2)/1e9*nx*ny*sizeof(Float64)  # effective main memory access per iteration [GB];
                                               # 5 read+write arrays (ϕ, dϕ_dτ, h, dh_dτ, m), 2 read arrays (ϕ_old, h_old) --> check!
    t_it  = t_toc/ittot                        # execution time per iteration, s
    T_eff = A_eff/t_it                         # effective memory throughput, GB/s
    @printf("Time = %1.3f sec, T_eff = %1.2f GB/s (iterations total = %d)\n", t_toc, round(T_eff, sigdigits=2), ittot)

    # calculate N, qx and qy as output parameters
    output_params!(N, ϕ, ρi, ρw, g, H, zb, qx, qy, dx, dy, k, h, α, β, small, qx_ice, qy_ice, qx_xlbound, qx_xubound, qy_ylbound, qy_yubound)

    return model_output(N=N, ϕ=ϕ, h=h, qx=qx, qy=qy, qx_ice=qx_ice, qy_ice=qy_ice,
            Err_ϕ=Δϕ, Err_h=Δh, Res_ϕ=Res_ϕ, Res_h=Res_h,
            ittot=ittot, iters=iters,
            errs_ϕ=errs_ϕ, errs_h=errs_h,
            errs_ϕ_rel=errs_ϕ_rel, errs_h_rel=errs_h_rel,
            errs_ϕ_res=errs_ϕ_res, errs_h_res=errs_h_res,
            errs_ϕ_resrel=errs_ϕ_resrel, errs_h_resrel=errs_h_resrel)
end

"""
Scale the parameters and call the model run function
"""
function runthemodel(input::Para, ϕ0, h0;
                    printit=10^5,         # error is printed after `printit` iterations of pseudo-transient time
                    printtime=10^5)       # time step and number of PT iterations is printed after `printtime` number of physical time steps
    params, ϕ0, h0, ϕ_, N_, h_, q_ = scaling(input, ϕ0, h0)
    output = runthemodel_scaled(params::Para, ϕ0, h0, printit, printtime)
    output_descaled = descaling(output, N_, ϕ_, h_, q_)
    return output_descaled
end
