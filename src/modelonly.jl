using LazyArrays: Diff
using Printf, Infiltrator

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

@views av(A)    = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end]) # average
@views av_xa(A) = 0.5.*(A[1:end-1,:].+A[2:end,:]) # average x-dir
@views av_ya(A) = 0.5.*(A[:,1:end-1].+A[:,2:end]) # average y-dir

"""
Apply Dirichlet boundary conditions to ϕ, at the moment pw(x=0) = 0
"""
function apply_bc(ϕ, h, H, ρw, g, zb) # TODO: shouldn't have any function with arrays as arguments
                                   # TODO: don't hard-wire, give bc as input parameters
    nx, ny = size(ϕ)
    #ϕ[H .== 0.0] .= ρw .* g .* zb[H .== 0.0] # zero water pressure outside of glacier domain
    h[H .== 0.0] .= 0.0                       # zero sheet thickness outside of glacier domain
    ϕ[2, :] .= ρw .* g .* zb[2, :]

    # ϕ[end-1,:] = ρw .* g .* zb[end-1,:]
    # ϕ[2, ny÷2+1] = ρw .* g .* zb[2, ny÷2+1]
    # if iseven(size(ϕ,2))
    #     ϕ[2, ny÷2] = ρw .* g .* zb[2, ny÷2]
    # end

    return ϕ, h
end

"""
Run the model with scaled parameters
"""
function runthemodel_scaled(params::Para, ϕ0, h0, printit, printtime)
    @unpack ev, g, ρw, ρi, n, A, Σ, Γ, Λ, calc_m_t, dx, dy, nx, ny, k, α, β, small,
            H, zb, ub, hr, lr, dt, ttot, tol, itMax, γ_ϕ, γ_h, dτ_ϕ_, dτ_h_ = params
    # Array allocation
    vo, vc, dϕ_dx, dϕ_dy, qx, qy, gp_ice, ix, iy, m, div_q, div_ϕ,
    dϕ_dτ, dh_dτ, Res_ϕ, Res_h, Err_ϕ, Err_h, d_eff, dτ_ϕ = array_allocation(params)

    # determine indices of glacier domain
    idx_ice = H .> 0.0
    qx_ice = zeros(Int, nx-1, ny)
    qx_ice  .= idx_ice[1:end-1, :]
    qx_ice .+= idx_ice[2:end, :]   # qx_ice = 1 for boundary, 2 for interior
    qy_ice = zeros(Int, nx, ny-1)
    qy_ice  .= idx_ice[:, 1:end-1]
    qy_ice .+= idx_ice[:, 2:end]
    qx_xlbound  = Array(Diff(idx_ice, dims=1) .== 1) # on qx grid
    qx_xubound  = Array(Diff(idx_ice, dims=1) .== -1)
    qy_ylbound  = Array(Diff(idx_ice, dims=2) .== 1) # on qy grid
    qy_yubound  = Array(Diff(idx_ice, dims=2) .== -1)

    # Apply boundary conditions
    ϕ0, h0 = apply_bc(ϕ0, h0, H, ρw, g, zb)

    ϕ_old   = copy(ϕ0)
    ϕ       = copy(ϕ0)
    h_old   = copy(h0)
    h       = copy(h0)

    # initiate time loop parameters
    t, it, ittot = 0.0, 0, 0

    # for iterations vs. error plot
    errs_ϕ     = Float64[]
    errs_h     = Float64[]
    errs_ϕ_rel = Float64[]
    errs_h_rel = Float64[]
    errs_ϕ_res = Float64[]
    errs_h_res = Float64[]
    errs_ϕ_resrel = Float64[]
    errs_h_resrel = Float64[]
    err_ϕ_ini = 0.0
    err_h_ini = 0.0

    # Physical time loop
    while t<ttot
        iter = 0
        err_ϕ_tol, err_h_tol = 2*tol, 2*tol

        m .= calc_m_t(t+dt)
        # Pseudo-transient iteration
        while !(max(err_ϕ_tol, err_h_tol) < tol) && iter<itMax # with the ! the loop also continues for NaN values of err
        #while !(err_ϕ < tol) && iter<itMax # only solving for ϕ
            # used indices:
            # - normal grid (i,j), size (nx, ny)
            #   ϕ, h, Res_ϕ, Res_h, ..
            # - qx grid (m,n) staggered in x-dir.: (nx-1, ny)
            #   qx, dϕ_dx[m,j]
            # - qy grid (m,n) staggered in y-dir.: (nx, ny-1)
            #   qy, dϕ_dy

            # save current ϕ and h for error calculation
            Err_ϕ   .= ϕ
            Err_h   .= h

            # hydraulic gradient
            dϕ_dx   .= Diff(ϕ, dims=1) ./ dx
            dϕ_dy   .= Diff(ϕ, dims=2) ./ dy

            # flux boundary conditions
            dϕ_dx[qx_ice .== 0] .= 0.0
            dϕ_dy[qy_ice .== 0] .= 0.0
            dϕ_dx[qx_xubound] .= 0.0
            dϕ_dx[qx_xlbound] .= 0.0
            dϕ_dy[qy_ylbound] .= 0.0
            dϕ_dy[qy_yubound] .= 0.0

            # effective diffusivity
            gradϕ = sqrt.(av_xa(dϕ_dx[:, 2:end-1]).^2 .+ av_ya(dϕ_dy[2:end-1, :]).^2) # on ϕ/h grid, size (nx-2, ny-2)
            d_eff .= k*h[2:end-1, 2:end-1].^α .* (gradϕ .+ small).^(β-2) # on ϕ/h grid, size (nx-2, ny-2)

            # calculate fluxes
            # upstream
            qx[2:end-1, 2:end-1] .= (dϕ_dx[2:end-1, 2:end-1] .>= 0.0) .* (.- d_eff[2:end, :]   .* dϕ_dx[2:end-1, 2:end-1]) .+
                                    (dϕ_dx[2:end-1, 2:end-1] .< 0.0)  .* (.- d_eff[1:end-1, :] .* dϕ_dx[2:end-1, 2:end-1])
            qy[2:end-1, 2:end-1] .= (dϕ_dy[2:end-1, 2:end-1] .>= 0.0) .* (.- d_eff[:, 2:end]   .* dϕ_dy[2:end-1, 2:end-1]) .+
                                    (dϕ_dy[2:end-1, 2:end-1] .< 0.0)  .* (.- d_eff[:, 1:end-1] .* dϕ_dy[2:end-1, 2:end-1])
            # central differences
            #qx[2:end-1, 2:end-1] .= 0.5 .* (.- d_eff[2:end, :]   .* dϕ_dx[2:end-1, 2:end-1]) .+
            #                        0.5 .* (.- d_eff[1:end-1, :] .* dϕ_dx[2:end-1, 2:end-1])
            #qy[2:end-1, 2:end-1] .= 0.5 .* (.- d_eff[:, 2:end]   .* dϕ_dy[2:end-1, 2:end-1]) .+
            #                        0.5 .* (.- d_eff[:, 1:end-1] .* dϕ_dy[2:end-1, 2:end-1])

            # set flux boundary conditions (not necessary here, it has to be imposed before calculating gradϕ and d_eff via dϕ_d...)
            #qx[qx_ice .== 0]   .= 0.0
            #qy[qy_ice .== 0]   .= 0.0
            #qx[qx_xubound] .= 0.0 # no flux boundary condition
            #qx[qx_xlbound] .= 0.0
            #qy[qy_ylbound] .= 0.0
            #qy[qy_yubound] .= 0.0

            vo     .= calc_vo.(h, ub, hr, lr)                 # opening rate
            vc     .= calc_vc.(ϕ, h, ρi, ρw, g, H, zb, n, A)  # closure rate
            div_q[2:end-1, 2:end-1]  .= Diff(qx, dims=1)[:, 2:end-1]/dx .+ Diff(qy, dims=2)[2:end-1, :]/dy .+ small

            # calculate residuals
            Res_ϕ   .=  idx_ice .* (
                            - ev/(ρw*g) * (ϕ .- ϕ_old)/dt .-         # dhe/dt
                            div_q .-                                 # div(q)
                            (Σ * vo .- Γ * vc)            .+         # dh/dt
                            Λ * m                                    # source term
                            )
            Res_h   .=  idx_ice .* (
                            - (h .- h_old) / dt  .+
                            (Σ * vo .- Γ * vc)
                            )

            # determine pseudo-time steps
            #dτ_ϕ[2:end-1, 2:end-1] .= dτ_ϕ_ .* (1.0 ./ (min(dx, dy)^2 ./ d_eff / 4.1) .+ 1.0 / dt) .^(-1)
            #dτ_ϕ[2:end-1, 2:end-1] .= dτ_ϕ_ .* (1.0 ./ (min(dx, dy)^2 ./ d_eff / 4.1)) .^(-1) # with this, if eq. are time-independent, #iterations is indep. of dt
            dτ_ϕ[2:end-1, 2:end-1] .= dτ_ϕ_ .* min.(min(dx, dy)^2 ./ d_eff / 4.1, dt)
            dτ_h   = dτ_h_   # pseudo-time step for h, scalar

            # damped rate of change
            dϕ_dτ .= Res_ϕ .+ γ_ϕ .* dϕ_dτ
            dh_dτ .= Res_h .+ γ_h .* dh_dτ

            # update fields
            ϕ .= ϕ .+ dτ_ϕ .* dϕ_dτ   # update ϕ
            h .= h .+ dτ_h .* dh_dτ   # update h

            # probably not going to use the following
            # ϕ, h = apply_bc(ϕ, h, H, ρw, g, zb)
            # vo     .= calc_vo.(h, ub, hr, lr)                 # opening rate
            # vc     .= calc_vc.(ϕ, h, ρi, ρw, g, H, zb, n, A)  # closure rate
            # h .= h_old .+ dt * (Σ * vo .- Γ * vc)

            # apply boundary conditions
            ϕ, h = apply_bc(ϕ, h, H, ρw, g, zb)

            # determine the errors (only consider points where the ice thickness is > 0)
            # error for ϕ

            # from Err - ϕ
            Err_ϕ .= abs.(Err_ϕ .- ϕ) # ./ dτ_ϕ
            #err_ϕ = norm(Err_ϕ[idx_ice]) ./ norm(Λ * m)
            err_ϕ = norm(Err_ϕ[idx_ice]) ./ sum(idx_ice)
            if (iter==1)  err_ϕ_ini = err_ϕ  end
            err_ϕ_rel = err_ϕ/err_ϕ_ini # relative error

            # from residual
            Res_ϕ[2, :] .= 0.0 # residual at b.c. should not be part of ϕ error
            #err_ϕ_res = norm(Res_ϕ[idx_ice]) / norm(Λ * m)
            err_ϕ_res = norm(Res_ϕ[idx_ice]) /sum(idx_ice)
            if (iter==1)  err_ϕ_ini = err_ϕ_res  end
            err_ϕ_resrel = err_ϕ_res/err_ϕ_ini # relative error

            # error for h

            # from Error - h
            Err_h .= abs.(Err_h .- h)
            err_h = norm(Err_h[idx_ice]) ./ norm(h0)
            #err_h = norm(Err_h[idx_ice]) ./ sum(idx_ice)
            if (iter==1)  err_h_ini = err_h  end
            err_h_rel = err_h/err_h_ini # relative error

            # from residual
            err_h_res   = norm(Res_h[idx_ice]) / norm(h0)
            #err_h_res   = norm(Res_h[idx_ice]) / sum(idx_ice)
            if (iter==1)  err_h_ini = err_h_res  end
            err_h_resrel = err_h_res/err_h_ini # relative error

            # decide which errors should be below the tolerance and be printed out
            err_ϕ_tol, err_h_tol = err_ϕ_res, err_h_res

            iter += 1

            if mod(iter, printit) == 0
                @printf("iterations = %d, error ϕ = %1.2e, error h = %1.2e \n", iter, err_ϕ_tol, err_h_tol)
            end

            # save error evolution in vector
            append!(errs_ϕ, err_ϕ)
            append!(errs_h, err_h)
            append!(errs_ϕ_rel, err_ϕ_rel)
            append!(errs_h_rel, err_h_rel)
            append!(errs_ϕ_res, err_ϕ_res)
            append!(errs_h_res, err_h_res)
            append!(errs_ϕ_resrel, err_ϕ_resrel)
            append!(errs_h_resrel, err_h_resrel)

        end
        ittot += iter; it += 1; t += dt
        if mod(it, printtime) == 0
            @printf("time step = %d, number of iterations = %d \n", it, iter)
        end

        ϕ_old .= ϕ
        h_old .= h
    end

    # give the effective pressure as output
    N = calc_N.(ϕ, ρi, ρw, g, H, zb)

    return model_output(N=N, ϕ=ϕ, h=h, qx=qx, qy=qy, qx_ice=qx_ice, qy_ice=qy_ice,
            Err_ϕ=Err_ϕ, Err_h=Err_h, Res_ϕ=Res_ϕ, Res_h=Res_h,
            ittot=ittot,
            errs_ϕ=errs_ϕ, errs_h=errs_h,
            errs_ϕ_rel=errs_ϕ_rel, errs_h_rel=errs_h_rel,
            errs_ϕ_res=errs_ϕ_res, errs_h_res=errs_h_res,
            errs_ϕ_resrel=errs_ϕ_resrel, errs_h_resrel=errs_h_resrel)
end