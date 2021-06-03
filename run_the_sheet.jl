using Pkg, PyPlot, Printf, LinearAlgebra, Parameters
Pkg.activate(@__DIR__)
using SheetModel
S = SheetModel

# enable plotting by default
if !@isdefined do_visu; do_visu = true end
pygui(true)


@views function runthemodel(; do_visu=true)


    Params = Para()
    @unpack ev, g, ρw, Σ, Γ, Λ, m, lx, ly, nx, ny, dx, dy, xc, yc, dt, ttot, tol, itMax, damp, dτ = Params

    # Array allocation
    qx, qy, dϕ_dτ, dh_dτ, Res_ϕ, Res_h = S.array_allocation(Params)

    # Initial condition
    #ϕ0    = exp.(- 1e-2*(xc.-lx/2).^2) * exp.(-1e-2*(yc.-ly/2).^2)'
     ϕ0 = 1/lx * xc * ones(ny)'
    # ϕ0 = rand(nx, ny)
    ϕ_old   = copy(ϕ0)
    ϕ      = copy(ϕ0)

    h0 = 0.5 * ones(nx, ny)
    h_old = copy(h0)
    h = copy(h0)

    t = 0.0; it = 0; ittot = 0

    # Physical time loop
    while t<ttot
        iter = 0; err_ϕ = 2*tol; err_h = 2*tol;

        # Pseudo-transient iteration
        while max(err_ϕ, err_h) >tol && iter<itMax
            dϕ_dx, dϕ_dy = diff(ϕ, dims=1) ./ dx, diff(ϕ, dims=2) ./ dy   # hydraulic gradient
            qx, qy     = S.calc_q(h[1:end-1, :], dϕ_dx, Params), S.calc_q(h[:, 1:end-1], dϕ_dy, Params) # TODO: which h values to take?!

            vo     = S.calc_vo(h, Params)                   # opening rate
            vc     = S.calc_vc(ϕ, h, Params)  # closure rate

            # calculate residuals
            Res_ϕ        =    - ev/(ρw*g) * (ϕ[2:end-1, 2:end-1] .- ϕ_old[2:end-1, 2:end-1])/dt .-    # dhe/dt
                                (diff(qx, dims=1)[:, 2:end-1]/dx .+ diff(qy, dims=2)[2:end-1, :]/dy)    .+    # div(q)
                                Σ * vo[2:end-1, 2:end-1] .- Γ * vc[2:end-1, 2:end-1]            .-    # dh/dt
                                Λ * m                                                                    # source term
            Res_h          =    (h[2:end-1, 2:end-1] .- h_old[2:end-1, 2:end-1]) / dt  .-
                                Σ * vo[2:end-1, 2:end-1] .+ Γ * vc[2:end-1, 2:end-1]

            # damped rate of change
            dϕ_dτ      = Res_ϕ .+ damp * dϕ_dτ
            dh_dτ        = Res_h .+ damp * dh_dτ

            # update fields
            ϕ[2:end-1, 2:end-1]  = ϕ[2:end-1, 2:end-1] + dτ * dϕ_dτ   # update rule, sets the BC as ϕ[1]=ϕ[end]=0
            h[2:end-1, 2:end-1]    = h[2:end-1, 2:end-1] + dτ * dh_dτ

            err_ϕ = norm(Res_ϕ)/length(Res_ϕ)
            err_h   = norm(Res_h)/length(Res_h)
            iter += 1
        end
        ittot += iter; it += 1; t += dt

        ϕ_old .= ϕ
        h_old .= h
    end

    # Visualize
    if do_visu
        figure()
        subplot(1, 2, 1)
        pcolor(xc, yc, ϕ0')
        subplot(1, 2, 2)
        pcolor(xc, yc, ϕ')
    end

end

runthemodel(; do_visu=do_visu);