include("SHMIP_cases.jl")
using KissMCMC, Printf

test_case = "A1"
nx, ny = 64, 32
ttot = 1000.0

function make_logpdf(test_case, nx, ny, ttot)
    return function(thetas)
        it = (nothing for i=1:100)
        damp_ϕ,  # damping parameter for ϕ update
        damp_h,  # damping parameter for h update
        tau_ϕ,   # scaling factor for dτ_ϕ
        tau_h =  # scaling factor for dτ_h
        thetas

        if any(thetas .< 0.0)
            return -Inf
        else
            try
                it = run_SHMIP(test_case, Nx=nx, Ny=ny, t_tot=ttot,
                                γ_ϕ=damp_ϕ, γ_h=damp_h, τ_ϕ_=tau_ϕ, τ_h_=tau_h);
            catch y
                if isa(y, DomainError)
                    return -Inf
                end
            end
            if it >= 10^3
                return - Inf
            else
                return - it
            end
        end
    end
end

logpdf = make_logpdf(test_case, nx, ny, ttot)
theta0 = [0.1, 0.78, 1e6, 41.0]
std0   = [0.05, 0.1, 1e5, 4.0]

thetas, accept_ratioe, logdensities  = emcee(logpdf, make_theta0s(theta0, std0, logpdf, 6), niter=100)
thetas, accept_ratioe, logdensities  = squash_walkers(thetas, accept_ratioe, logdensities) # puts all walkers into one

damp_ϕ = [thetas[i][1] for i in 1:length(thetas)]
damp_h = [thetas[i][2] for i in 1:length(thetas)]
tau_ϕ  = [thetas[i][3] for i in 1:length(thetas)]
tau_h  = [thetas[i][4] for i in 1:length(thetas)]
iterations = -1 * logdensities

# compare a certain tau against number of iterations
scatter(damp_ϕ, iterations)

@printf("Minimal number of iterations = %d for parameters \n
         γ_ϕ  = %f, \n
         γ_h  = %f, \n
         τ_ϕ_ = %f, \n
         τ_h_ = %f",
         minimum(iterations),
         damp_ϕ[iterations .== minimum(iterations)][1],
         damp_h[iterations .== minimum(iterations)][1],
         tau_ϕ[iterations .== minimum(iterations)][1],
         tau_h[iterations .== minimum(iterations)][1]
         )
