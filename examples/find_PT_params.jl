include("SHMIP_cases.jl")
using KissMCMC, Printf, StatsPlots

test_case = "A5"
nx, ny = 64, 32

function make_logpdf(test_case, nx, ny)
    return function(thetas)
        it = (nothing for i=1:100)
        damp_ϕ,  # damping parameter for ϕ update
        damp_h,  # damping parameter for h update
        dtau_ϕ,   # scaling factor for dτ_ϕ
        dtau_h =  # scaling factor for dτ_h
        thetas

        if 0.5 < damp_ϕ < 1.0 && 0.5 < damp_h < 1.0 && 0.0 < dtau_ϕ < 5.0 && 0.0 < dtau_h < 5.0
            try
                it, ϕ = run_SHMIP(test_case, Nx=nx, Ny=ny,
                                γ_ϕ=damp_ϕ, γ_h=damp_h, dτ_ϕ_=dtau_ϕ, dτ_h_=dtau_h);
            catch y
                if isa(y, DomainError)
                    return typemin(Int)
                end
            end
            if it >= 2*10^4
                return typemin(Int)
            else
                return - it/100.0
            end
        else
            return typemin(Int)
        end
    end
end

logpdf = make_logpdf(test_case, nx, ny)
theta0 = [0.9, 0.9, 1.0, 1.0]
std0   = [0.05, 0.05, 0.1, 0.1]

thetas, accept_ratioe, logdensities  = emcee(logpdf, make_theta0s(theta0, std0, logpdf, 6), niter=100)
thetas, accept_ratioe, logdensities  = squash_walkers(thetas, accept_ratioe, logdensities) # puts all walkers into one

damp_ϕ = [thetas[i][1] for i in 1:length(thetas)]
damp_h = [thetas[i][2] for i in 1:length(thetas)]
dtau_ϕ  = [thetas[i][3] for i in 1:length(thetas)]
dtau_h  = [thetas[i][4] for i in 1:length(thetas)]
iterations = - logdensities
iterations[iterations .== typemin(Int)] .= typemax(Int) # typemin() are not affected by the minus sign one line above

@printf("Minimal number of iterations = %d for parameters \n
         γ_ϕ  = %f,
         γ_h  = %f,
         dτ_ϕ_ = %f,
         dτ_h_ = %f",
         minimum(iterations),
         damp_ϕ[iterations .== minimum(iterations)][1],
         damp_h[iterations .== minimum(iterations)][1],
         dtau_ϕ[iterations .== minimum(iterations)][1],
         dtau_h[iterations .== minimum(iterations)][1]
         )

# compare a certain tau against number of iterations
gr(size=(1300,700));
cornerplot([damp_ϕ damp_h dtau_ϕ dtau_h], label = ["γ_ϕ" "γ_h" "dτ_ϕ_" "dτ_h_"], leftmargin = 0.8Plots.cm, bottommargin = 0.8Plots.cm)