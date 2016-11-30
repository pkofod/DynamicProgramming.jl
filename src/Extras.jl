module Extras
using StatsFuns
function tauchen(ρ, σₛ, m, N)
    const Φ = normcdf # CDF of standard normal
    s̃₁, s̃ₙ = -m*σₛ, m*σₛ # end points
    s̃ = linspace(s̃₁, s̃ₙ, N) # grid
    w = (s̃[2]-s̃[1])/2 # half distance between grid points
    F = zeros(N, N) # empty transition matrix
    F[:, 1] = Φ.((s̃[1]-ρ.*s̃+w)/sqrt(σₛ))
    F[:, N] = 1-Φ.((s̃[end]-ρ.*s̃-w)/sqrt(σₛ))
    for j = 2:N-1
        for i = 1:N
            F[i, j] = Φ.((s̃[j]-ρ*s̃[i]+w)/sqrt(σₛ))-Φ.((s̃[j]-ρ*s̃[i]-w)/sqrt(σₛ))
        end
    end
    F
end
end
