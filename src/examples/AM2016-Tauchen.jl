using StatsFuns
function am_tauchen(;β=0.95, N = 500)
    # Aguirregabiria & Magesan 2016
    srand(122)
    bdisc = 0.95 # discount factor
    truepar = [-1.9;1.;2.] # fixed costs, ..., entry costs
    # State 1: Exogenous state
    nX1 = 2
    X1 = [0.;1.]
    F1 = [[1. 0.; 1. 0.], [0. 1.; 0. 1.]]

    nX2 = N
    X2 = 1:nX2
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
    F2 = tauchen(0.1,0.01,2, nX2)


    S = States(State(X1, F1),
               CommonState(X2, F2))

    Z = [zeros(nX2*2, 3), # don't buy
                 [-ones(nX2) log(X2) -ones(nX2); # buy
    			  -ones(nX2) log(X2)  zeros(nX2)]]

    U = LinearUtility(Z, bdisc, truepar)

    return U, S
end
