using Stocas
function am_tauchen(;β=0.95, # discount factor,
                     N = 500)
    # Aguirregabiria & Magesan 2016
    truepar = [-1.9;1.;2.] # fixed costs, ..., entry costs

    nX2 = N
    X2 = 1:nX2
    F2 = Extras.tauchen(0.1,0.01,2, nX2)

    S = States(EntryState(),
               CommonState(X2, F2))

    Z = [zeros(nX2*2, 3), # don't buy
        [-ones(nX2) log(X2) -ones(nX2); # buy
         -ones(nX2) log(X2)  zeros(nX2)]]

    U = LinearUtility(Z, β, truepar)

    return U, S
end
solve!(am_tauchen()...)
