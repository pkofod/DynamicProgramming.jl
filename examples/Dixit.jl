function dixit(;N = 5, β = 0.99)
    # State 1
    X1 = 0:1
    F1 = [[1. 0.; 1. 0.], [0. 1.; 0. 1.]]
    # State 2
    nX2 = N
    X2 = 1:nX2
    F2 = 1./(1+abs(ones(length(X2),1)*X2'-X2*ones(1, length(X2))))
    F2 = F2./sum(F2,1)'
    # States
    S = States(State(X1, F1), CommonState(X2, F2))

    # Utility variables
    Z1 = zeros(nX2*2, 3)
    Z2 = [ones(nX2) X2 -ones(nX2);
          ones(nX2) X2 zeros(nX2)]
    U = LinearUtility([Z1, Z2], β, [-.1;.2;1])

    return U, S
end
