using MDPTools
function hitsch(;N = 20)
    n = [0, 2, 5]
    next_i(i, k) = max(0, min(I, i+n[k]-1))
    F_X = ["low", "high"]
    F_P = [0.16 0.84; 0.16 0.84]
    SP = CommonState(F_X, F_P)

    I = N
    K = length(n)
    F_i = [spzeros(1+I, 1+I) for k = 1:K]
    for i = 0:I, j = 1:K
            F_i[j][i+1, next_i(i, j)+1]=1
    end
    Si = State(0:I, F_i);

    S = States(SP, Si)

    delta, alpha = 4., 4.

    Z1 = [[zeros(1); ones(I)] zeros(I+1) -next_i(0:I, 1);
          [zeros(1); ones(I)] zeros(I+1) -next_i(0:I, 1)]
    Z2 = [ones(2*(I+1)) -[1.2*ones(I+1);2*ones(I+1)] -kron(ones(2), next_i(0:I, 2))]
    Z3 = [ones(2*(I+1)) -[3.0*ones(I+1);5*ones(I+1)] -kron(ones(2), next_i(0:I, 3))]
    U = LinearUtility(["Buy 0", "Buy 2", "Buy 5"], [Z1, Z2, Z3], 0.998, [delta; alpha; 0.05])
    return U, S
end
U, S = hitsch()
V, iter = solve!(U, S)<s
