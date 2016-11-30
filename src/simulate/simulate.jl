"""
    simulate(X, P, M, N, T, p1)
"""
function simulate(P::AbstractVector, X::DiscreteStates, N, T, p1)
    NT = N*T
    id = kron(collect(1:1:N), ones(Int64, T))
    a = fill(1, NT)
    x = fill(0, NT)
    xs = zeros(Int64, NT, length(X.Fs))
    for i = 1:T:NT
        for j = 1:length(X.X)
            if i == 1 && j in X.market
                temp = rand(Multinomial(1, ones(X.X[j])/length(X.X[j])))
                xs[1:T:(NT), j] = findfirst(temp)
            elseif !(j in X.market)
                temp = rand(Multinomial(1, ones(X.X[j])/length(X.X[j])))
                xs[i, j] = findfirst(temp)
            end
        end
    end
    for i in 1:length(X.X)
        if i in X.market
            for t = 1:T-1
                # this should be a draw function and I should use foreach(draw, 1:T-1)
                xs[(1:T:(NT))+t, i] = findfirst(rand(Multinomial(1, vec(full(X.Fs[i][1][xs[t,i],:])))))
            end
        end
    end

    for i = 1:NT
        a[i] = findfirst(rand(Multinomial(1, [P[ip][X.val_to_ind[xs[i, :]...]] for ip = 1:length(P)])))
        if i%T != 0
            for j = 1:length(X.X)
                if !(j in X.market)
                    xs[i+1, j] = findfirst(rand(Multinomial(1, vec(full(X.Fs[j][a[i]][xs[i,j],:])))))
                end
            end
        end
    end
    for i = 1:NT
        x[i]  = X.val_to_ind[(xs[i, :]...)]
    end
    id, a, xs, x
end

function simulate(P::AbstractVector, X::DiscreteState, N, T, p1)
    NT = N*T
    id = kron(collect(1:1:N), ones(Int64, T))
    a = fill(1, NT)
    x = fill(0, NT)
    xs = zeros(Int64, NT, 1)
    for i = 1:T:NT
        temp = rand(Multinomial(1, ones(X.X)/length(X.X)))
        x[i] = findfirst(temp)
    end
    for i = 1:NT
        if P[2][x[i]] > rand()
            a[i] = 2
        end
        if i%T != 0
            mm=Multinomial(1, vec(full(X.F[a[i]][x[i],:])))
            x[i+1] = findfirst(rand(mm))
        end
    end
    id, a, xs, x
end

# TODO populate x, and remove ds
function simulate(U::AbstractUtility, X, M, N, T; max_iterations = 10_000)
    V, iter = solve(U, X, Poly())
    p0 = stationary(policy(U), X.F, max_iterations = max_iterations) #API solP into CCP(::sol)
    ds = [simulate(policy(U), X, N, T, p0)...  ones(Int64, N*T)]
    if M > 1
        for m = 2:M
            ds = [ds; [simulate(policy(U), X, N, T, p0)... fill(Int64(m), N*T)]]
        end
    end
    Data(ds[:,1], ds[:, 2], ds[:, end-1], ds[:, 3:end-2], ds[:, end], M*N*T)
end

function simulate(U::AbstractUtility, X, M, N, T, p0)
    sum(p0) == 1.0 || throw(error("Initial distribution did not sum to $(sum(p0)), not 1.0."))
    length(p0) == X.nX || throw(error("Initial distribution was of different dimension ($(length(p0))) than the state space ($(X.nX))."))
    V, iter = solve(U, X)
    simulate(policy(U), X, M, N, T, p0)
end

function stationary(P, F; max_iterations = 10_000)
    Fᵁ =  sum(P[ia].*F[ia] for ia = 1:length(P))
    Eq0 = copy(Fᵁ)
    Eq1 = Eq0^2
    for i = 1:10:max_iterations
        copy!(Eq0, Eq1)
        Eq1[:] = Fᵁ^i
        if norm(Eq0-Eq1, Inf) < 1e-12
            break
        end
        if i == max_iterations
            throw(error("Did not reach convergence in equilibrium probability matrix. Current sequential error is $(norm(Eq0-Eq1, Inf))"))
            # TODO add an alternative method to solve for the stationary distribution
        end
    end
    full(Eq1[1,:])
end
