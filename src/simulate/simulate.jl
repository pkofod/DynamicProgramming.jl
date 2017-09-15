"""
    simulate(X, P, M, N, T, p1)
"""
function simulate(P::NTuple, X::DiscreteStates, N::Integer, T, p1)
    NT = N*T
    id = kron(collect(1:1:N), ones(Int64, T))
    a = fill(1, NT)
    x = fill(0, NT)
    xs = zeros(Int64, NT, length(X.Fs))
    for i = 1:T:NT
        for j = 1:length(X.X)
            xj = X.X[j]
            nj = length(xj)
            if i == 1 && j in X.market
                temp = draw(ones(nj)/nj)
                xs[1:T:(NT), j] = findfirst(temp)
            elseif !(j in X.market)
                temp = draw(ones(nj)/nj)
                xs[i, j] = findfirst(temp)
            end
        end
    end
    for i in 1:length(X.X)
        if i in X.market
            for t = 1:T-1
                xs[(1:T:(NT))+t, i] = draw(X.Fs[i][1][xs[t,i],:])
            end
        end
    end

    for i = 1:NT
        x[i]  = X.val_to_ind[(xs[i, :]...)]
        a[i] = draw([P[ip][x[i]] for ip = 1:length(P)])
        if i%T != 0
            for j = 1:length(X.X)
                if !(j in X.market)
                    xs[i+1, j] = draw(X.Fs[j][a[i]][xs[i,j],:])
                end
            end
        end
    end
    id, a, xs, x
end

function simulate(P::Tp, X::DiscreteState, N::Integer, T::Integer, p1::AbstractArray) where {Tp<:NTuple}
    NT = N*T
    id = kron(collect(1:1:N), ones(Int64, T))
    a = fill(1, NT)
    x = fill(0, NT)
    xs = zeros(Int64, NT, 1)

    for i = 1:T:NT
        temp = rand(Multinomial(1, ones(X.nX)/X.nX))
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
function simulate(U::AbstractUtility, X::AbstractState, M::Integer, N::Integer, T::Integer; maxiters = 10_000)
    V, iter = solve!(U, X, Poly())
    p0 = stationary(policy(U), X.F, maxiters = maxiters) #API solP into CCP(::sol)
    ds = [simulate(policy(U), X, N, T, p0)...  ones(Int64, N*T)]
    if M > 1
        for m = 2:M
            ds = [ds; [simulate(policy(U), X, N, T, p0)... fill(Int64(m), N*T)]]
        end
    end
    Data(ds[:,1], ds[:, 2], ds[:, end-1], ds[:, 3:end-2], ds[:, end], M*N*T)
end

function simulate(U::AbstractUtility, X, M::Integer, N::Integer, T::Integer, p0)
    sum(p0) == 1.0 || throw(error("Initial distribution did not sum to $(sum(p0)), not 1.0."))
    length(p0) == X.nX || throw(error("Initial distribution was of different dimension ($(length(p0))) than the state space ($(X.nX))."))
    V, iter = solve!(U, X)
    simulate(policy(U), X, M, N, T, p0)
end

function stationary(P, F; maxiters = 10_000)
    T = eltype(F[1])
    if T<:AbstractSparseMatrix
        Fᵁ = mapreduce(a->P[a].*F[a], +, 1:length(P))
    else
        Fᵁ = sum(P[ia].*F[ia] for ia = 1:length(P))
    end
    n = size(F[1], 1)
    πᵏ = fill(1/n, n)
    πˡ = fill(1/n, n)
    # l = k + 1
    for i = 1:10:maxiters
        copy!(πᵏ, πˡ)
        πˡ .= Fᵁ*πᵏ
        if norm(πᵏ - πˡ, Inf) < 1e-12
            break
        end
        if i == maxiters
            throw(error("Did not reach convergence in equilibrium probability matrix. Current sequential error is $(norm(πᵏ - πˡ, Inf))"))
            # TODO add an alternative method to solve for the stationary distribution
        end
    end
    return πˡ
end

function draw(p)
    u = rand()
    c = zero(eltype(p))
    a = 0
    @inbounds while a < length(p)
        a += 1
        c += p[i]
        if c > u
            break
        end
    end
    return a
end
