policy(U::AbstractUtility) = U.P
policy(U::AbstractUtility, a) = U.P[a]

"""
Returns the current value function iterate.
"""
current(V::AbstractValueFunction) = V.Vᵏ⁺¹

"""
Returns the previous value function iterate.
"""
previous(V::AbstractValueFunction) = V.Vᵏ

"""
    βEV!(β, F, V)
Calculates the discounted expected value in-place. This function is used in the
bellman operator.
"""
function βEV!{T<:Real}(βEV, β::T, F, V)
    for i = 1:length(F)
        A_mul_B!(β, F[i], V, zero(T), βEV[i])
    end
end

"""
    βFP!(V, β, P, F)
Calculates the discounted unconditional transition probabilities induced by P.
"""
βFP!{TM <: AbstractMatrix}(V::IntegratedValueFunction, β, P, F::Vector{TM}) = βFP!(V.βFP, β, F, P)
βFP!(βFP::AbstractMatrix, β, F, P) = _βFP!(βFP, β, F, P)
function βFP!{T}(βFP::AbstractMatrix, β::T, P::Vector{Vector{T}}, F::Vector{Matrix{T}})
    # For small problems or many choices, this unrolled loop seems to be faster than the fallback _βFP!
    if nX < 500 || length(P) > 10
        nX = length(P[1])
        fill!(βFP, zero(T))
        @inbounds for j = eachindex(P)
            bpj = β*P[j]
            Fj = F[j]
            for x′ = 1:nX
                for x = 1:nX
                    βFP[x, x′] += bpj[x]*Fj[x, x′]
                end
            end
        end
    else
        _βFP!(βFP, β, F, P)
    end
end

function _βFP!(βFP, β, F, P)
    J = length(P)
    βFP .= β*P[1].*F[1]
    for a = 2:J
        βFP .+= β*P[a].*F[a]
    end
end
# For sparse matrices this is faster
function βFP!{T<:AbstractSparseMatrix}(V::IntegratedValueFunction, β, P, F::Vector{T})
    V.βFP = mapreduce(a->β*P[a].*F[a],+, eachindex(P))
end

"""
    P!(args...)
Updates the choice probabilities given the provided value function.
"""
P!(U, F, V) = P!(U.P, U.U, U.β, F, V)
function P!(P, U, β, F, V)
    if length(P) == 2
        P[1] .= 1./(1 .+ exp.(U[2] .+ β*(F[2]*V) .- (U[1] .+ β*(F[1]*V))))
        P[2] .= 1 .- P[1]
    else
        _maximum = maximum(V)
        denominator = sum(exp.(U[a]+ β*(F[a]*V)-_maximum) for a = 1:length(P))
        for a = 1:length(P)-1
            P[a][:] = exp.(U[a] .+ β*(F[a]*V).-_maximum)./denominator
        end
        P[end][:] = 1 - sum(P[a] for a = 1:length(P)-1)
    end
end

P!(U, EV) = P!(U.P, U.U, U.β, EV)
P!{Tm<:AbstractVector}(P, U::Vector, β::Real, EV::Vector{Tm}) = P!(P, U, map(ev->β*ev,EV))
function P!{Tm<:AbstractVector}(P, U::Vector, βEV::Vector{Tm})
    if length(P) == 2
        P[1][:] = 1./(1 + exp.(U[2] + βEV[2] - (U[1] + βEV[1])))
        P[2][:] = 1 - P[1]
    else
        _maximum = mapreduce(maximum, max, βEV)
        denominator = sum(exp.(U[a]+βEV[a]-_maximum) for a = 1:length(P))
        for a = 1:length(P)-1
            P[a][:] = exp.(U[a]+βEV[a]-_maximum)./denominator
        end
        P[end][:] = 1 - sum(P[a] for a = 1:length(P)-1)
    end
end

# These are just short-hand notation for the method above. Makes it easier to
# update U.P with a U, S, and value function.
P!{Ts<:AbstractState}(U, S::Ts, V) = P!(U, S.F, current(V))
