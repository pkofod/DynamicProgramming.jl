abstract AbstractApproximation
abstract Approximation

immutable Sieve{T<:Real} <: Approximation
    η::Vector{T} # Basis function
    ξ::Vector{T} # Evaluation grid
    ω::Vector{T} # Weights
    o::Integer # Order
    TpT::Matrix{T} # Projection matrix
end

abstract Chebyshev
function Chebyshev(order)
    η = nodes(Chebyshev, order)
    ω = weights(Chebyshev, order)
    TpT = TpT(Chebyshev, order)
    Sieve(η, ω, order,TpT)
end

function nodes(t::Type{Chebyshev}, order)
    generate chebyshev on unigrid
    return nodes
end

function weights(t::Type{Chebyshev}, order)
    caluclate widhgts
    return weights
end

function TpT(t::Type{Chebyshev}, order)
    generate projection
    return evaluation at nodes
