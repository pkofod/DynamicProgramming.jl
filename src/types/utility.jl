abstract AbstractUtility
type LinearUtility{T<:Real} <: AbstractUtility
    Z::Vector{Matrix{T}}
    β::T
    θ::Vector{T}
    U::Vector{Vector{T}}
    P::Vector{Vector{T}}
    nvar::Int64
end

Base.copy(m::LinearUtility) = LinearUtility([ copy(getfield(m, k)) for k = 1:length(fieldnames(m)) ]...)

function LinearUtility{T<:Real}(Z::Vector{Vector{T}}, β::T)
    nrow, ncol = size(Z[1])
    LinearUtility(Z, β, zeros(T, ncol),
                        [zeros(T, nrow) for i = 1:length(Z)],
                        [zeros(T, nrow) for i = 1:length(Z)], ncol)
end

function LinearUtility{T<:Real}(Z, β::T, θ::Vector{T})
    U = [Z[i]*θ for i in 1:length(Z)]
    P = [zeros(T, size(Z[iA], 1)) for iA = 1:length(Z)]

    LinearUtility(Z, β, θ, U, P, size(Z[1], 2))
end
