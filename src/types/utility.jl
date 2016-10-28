abstract AbstractUtility
type LinearUtility{T<:Real} <: AbstractUtility
    Z::Vector{Matrix{T}}
    β::T
    Θ::Vector{T}
    U::Vector{Vector{T}}
    P::Vector{Vector{T}}
    nvar::Int64
    names::Vector
end

Base.copy(m::LinearUtility) = LinearUtility([ copy(getfield(m, k)) for k = 1:length(fieldnames(m)) ]...)
_makenames(Z) = Symbol[Symbol("choice$i") for i = 1:length(Z)]
LinearUtility{T<:Real}(Z::Vector{Vector{T}}, β::T) =
    LinearUtility(Z, β, _makenames(Z))
function LinearUtility{T<:Real}(Z::Vector{Vector{T}}, β::T, names)
    nrow, ncol = size(Z[1])
    LinearUtility(Z, β, zeros(T, ncol),
                        [zeros(T, nrow) for i = 1:length(Z)],
                        [zeros(T, nrow) for i = 1:length(Z)],
                        ncol,
                        names)
end

LinearUtility{T<:Real}(Z, β::T, Θ::Vector{T}) =
    LinearUtility(Z, β, Θ, _makenames(Z))
function LinearUtility{T<:Real}(names, Z, β::T, Θ::Vector{T})
    U = [Z[i]*Θ for i in 1:length(Z)]
    P = [zeros(T, size(Z[iA], 1)) for iA = 1:length(Z)]

    LinearUtility(Z, β, Θ, U, P, size(Z[1], 2), names)
end
