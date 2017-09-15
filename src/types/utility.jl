abstract type AbstractUtility end
struct LinearUtility{N, T} <: AbstractUtility
    names::NTuple{N, String}
    Z::NTuple{N, Matrix{T}}
    β::T
    Θ::Vector{T}
    U::NTuple{N, Vector{T}}
    P::NTuple{N, Vector{T}}
    nvar::Int
end

Base.copy(m::LinearUtility) = LinearUtility([ copy(getfield(m, k)) for k = 1:length(fieldnames(m)) ]...)
_makenames{N,T}(Z::NTuple{N, T}) = Tuple("choice $i" for i = 1:N)

LinearUtility{N, T<:Real}(Z::NTuple{N, AbstractMatrix{T}}, β::T) = LinearUtility(_makenames(Z), Z, β)
LinearUtility{N, T<:Real}(Z::NTuple{N, AbstractMatrix{T}}, β::T, Θ::Vector{T}) = LinearUtility(_makenames(Z), Z, β, Θ)

function LinearUtility{N, T}(names, Z::NTuple{N, AbstractMatrix{T}}, β::T)
    nrow, ncol = size(Z[1])
    LinearUtility(names, Z, β, zeros(T, ncol),
                        Tuple(zeros(T, nrow) for i = 1:N),
                        Tuple(zeros(T, nrow) for i = 1:N),
                        ncol)
end

function LinearUtility{N, T}(names, Z::NTuple{N, AbstractMatrix{T}}, β::T, Θ::Vector{T})
    U = Tuple(Z[i]*Θ for i = 1:N)
    P = Tuple(zeros(T, size(Z[iA], 1)) for iA = 1:N)

    LinearUtility(names, Z, β, Θ, U, P, size(Z[1], 2))
end

function Base.display(U::LinearUtility)
    @printf "Linear utility\n"
    @printf " * Number of choices: %s\n" length(U.U)
    for ix = 1:length(U.U)
        @printf "   %d) %s\n" ix U.names[ix]
    end
    @printf " * Number of variables:  %s\n" U.nvar
    @printf " * Discount rate:        %s\n" U.β
end
