abstract type AbstractValueFunction end
difference(V::AbstractValueFunction) = current(V)-previous(V)
maximum(AVF::AbstractValueFunction) = AVF.maximum[1]

# Actual value function
immutable ValueFunction{T<:Real}<:AbstractValueFunction
    Vᵏ⁺¹::Vector{T}
    Vᵏ::Vector{T}
    βFP::Matrix{T}
    βEV::Vector{Vector{T}} # maybe even beta ev?
end
ValueFunction(nX, nA, T=Float64) = ValueFunction(zeros(T, nX), zeros(T, nX), zeros(T, nX, nX), [zeros(T, nX) for i = 1:nA])
ValueFunction(S) = ValueFunction(S.nX, length(S.F))
supnorm(V) = Base.LinAlg.generic_vecnormInf(V)
supnorm(V::AbstractValueFunction) = Base.LinAlg.generic_vecnormInf(V.Vᵏ - V.Vᵏ⁺¹)
maximum(V::ValueFunction) = err("maximum not defined for ValueFunction.")

# Integrated value function
mutable struct IntegratedValueFunction{T<:Real, Tm<:AbstractMatrix}<:AbstractValueFunction
    Vᵏ⁺¹::Vector{T}
    Vᵏ::Vector{T}
    βFP::Tm
    βEV::Vector{Vector{T}}
end
IntegratedValueFunction(nX, nA, T=Float64) = IntegratedValueFunction(zeros(T, nX), zeros(T, nX), zeros(T, nX, nX), [zeros(T, nX) for i = 1:nA])
function IntegratedValueFunction(S::AbstractState, T=Float64)
    nX, nA = S.nX, length(S.F)
    βFP = similar(S.F[1])
    IntegratedValueFunction(zeros(T, nX), zeros(T, nX), βFP, [zeros(T, nX) for i = 1:nA])
end
get_maximum!(V::AbstractValueFunction) = copy!(V.maximum, maximum(V.Vᵏ))
size(V::IntegratedValueFunction) = size(V.Vᵏ)
#=
# Expected value function
immutable ExpectedValueFunction{T<:Real}<:AbstractValueFunction
    EVᵏ::Vector{Vector{T}}
    EVᵏ⁺¹::Vector{Vector{T}}
    βFP::Vector{Matrix{T}}
    logsum::Vector{T}
    maximum::Vector{T}
end

ExpectedValueFunction(nX, nA, T=Float64) = ExpectedValueFunction([zeros(T, nX) for i = 1:nA],
                                                         [zeros(T, nX) for i = 1:nA],
                                                         [zeros(T, nX, nX) for i = 1:nA],
                                                         zeros(T, nX), [zero(T),])
ExpectedValueFunction(S) = ExpectedValueFunction(S.nX, length(S.F))
get_maximum!(EV::ExpectedValueFunction) = copy!(EV.maximum, mapreduce(maximum, max, EV.EVᵏ))
supnorm(EV::ExpectedValueFunction) = maximum([Base.LinAlg.generic_vecnormInf(EV.EVᵏ[iA] - EV.EVᵏ⁺¹[iA]) for iA = 1:length(EV.EVᵏ)])
=#
