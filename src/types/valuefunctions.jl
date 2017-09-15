abstract type AbstractValueFunction end
difference(V::AbstractValueFunction) = current(V)-previous(V)
maximum(AVF::AbstractValueFunction) = AVF.maximum[1]

# Actual value function
struct ValueFunction{Tv, Tm, Tvv}<:AbstractValueFunction
    Vᵏ⁺¹::Tv
    Vᵏ::Tv
    βFP::Tm
    βEV::Tvv
end
ValueFunction(nX, nA, T=Float64) = ValueFunction(zeros(T, nX), zeros(T, nX), zeros(T, nX, nX), [zeros(T, nX) for i = 1:nA])
ValueFunction(S) = ValueFunction(S.nX, length(S.F))
supnorm(V) = Base.LinAlg.generic_vecnormInf(V)
supnorm(V::AbstractValueFunction) = Base.LinAlg.generic_vecnormInf(V.Vᵏ - V.Vᵏ⁺¹)
maximum(V::ValueFunction) = err("maximum not defined for ValueFunction.")

# Integrated value function
mutable struct IntegratedValueFunction{Tv, Tm, Tvv} <: AbstractValueFunction
    Vᵏ⁺¹::Tv
    Vᵏ::Tv
    βFP::Tm
    βEV::Tvv
end
IntegratedValueFunction(nX, nA, T=Float64) = IntegratedValueFunction(zeros(T, nX), zeros(T, nX), zeros(T, nX, nX), [zeros(T, nX) for i = 1:nA])
function IntegratedValueFunction(S::AbstractState, T=Float64)
    nX, nA = S.nX, length(S.F)
    βFP = similar(S.F[1])
    IntegratedValueFunction(zeros(T, nX), zeros(T, nX), βFP, [zeros(T, nX) for i = 1:nA])
end
get_maximum!(V::AbstractValueFunction) = copy!(V.maximum, maximum(V.Vᵏ))
size(V::IntegratedValueFunction) = size(V.Vᵏ)
