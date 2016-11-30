abstract AbstractValueFunction
maximum(AVF::AbstractValueFunction) = AVF.maximum[1]

# Actual value function
immutable ValueFunction{T<:Real}<:AbstractValueFunction
    Vᵏ⁺¹::Vector{T}
    Vᵏ::Vector{T}
    βFP::Matrix{T}
    βEV::Vector{Vector{T}} # maybe even beta ev?
end
ValueFunction(nX, nA) = ValueFunction(zeros(nX), zeros(nX), zeros(nX, nX), [zeros(nX) for i = 1:nA])
supnorm(V::AbstractValueFunction) = norm(V.Vᵏ - V.Vᵏ⁺¹, Inf)
maximum(V::ValueFunction) = err("maximum not defined for ValueFunction.")

# Integrated value function
immutable IntegratedValueFunction{T<:Real}<:AbstractValueFunction
    Vᵏ⁺¹::Vector{T}
    Vᵏ::Vector{T}
    βFP::Matrix{T}
    βEV::Vector{Vector{T}}
    maximum::Vector{T}
end

IntegratedValueFunction(nX, nA) = IntegratedValueFunction(nX, nA, Float64)
IntegratedValueFunction(nX, nA, T) = IntegratedValueFunction(zeros(T, nX), zeros(T, nX), zeros(T, nX, nX), [zeros(T, nX) for i = 1:nA], [zero(T),])
get_maximum!(V::AbstractValueFunction) = copy!(V.maximum, maximum(V.Vᵏ))

# Expected value function
immutable ExpectedValueFunction{T<:Real}<:AbstractValueFunction
    EVᵏ::Vector{Vector{T}}
    EVᵏ⁺¹::Vector{Vector{T}}
    βFP::Vector{Matrix{T}}
    logsum::Vector{T}
    maximum::Vector{T}
end

ExpectedValueFunction(nX, nA) = ExpectedValueFunction(nX, nA, Float64)
ExpectedValueFunction(nX, nA, T) = ExpectedValueFunction([zeros(T, nX) for i = 1:nA],
                                                         [zeros(T, nX) for i = 1:nA],
                                                         [zeros(T, nX, nX) for i = 1:nA],
                                                         zeros(T, nX), [zero(T),])
supnorm(EV::ExpectedValueFunction) = maximum([norm(EV.EVᵏ[iA] - EV.EVᵏ⁺¹[iA], Inf) for iA = 1:length(EV.EVᵏ)])
get_maximum!(EV::ExpectedValueFunction) = copy!(EV.maximum, mapreduce(maximum, maximum, EV.EVᵏ))

# Relative value function
immutable RelativeValueFunction{T<:Real}<:AbstractValueFunction
    RVᵏ::Vector{Vector{T}}
    RVᵏ⁺¹::Vector{Vector{T}}
    βFP::Vector{Matrix{T}}
    logsum::Vector{T}
    maximum::Vector{T}
end
RelativeValueFunction(nX, nA) = RelativeValueFunction(nX, nA, Float64)
RelativeValueFunction(nX, nA, T) = RelativeValueFunction([zeros(T, nX) for i = 1:nA],
                                                         [zeros(T, nX) for i = 1:nA],
                                                         [zeros(T, nX, nX) for i = 1:nA],
                                                         zeros(T, nX), [zero(T),])
supnorm(RV::RelativeValueFunction) = maximum([norm(RV.RVᵏ[iA] - RV.RVᵏ⁺¹[iA], Inf) for iA = 1:length(RV.RVᵏ)])
get_maximum!(RV::RelativeValueFunction) = copy!(RV.maximum, mapreduce(maximum, max, RV.RVᵏ))
