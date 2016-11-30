# In-place Discounted Expected Value Function Calculation
function βEV!{T<:Real}(β::T, F, V)
    for i = 1:length(F)
        A_mul_B!(β, F[i], V.Vᵏ, zero(T), V.βEV[i])
    end
end

function βEV!{T<:Real}(β::T, F, EV::ExpectedValueFunction)
    for iA = 1:length(F)
        A_mul_B!(β, F[iA], EV.EVᵏ, zero(T), EV.βEV[iA])
    end
    mapreduce(maximum, max, EV.EVᵏ)
end

# In-place Bellman operators
function Γ!(U, F, β, V::ValueFunction)
    copy!(V.Vᵏ,V.Vᵏ⁺¹)
    βEV!(β, F, V)
    V.Vᵏ⁺¹ .= maximum([U[i] + V.βEV[i] for i = 1:length(U)]...)
    supnorm(V)
end

# NOTE
# mapreduce(i->exp.(U.U[i] + V.βEV[i] - maximum(V)),+, 1:length(U.U)) is faster than
# V.Vᵏ⁺¹ .= sum(exp.(U[i] + V.βEV[i] - maximum(V)) for i = 1:length(U)), but less readable

function Γ!(U, F, β, V::IntegratedValueFunction)
    copy!(V.Vᵏ,V.Vᵏ⁺¹)
    get_maximum!(V)
    βEV!(β, F, V)
    V.Vᵏ⁺¹ .= sum(exp.(U[i] + V.βEV[i] - maximum(V)) for i = 1:length(U))
    map!(x->maximum(V)+log(x), V.Vᵏ⁺¹)
    supnorm(V)
end

function Γ!(U, F, β, EV::ExpectedValueFunction)
    for i = 1:length(EV.EVᵏ⁺¹)
        copy!(EV.EVᵏ[i],EV.EVᵏ⁺¹[i])
    end
    get_maximum!(EV)
    EV.logsum .= sum(exp.(U[iA] + β*EV.EVᵏ[iA] - maximum(EV)) for iA = 1:length(U))
    map!(log, EV.logsum)
    for iA = 1:length(U)
        EV.EVᵏ⁺¹[iA].=maximum(EV)
        A_mul_B!(1.0, F[iA], EV.logsum, 1.0, EV.EVᵏ⁺¹[iA])
    end
    supnorm(EV)
end

function Γ!(U, F, β, RV::RelativeValueFunction)
    for i = 1:length(RV.RVᵏ⁺¹)
        copy!(RV.RVᵏ[i],RV.RVᵏ⁺¹[i])
    end
    get_maximum!(RV)
    RV.logsum .= sum(exp.(U[iA] + β*RV.RVᵏ[iA] - maximum(RV)) for iA = 1:length(U))
    map!(log, RV.logsum)
    for iA = 1:length(U)
        RV.RVᵏ⁺¹[iA].=maximum(RV)
        A_mul_B!(1.0, F[iA], RV.logsum, 1.0, RV.RVᵏ⁺¹[iA])
    end
    supnorm(RV)
end
