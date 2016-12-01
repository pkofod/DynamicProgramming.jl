"""
    bellman!(U, F, β, V)

Applies the Bellman operator. Uses the previous value function to generate
a new value function in-place.
"""
function bellman!(U, F, β, V::ValueFunction)
    copy!(V.Vᵏ,V.Vᵏ⁺¹)
    βEV!(β, F, V)
    V.Vᵏ⁺¹ .= max([U[i] + V.βEV[i] for i = 1:length(U)]...)
    supnorm(V)
end

# NOTE
# mapreduce(i->exp.(U.U[i] + V.βEV[i] - maximum(V)),+, 1:length(U.U)) is faster than
# V.Vᵏ⁺¹ .= sum(exp.(U[i] + V.βEV[i] - maximum(V)) for i = 1:length(U)), but less readable

function bellman!(U, F, β, V::IntegratedValueFunction)
    copy!(previous(V),current(V))
    get_maximum!(V)
    βEV!(β, F, V)
    V.Vᵏ⁺¹ .= sum(exp.(U[i] + V.βEV[i] - maximum(V)) for i = 1:length(U))
    map!(x->maximum(V)+log(x), V.Vᵏ⁺¹)
    supnorm(V)
end

function bellman!(U, F, β, EV::ExpectedValueFunction)
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
