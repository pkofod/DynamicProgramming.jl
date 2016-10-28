"""
    bellman!(U, F, β, V)

Applies the Bellman operator. Uses the previous value function to generate
a new value function in-place.
"""
function bellman!(U, F, β, V::ValueFunction)
    copy!(previous(V),current(V))
    βEV!(β, F, V)
    V.Vᵏ⁺¹ .= max([U[i] + V.βEV[i] for i = 1:length(U)]...)
    supnorm(V)
end

bellman!(U, F, β, V::IntegratedValueFunction) = bellman!(current(V), previous(V), U, F, β, V.βEV)
# should have an evaluate only one
function bellman!(Vᵏ⁺¹, Vᵏ, U, F, β, βEV)
    copy!(Vᵏ, Vᵏ⁺¹)
    maxV = maximum(Vᵏ)
    βEV!(βEV, β, F, Vᵏ)
    Vᵏ⁺¹ .= maxV .+ log.(sum(exp.(U[i] .+ βEV[i] .- maxV) for i = 1:length(U)))
    supnorm(Vᵏ - Vᵏ⁺¹)
end
