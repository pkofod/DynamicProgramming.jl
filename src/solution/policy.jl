"""
    policy!(U, β, P, F, V)
Applies a policy iteration.
"""
policy!{T<:AbstractValueFunction}(U, β, P, F, V::T) = error("policy! is not implemented for $(T).")
policy!(U, S, V) = policy!(U.U, U.β, U.P, S.F, V::IntegratedValueFunction)
function policy!(U, β, P, F, V::IntegratedValueFunction)
    copy!(previous(V),current(V))
    P!(P, U, β, F, current(V))
    βFP!(V, β, P, F)
    current(V) .= (I-V.βFP)\sum(P[a].*(U[a]-log.(P[a])) for a in 1:length(P))
end
