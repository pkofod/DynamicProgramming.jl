newton!(U, S, V::AbstractValueFunction, d) = newton!(U.U, U.β, U.P, S.F, V, d)
function newton!(U, β, P, F, V::AbstractValueFunction, d)
    bellman!(U, F, β, V)
    V_diff = V.Vᵏ⁺¹-V.Vᵏ
    βEV!(V.βEV, β, F, V.Vᵏ + d*V_diff)
    P!(P, U, V.βEV)
    βFP!(V, β, P, F)
    V.Vᵏ⁺¹ .= V.Vᵏ + (I-V.βFP)\V_diff
end
