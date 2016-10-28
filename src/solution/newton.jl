newton!(U, S, V::AbstractValueFunction, d) = newton!(U.U, U.β, U.P, S.F, V, d)
function newton!(U, β, P, F, V::AbstractValueFunction, d)
    bellman!(U, F, β, V)
    V_diff = V.Vᵏ⁺¹-V.Vᵏ
    P!(P, U, β, F, V.Vᵏ + d*V_diff)
    βFP!(V, β, P, F)
    V.Vᵏ⁺¹ .= V.Vᵏ + (I-V.βFP)\V_diff
end

newton!(U, S, EV::ExpectedValueFunction, d) = newton!(U.U, U.β, U.P, S.F, EV, d)
function newton!(U, β, P, F, EV::ExpectedValueFunction, d)
    nA, nX = length(U), length(U[1])
    bellman!(U, F, β, EV)
    EVᵏ⁺¹, EVᵏ = EV.EVᵏ⁺¹, EV.EVᵏ
    EV_diff = [EVᵏ⁺¹[ia]-EVᵏ[ia] for ia = 1:nA]
    P!(P, U, β, [EVᵏ[ia] + d*EV_diff[ia] for ia = 1:nA])
    βFP!(EV, P, β, F)
    ∂Γ∂EV = zeros(nA*size(F[1],1),nA*size(F[1],1))
    for ia = 1:nA
        for ia′ = 1:nA
            ∂Γ∂EV[(nX*(ia-1)+1):nX*ia, (nX*(ia′-1)+1):nX*ia′] = F[ia].*P[ia′]'
        end
    end
    ev = vcat(EVᵏ...) + (I-β*∂Γ∂EV)\vcat(EV_diff...)
    for ia = 1:nA
        EV.EVᵏ⁺¹[ia] .= ev[(nX*(ia-1)+1):nX*ia]
    end
end
