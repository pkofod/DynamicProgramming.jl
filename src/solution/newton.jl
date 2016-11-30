function P!(U, F, V)
    if length(U.P) == 2
        U.P[1][:] = 1./(1 + exp.(U.U[2] + U.β*(F[2]*V) - (U.U[1] + U.β*(F[1]*V))))
        U.P[2][:] = 1 - U.P[1]
    else
        _maximum = maximum(V)
        denominator = sum(exp.(U.U[ia]+U.β*(F[ia]*V)-_maximum) for ia = 1:length(U.P))
        for ia = 1:length(U.P)-1
            U.P[ia][:] = exp.(U.U[ia]+U.β*(F[ia]*V)-_maximum)./denominator
        end
        U.P[end][:] = 1 - sum(U.P[ia] for ia = 1:length(U.P)-1)
    end
end

function P!(U, EV)
    if length(U.P) == 2
        U.P[1][:] = 1./(1 + exp.(U.U[2] + U.β*(EV[2]) - (U.U[1] + U.β*(EV[1]))))
        U.P[2][:] = 1 - U.P[1]
    else
        _maximum = mapreduce(maximum, maximum, EV)
        denominator = sum(exp.(U.U[ia]+U.β*EV[ia]-_maximum) for ia = 1:length(U.P))
        for ia = 1:length(U.P)-1
            U.P[ia][:] = exp.(U.U[ia]+U.β*(EV[ia])-_maximum)./denominator
        end
        U.P[end][:] = 1 - sum(U.P[ia] for ia = 1:length(U.P)-1)
    end
end

function βFP!(V, U, F)
    V.βFP .= sum(U.β*U.P[ia].*F[ia] for ia = 1:length(U.P))
end

function βFP!(EV::ExpectedValueFunction, U, F)
    for ia = 1:length(F)
        EV.βFP[ia] .= U.β*U.P[ia].*F[ia]
    end
end

function Ψ!(U, S, V::AbstractValueFunction, d)
    Γ!(U.U, S.F, U.β, V)
    Vᵏ⁺¹, Vᵏ = current(V), previous(V)
    V_diff = Vᵏ⁺¹-Vᵏ
    P!(U, S.F, Vᵏ + d*V_diff)
    βFP!(V, U, S.F)
    V.Vᵏ⁺¹ .= Vᵏ + (I-V.βFP)\V_diff
end

function Ψ!(U, S, EV::ExpectedValueFunction, d)
    Γ!(U.U, S.F, U.β, EV)
    EVᵏ⁺¹, EVᵏ = EV.EVᵏ⁺¹, EV.EVᵏ
    EV_diff = [EVᵏ⁺¹[ia]-EVᵏ[ia] for ia = 1:length(U.U)]
    P!(U, [EVᵏ[ia] + d*EV_diff[ia] for ia = 1:length(U.U)])
    βFP!(EV, U, S.F)
    ∂Γ∂EV = zeros(length(U.U)*size(S.F[1],1),length(U.U)*size(S.F[1],1))
    for ia = 1:length(U.U)
        for ia′ = 1:length(U.U)
            ∂Γ∂EV[(S.nX*(ia-1)+1):S.nX*(ia), (S.nX*(ia′-1)+1):S.nX*(ia′)] = S.F[ia].*U.P[ia′]'
        end
    end
    ev = vcat(EVᵏ...) + (I-U.β*∂Γ∂EV)\vcat(EV_diff...)
    for ia = 1:length(U.U)
        EV.EVᵏ⁺¹[ia] .= ev[(S.nX*(ia-1)+1):S.nX*ia]
    end
end
