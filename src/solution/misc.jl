current(V::AbstractValueFunction) = V.Vᵏ⁺¹
previous(V::AbstractValueFunction) = V.Vᵏ
current(EV::ExpectedValueFunction) = EV.EVᵏ⁺¹
previous(EV::ExpectedValueFunction) = EV.EVᵏ
