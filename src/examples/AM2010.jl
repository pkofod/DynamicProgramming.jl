function am(;N = 5, β = 0.95)
	# Aguirregabiria & Mira
	# 2 Model
	truepar = [-1.9;1.;2.] # fixed costs, ..., entry costs
	# State 1: Exogenous state
	nX1 = 2
	X1 = [0.;1.]
	F1 = [sparse([1. 0.; 1. 0.]), sparse([0. 1.; 0. 1.])]

	nX2 = N
	X2 = 1:nX2
	M = zeros(nX2, nX2)
	p = [0.2, 0.6, 0.2]
	M[1, 1:2] = [sum(p[1:2]), p[3]]
	M[end, end-1:end] = [p[1], sum(p[2:3])]
	M[2:end-1, :] = vcat([vcat(zeros(i-2), p, zeros(nX2-1-i))' for i = 2:nX2-1]...)

	S = States(State(X1, F1),
	           CommonState(X2, M))

	Z = [zeros(nX2*2, 3), # don't buy
	             [-ones(nX2) log(X2) -ones(nX2); # buy
				  -ones(nX2) log(X2)  zeros(nX2)]]

	U = LinearUtility(Z, β, truepar)

    return U, S
end
