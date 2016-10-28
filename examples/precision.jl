using Stocas
T = BigFloat
# State space
nX = 175
n = nX
o_maximum = 450
X = linspace(0, o_maximum, nX)
p = T[0.0937; 0.4475; 0.4459; 0.0127; 0.0002];
F1 = zeros(T, nX, nX)
offset = -1
for i = 1:nX, j = 1:length(p)
    i+offset+j > nX && continue
    F1[i,i+offset+j] = i+offset+j == nX ? sum(p[j:end]) : p[j]
end
F = [F1, F1[ones(Int64, nX), :]]
S = State(X, F)

Z1 = [zeros(T, nX) -0.001*T.(X)]
Z2 = [-ones(T, nX) zeros(T, nX)]
U = LinearUtility([Z1, Z2], T(0.9), T[11.;2.5])

V = IntegratedValueFunction(nX, 2, T)
EV = ExpectedValueFunction(nX, 2, T)
solve(U, S, V, Newton(1e-100, 30,0.))
solve(U, S, EV, Newton(1e-100, 30,0.))




EV = ExpectedValueFunction(nX, 2, T)
solve(U, S, EV, Newton(1e-100, 60,1.))
