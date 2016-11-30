using DynamicProgramming, Plots
pgfplots()
srand(124)

# State space
nX = 175
n = nX
o_maximum = 450
X = linspace(0, o_maximum, nX)
p = [0.0937; 0.4475; 0.4459; 0.0127; 0.0002];
F1 = zeros(nX, nX)
offset = -1
for i = 1:nX, j = 1:length(p)
    i+offset+j > nX && continue
    F1[i,i+offset+j] = i+offset+j == nX ? sum(p[j:end]) : p[j]
end
F = [F1, F1[ones(Int64, nX), :]]
S = State(X, F)

Z1 = [zeros(nX) -0.001*X]
Z2 = [-ones(nX) zeros(nX)]
U = LinearUtility([Z1, Z2], 0.9999, [11.;2.5])

v = ValueFunction(nX, 2)
V = IntegratedValueFunction(nX, 2)
solve(U, S, V, Poly())

EV = ExpectedValueFunction(nX, 2)
solve(U, S, EV, Poly())

plot(X, V.βEV./U.β, labels = ["EV(1)" "EV(2)"], ylims= (-3664,-3650), xlabel = "mileage", ylabel= "EV")
savefig("/home/pkm/Dropbox/phd/Projekter/Structural Econometrics/right.pdf")

F_wrong = [F1, [ones(Int64, nX) zeros(nX, nX-1)]]
S_wrong = State(X, F_wrong)
V_wrong = IntegratedValueFunction(nX, 2)
solve(U, S_wrong, V_wrong, Poly())
plot(X, V_wrong.βEV./U.β, labels = ["EV(1)" "EV(2)"], ylims= (-3576,-3562), xlabel = "mileage", ylabel= "EV")
savefig("/home/pkm/Dropbox/phd/Projekter/Structural Econometrics/wrong.pdf")
