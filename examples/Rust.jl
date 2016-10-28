using Stocas
function rust(;N = 175, β = 0.95, sparse = false)
    # State space
    nX = N
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

    # We can handle sparse matrices
    F = [F1, F1[ones(Int64, nX), :]]
    sparse == true && map!(sparse, F)
    S = State(:mileage, X, F)

    Z1 = [zeros(nX) -0.001*X]
    Z2 = [-ones(nX) zeros(nX)]

    U = LinearUtility(["replace", "repair"],[Z1, Z2], β, [11.;2.5])

    return U, S
end
U, S = rust(N=1500)
V, iter = solve!(U, S)
D = simulate(U, S, 10, 1000, 10)
v = U.U + V.βEV
ṽ = v[2]-v[1]
Fe = locallinear(ṽ[D.x], D.a-1)
using UnicodePlots
scatterplot(D.x, Fe)
Dx_permutation = sortperm(D.x)
Fe_sorted = Fe[Dx_permutation]
lineplot(D.x[Dx_permutation][2:end], (D.x[Dx_permutation][2:end]-D.x[Dx_permutation][1:end-1]).*(Fe_sorted[2:end]-Fe_sorted[1:end-1]))
