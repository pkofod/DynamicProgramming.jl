using Stocas, ProgressMeter, Plots
cd(Pkg.dir("Stocas")*"/src/examples/")
function add_model!(m, s)
    f = include(s)
    push!(m, f)
end

models = []
model_paths = ["AM2010.jl", "AM2016-Tauchen.jl", "Dixit.jl", "Hitsch.jl", "Rust.jl"]
for mp in model_paths
    add_model!(models, mp)
end

N_range = 3:50:500
β_range = [0.3, 0.8, (0.89:0.10:0.99)...]
X1 = []
X2 = []
Y = []

for model_number = [1,3,5]
    push!(X1, zeros(length(N_range), length(β_range)))
    push!(X2, zeros(length(N_range), length(β_range)))
    push!(Y, zeros(length(N_range), length(β_range)))
    @showprogress 1 "Solving $(models[model_number])..." for (i, N) = enumerate(N_range)
        if model_number == 5
            N *= 2
        end
        for (j, β) = enumerate(β_range)
            if N > 0
                Y[end][i,j] = @elapsed minimum([solve(models[model_number](N=N, β = β)...,Newton()) for i = 1:30])
                X1[end][i,j] = models[model_number](N=N, β = β)[2].nX
            else
                Y[end][i,j] = @elapsed minimum([solve(models[model_number]()...,Newton()) for i = 1:30])
                X1[end][i,j] = models[model_number](N=N, β = β)[2].nX
            end
            X2[end][i,j] =  models[model_number](N=N, β = β)[1].β
        end
    end
end

plot(_x1, _y, c=:black)
plot!(__x1, __y, c=:red)
plot!(x1, y, c=:blue)
p1=plot(x1, y, labels = hcat(β_range...))
p2=plot(x2, y, labels = hcat(N_range...))
p=plot(p1,p2)
coef =
coef = [ones(length(x1)) x1 x1.^2 x2]\Float64.(y)

predict(x1, x2) = dot([1 x1 x1.^2 x2], coef)
predict(500,0.95)

# We basically see that AM and dixit have VERY simlar solution times as function of
#   x1 and x2!
