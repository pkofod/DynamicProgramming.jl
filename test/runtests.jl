using DynamicProgramming
using Plots
using Base.Test

# write your own tests here
problems = []
examples_path = "$(Pkg.dir("DynamicProgramming"))"*"/examples/"
for filename in ("7Days.jl", "AM2010.jl", "AM2016-Tauchen.jl", "Dixit.jl", "Hitsch.jl", "Rust.jl")
    push!(problems, include(examples_path*filename))
end

problem_N = []
vfi_time = []
newton_time = []
vfi_iter = []
newton_iter = []
for problem in problems
    U, S = problem()
    push!(problem_N, S.nX)
    W = IntegratedValueFunction(S)
    println("Problem is: ", problem)
    the_time = @elapsed iter = solve!(U, S, W)
    push!(vfi_time, the_time)
    push!(vfi_iter, iter)
    println("* VFIs: ", iter)
    println("* Elapsed: ", the_time)
    println("* Stopping criterion: ", DynamicProgramming.supnorm(W))
    println()
    W = IntegratedValueFunction(S)
    the_time = @elapsed iter = solve!(U, S, W, Newton())
    push!(newton_time, the_time)
    push!(newton_iter, iter)
    println("* Newton Iterations: ", iter)
    println("* Elapsed: ", the_time)
    println("* Stopping criterion: ", DynamicProgramming.supnorm(W))
    EV = ExpectedValueFunction(S)
    solve!(U, S, EV)
    EV = ExpectedValueFunction(S)
    solve!(U, S, EV, Newton())
    println()
end
