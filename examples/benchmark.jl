using BenchmarkTools, ProgressMeter, MDPTools
cd(Pkg.dir("MDPTools")*"/examples/")
function add_model!(m, s)
    f = include(s)
    push!(m, f)
end

models = []
model_paths = ["AM2010.jl", "AM2016-Tauchen.jl", "Dixit.jl", "Hitsch.jl", "Rust.jl"]
for mp in model_paths
    add_model!(models, mp)
end

run(`clear`)
method_text = ["VFI   ", "Newton", "Poly  "]
@showprogress 1 "Benchmarking..." for m in models
    U, S = m()
    IV_VFI = @benchmark solve($U, $S, VFI())
    EV_VFI = @benchmark solve!($U, $S, ExpectedValueFunction($S.nX, length($U.U)), VFI())
    IV_Newton = @benchmark solve($U, $S, Newton())
    EV_Newton = @benchmark solve!($U, $S, ExpectedValueFunction($S.nX, length($U.U)), Newton())
    IV_Poly = @benchmark solve($U, $S, Poly())
    EV_Poly = @benchmark solve!($U, $S, ExpectedValueFunction($S.nX, length($U.U)), Poly())
#    IV_best = method_text[indmin((minimum(IV_VFI).time, minimum(IV_Newton).time, minimum(IV_Poly).time))]
#    EV_best = method_text[indmin((minimum(EV_VFI).time, minimum(EV_Newton).time, minimum(EV_Poly).time))]
    @printf "Dimension of state space %i\n" S.nX
    @printf "Time to run\n"
    @printf "* IntegratedValueFunction   ExpectedValueFunction\n"
    @printf "  VFI:           %f   VFI:         %f\n" minimum(IV_VFI).time/1e6 minimum(EV_VFI).time/1e6
    @printf "  Newton:        %f   Newton:      %f\n" minimum(IV_Newton).time/1e6 minimum(EV_Newton).time/1e6
    @printf "  Poly:          %f   Poly:        %f\n" minimum(IV_Poly).time/1e6 minimum(EV_Poly).time/1e6
#    @printf "  Best: %s      Best: %s\n" minimum(IV_best).time/1e6 minimum(EV_best).time/1e6
end
