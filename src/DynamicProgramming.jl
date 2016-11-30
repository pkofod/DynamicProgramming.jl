# states need a shock type
# then a state or states without shock in combination
# with an integrated value function would add it.

# should also add a state which takes a function
# and some sort of "Parametric" type that
# specifies stuff

# should FP be in S?

# should maximum be maximum?

module DynamicProgramming
export State, States, CommonState,
        LinearUtility,
        ValueFunction,
        IntegratedValueFunction,
        ExpectedValueFunction,
        VFI, Policy, Newton, Poly,
        Γ!, Ψ!, solve

        import Base: maximum
        import Base.LinAlg: BlasFloat, A_mul_B!, Ac_mul_B!
        import Base.BLAS.gemv!
        A_mul_B!{T<:BlasFloat}(α::T, A::StridedMatrix{T}, x::StridedVector{T}, β::T, y::StridedVector{T}) = gemv!('N', α, A, x, β, y)
        Ac_mul_B!{T<:BlasFloat}(α::T, A::StridedMatrix{T}, x::StridedVector{T}, β::T, y::StridedVector{T}) = gemv!('C', α, A, x, β, y)
        A_mul_B!(α, A, x, β, y) = copy!(y, β*y+α*(A*x))
# types
include("types/discrete.jl")
include("types/states.jl")
include("types/valuefunctions.jl")

# solution
include("solution/misc.jl")
include("solution/bellman.jl")
include("solution/newton.jl")
include("solution/solve.jl")
end # module
