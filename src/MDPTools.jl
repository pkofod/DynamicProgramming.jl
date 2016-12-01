# states need a shock type
# then a state or states without shock in combination
# with an integrated value function would add it.

module MDPTools
using StatsFuns, Distributions
        # Generic state(s)
export  State, States,
        # Discrete state(s)
        DiscreteState, DiscreteStates,
        # Special state(s)
        CommonState, EntryState,
        # Utility and policy related functions
        LinearUtility, Utility, policy,
        # Different value functions
        ValueFunction, IntegratedValueFunction, ExpectedValueFunction,
        # Solution algorithms
        VFI, Policy, Newton, Poly,
        # Solution relation functions
        bellman!, newton!, solve, solve!,
        # Simulation function and Data type for output
        simulate, Data,
        # Module with extras
        Extras

# Import here so we can extend them later
import Base: maximum, size
import Base.LinAlg: BlasFloat, A_mul_B!, Ac_mul_B!
import Base.BLAS.gemv!
# In-place gemv! using BLAS calls
A_mul_B!{T<:BlasFloat}(α::T, A::StridedMatrix{T}, x::StridedVector{T}, β::T, y::StridedVector{T}) = gemv!('N', α, A, x, β, y)
Ac_mul_B!{T<:BlasFloat}(α::T, A::StridedMatrix{T}, x::StridedVector{T}, β::T, y::StridedVector{T}) = gemv!('C', α, A, x, β, y)
# Fallback for other types of matrices
A_mul_B!(α, A, x, β, y) = copy!(y, β*y+α*(A*x))

# types
include("types/utility.jl")
include("types/states.jl")
include("types/valuefunctions.jl")

# solution
include("solution/misc.jl")
include("solution/bellman.jl")
include("solution/newton.jl")
include("solution/solve.jl")

# data

# simulation
include("simulate/data.jl")
include("simulate/simulate.jl")

# misc
include("Extras.jl")
end # module
