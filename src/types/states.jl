abstract AbstractState

"""
Represents a discrete, finite state variable.
"""
type DiscreteState{Ta <: AbstractMatrix} <: AbstractState
    X::AbstractVector
    F::Vector{Ta}
    nX::Int64
end
DiscreteState(X::AbstractVector, F::AbstractVector) = DiscreteState(X, F, length(X))

# Helpers
State(X::AbstractVector, F::AbstractVector) = DiscreteState(X, F)
EntryState(;dense=true) = State([0.;1.],
                                  dense ? [sparse([1. 0.; 1. 0.]), sparse([0. 1.; 0. 1.])]:
                                           [[1. 0.; 1. 0.], [0. 1.; 0. 1.]])

# Exogenous states that are common to all agents
type CommonState{Ta <: AbstractMatrix} <: AbstractState
    X::AbstractVector
    F::Vector{Ta}
    nX::Int64
end
CommonState(X, F) = CommonState(X, [F,], length(X)) # should be fixed!

# Experimental
immutable Terminator
    choices::AbstractVector
end
immutable TerminatorMatrix <: AbstractArray
end

# The multivariate version of the individual states
type DiscreteStates{Tm <: AbstractMatrix, Tv<:AbstractVector} <: AbstractState
    val_to_ind::Dict
    ind_to_val::Dict
    X::Vector{Tv}
    Fs::Vector{Union{Vector{Tm}}}
    F::Union{Vector{SparseMatrixCSC{Float64,Int64}}, Vector{Matrix{Float64}}}
    nX::Int64
    tup::Tuple
    market::Vector{Int64} # index of CommonStates #FIXME change to common
end

function DiscreteStates(terminators::Union{Void, Terminator}, states::Union{DiscreteState, CommonState}...)
    nA = maximum(ifelse(typeof(s) <: AbstractState, length(s.F), 1) for s in states)
    F = [ones(1,1) for a = 1:nA]
    Fs = Vector{AbstractMatrix}[]
    lens = []
    n = 1
    ind_market = Int64[]
    for (is, state) in enumerate(states)
        if typeof(state) <: CommonState
            push!(ind_market, is)
        end
        n *= length(state.X)
        for ia = 1:nA
            F[ia] = kron(F[ia], get_F(state,ia))
        end
        push!(Fs, state.F)
        push!(lens, length(state.X))
    end
    tup = (reverse(lens)...)
    val_to_ind = Dict(reverse(ind2sub(tup, i)) => i for i = 1:n)
    ind_to_val = Dict(i => reverse(ind2sub(tup, i)) for i = 1:n)
    DiscreteStates(val_to_ind, ind_to_val, Vector{Float64}[state.X for state in states], Fs, F, n, tup, ind_market)
end

get_F(state::AbstractState, ia) = state.F[ia]
get_F(state::CommonState, ia) = state.F[1]

# Helper function
States(states::Union{DiscreteState, CommonState}...) = DiscreteStates(nothing, states...)

#

# No-op States "constructor"
States(S::DiscreteStates) = S
