abstract type AbstractState end

"""
Represents a discrete, finite state variable.
"""
type DiscreteState{Ta <: AbstractMatrix, Tn<:Any} <: AbstractState
    name::Tn
    X::AbstractVector
    F::Vector{Ta}
    nX::Int64
end
DiscreteState(X::AbstractVector, F::AbstractVector) = DiscreteState("state", X, F, length(X))
DiscreteState(name, X::AbstractVector, F::AbstractVector) = DiscreteState(name, X, F, length(X))

# Helpers
State(name, X::AbstractVector, F::AbstractVector) = DiscreteState(name, X, F, length(X))
State(X::AbstractVector, F::AbstractVector) = DiscreteState(X, F)
EntryState(;dense=true) = State("incumbancy", ["don't enter", "enter"],
                                  dense ? [sparse([1. 0.; 1. 0.]), sparse([0. 1.; 0. 1.])]:
                                           [[1. 0.; 1. 0.], [0. 1.; 0. 1.]])

# Exogenous states that are common to all agents
type CommonState{Ta <: AbstractMatrix, Tn<:AbstractString} <: AbstractState
    name::Tn
    X::AbstractVector
    F::Vector{Ta}
    nX::Int64
end
CommonState(X, F) = CommonState("common state", X, [F,], length(X))
CommonState(name, X, F) = CommonState(name, X, [F,], length(X))

# The multivariate version of the individual states
type DiscreteStates{Tm <: AbstractMatrix} <: AbstractState
    names::Vector
    val_to_ind::Dict
    ind_to_val::Dict
    X::Vector
    Fs::Vector{Union{Vector{Tm}}}
    F::Union{Vector{SparseMatrixCSC{Float64,Int64}}, Vector{Matrix{Float64}}}
    nX::Int64
    tup::Tuple
    market::Vector{Int64} # index of CommonStates #FIXME change to common
end

function DiscreteStates(states::Union{DiscreteState, CommonState}...)
    nA = maximum(ifelse(typeof(s) <: AbstractState, length(s.F), 1) for s in states)
    names = Symbol[]
    for (i, s) in enumerate(states)
        if typeof(s.name) <: String
            push!(names, s.name)
        else
            push!(names, Symbol("state$i"))
        end
    end
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
    DiscreteStates(names, val_to_ind, ind_to_val, Vector[state.X for state in states], Fs, F, n, tup, ind_market)
end

get_F(state::AbstractState, ia) = state.F[ia]
get_F(state::CommonState, ia) = state.F[1]

States(S::Union{DiscreteState, CommonState}...) = DiscreteStates(S...)
# No-op States "constructor"
States(S::DiscreteStates) = S

names(S::AbstractState) = S.names
names(S::Union{CommonState, DiscreteState}) = S.name


function Base.display(S::DiscreteStates)
    @printf "Discrete states\n"
    @printf " * Number of state variables: %s\n" length(S.tup)
    @printf " * Total number of states:    %s\n" S.nX
    @printf " * Individual states:           \n"
    for ix = 1:length(S.X)
        @printf "  %d) %s (n: %s)\n" ix S.names[ix] length(S.X[ix])
    end
end
