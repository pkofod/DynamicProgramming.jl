abstract AbstractState
immutable Terminator
    choices::AbstractVector
end
immutable TerminatorMatrix <: AbstractArray
end
"""
    State(Xval, F)

Constructs a State type instance with the fields

- X
- Xval
- F
- nX

"""

type State{Ta <: AbstractMatrix} <: AbstractState
    X::Vector{Int64}
    Xval::AbstractVector
    F::Vector{Ta}
    nX::Int64
end

State(Xval, F) = State([1:length(Xval);], Xval, F, length(Xval))


type CommonState{Ta <: AbstractMatrix} <: AbstractState
    X::Vector{Int64}
    Xval::AbstractVector
    F::Vector{Ta}
    nX::Int64
end

CommonState(Xval, F) = CommonState([1:length(Xval);], Xval, [F,], length(Xval)) # should be fixed!


type States{Tm <: AbstractMatrix, Tv<:AbstractVector} <: AbstractState
    X::Vector{Int64}
    val_to_ind::Dict
    ind_to_val::Dict
    Xval::Vector{Tv}
    Fs::Vector{Union{Vector{Tm}}}
    F::Union{Vector{SparseMatrixCSC{Float64,Int64}}, Vector{Matrix{Float64}}}
    nX::Int64
    tup::Tuple
    market::Vector{Int64}
end
#States(S::State) = States(S.X, Vector[S.Xval], [S.F], S.F, S.nX, (S.nX,))

"""
helper function for States generator
"""
get_F(state::AbstractState, ia) = state.F[ia]
get_F(state::CommonState, ia) = state.F[1]

States(states::Union{State, CommonState}...) = States(nothing, states...)
function States(terminators::Union{Void, Terminator}, states::Union{State, CommonState}...)
    nA = maximum(ifelse(typeof(s) <: State, length(s.F), 1) for s in states)
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
        push!(lens, length(state.Xval))
    end
    tup = (reverse(lens)...)
    val_to_ind = Dict(reverse(ind2sub(tup, i)) => i for i = 1:n)
    ind_to_val = Dict(i => reverse(ind2sub(tup, i)) for i = 1:n)
    States([1:n;], val_to_ind, ind_to_val, Vector{Float64}[state.Xval for state in states], Fs, F, n, tup, ind_market)
end


States(S::States) = S
