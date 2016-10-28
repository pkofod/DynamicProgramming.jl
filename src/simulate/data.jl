immutable Data{T<:Integer}
    id::Vector{T} # Agent IDs
    a::Vector{T} # Choices
    x::Vector{T} # Observed flattened state
    xs::Matrix{T} # Observed states
    m::Vector{T} # Market ID
    #Y::PayoffVariable
    nobs::T # Number of observations
end
Data(id, A_obs::Vector, X_obs::Vector) = Data(id, A_obs, X_obs, ones(Int, length(A_obs)))
function Data(id, A_obs::Vector, X_obs::Vector, M)
    n_obs = length(A_obs)
    n_obs == length(X_obs) || throw(error("Number of observed states and decisions does not match."))
    # Construct Data type
    Data(id, A_obs, X_obs, fill(0, n_obs, 1), M, n_obs)
end

# Automatic flattening of states given matrix X of
# observed states in the same order they are used
# to construct the States type.
function Data(id, A_obs, X_obs::Matrix, X)
    length(A_obs) == size(X_obs, 1) || throw(error("Number of observed states and decisions does not match."))
    n_obs = length(A_obs)
    #run through X and construnct x according to Xva i guess
    maxsub =  ([maximum(X.X[i]) for i = 1:length(X.X)]...)
    x = Int64[sub2ind(maxsub, X_obs[i,:]...) for i = 1:n_obs]
    Data(id, A_obs, X_obs, ones(n_obs), n_obs)
end

function Base.display(data::Data)
    @printf "Data instance with\n"
    @printf "* %i agents\n" length(unique(data.id))
    @printf "* %i observations\n" data.nobs
    @printf "* %i group(s)\n" length(unique(data.m))
end
