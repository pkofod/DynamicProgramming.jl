parsedims(dims, S) = (dims[i] == Colon() ? Colon() : findfirst(S.X[i], dims[i]) for i = length(dims):-1:1)

@recipe function f(dims, S::Stocas.DiscreteStates, V::Stocas.AbstractValueFunction)
        xind = findfirst(dims, Colon())
        X = S.X[xind]
        Vv = reshape(Stocas.current(V), S.tup...)[parsedims(dims, S)...]
        linetype --> :line
        xlabel --> S.names[xind]
        ylabel --> "value"
        X, Vv
end

@recipe function f(S::Stocas.DiscreteState, V::Stocas.AbstractValueFunction)
        X = S.X
        Vv = Stocas.current(V)
        linetype --> :line
        xlabel --> S.names[xind]
        ylabel --> "value"
        X, Vv
end

@recipe function f(dims, S::Stocas.DiscreteStates, V)
        xind = findfirst(dims, Colon())
        X = S.X[xind]
        Vv = reshape(V, S.tup...)[parsedims(dims, S)...]
        linetype --> :line
        xlabel --> S.names[xind]
        X, Vv
end

@recipe function f{T<:Real}(dims, S::Stocas.DiscreteStates, V::Vector{Vector{T}})
        xind = findfirst(dims, Colon())
        X = S.X[xind]
        Y = []
        for i = 1:length(V)
            push!(Y, reshape(V[i], S.tup...)[parsedims(dims, S)...])
        end
        linetype --> :line
        xlabel --> S.names[xind]
        X, Y
end

@recipe function f(dims, S::Stocas.DiscreteStates, U::LinearUtility, y = :utility)
        n = length(U.U)
        if y == :utility
            obj = U.U
        elseif y == :policy
            obj = policy(U)
        else
            error("You need to choose :utility or :policy.")
        end
        xind = findfirst(dims, Colon())
        X = S.X[xind]
        Y = []
        for i = 1:n
            push!(Y, reshape(obj[i], S.tup...)[parsedims(dims, S)...])
        end
        linetype --> :line
        xlabel --> S.names[xind]
        ylab = y == :utility ? "utility" : "probability"
        ylabel --> ylab
        labels --> hcat(U.names...)
        X, Y
end

@recipe function f(S::Stocas.DiscreteState, U::LinearUtility, y = :utility)
        n = length(U.U)
        if y == :utility
            obj = U.U
        elseif y == :policy
            obj = policy(U)
        else
            error("You need to choose :utility or :policy.")
        end
        X = S.X
        Vv = []
        for i = 1:n
            push!(Vv, obj[i])
        end
        linetype --> :line
        xlabel --> S.name
        ylab = y == :utility ? "utility" : "probability"
        ylabel --> ylab
        labels --> hcat(U.names...)
        X, Vv
end

@recipe function f(S::Stocas.DiscreteState, V)
        X = S.X
        linetype --> :line
        xlabel --> S.name
        X, V
end
