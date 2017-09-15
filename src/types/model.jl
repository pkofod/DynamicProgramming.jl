abstract type Horizon end
struct InifiniteHorizon
end
struct FiniteHorizon
    T
end
type DCModel{H<:Horizon}
    utilities
    states
    data
end
