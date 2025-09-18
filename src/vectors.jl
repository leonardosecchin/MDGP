# Z = X
@inline function cp!(
            idxs,
            Z::AbstractVecOrMat{Float64},
            X::AbstractVecOrMat{Float64}
            )
    @inbounds for i in idxs
        Z[i] = X[i]
    end
end

# Z = a*X
@inline function aX!(
            idxs,
            Z::AbstractVecOrMat{Float64},
            a::Float64,
            X::AbstractVecOrMat{Float64}
            )
    @inbounds for i in idxs
        Z[i] = a*X[i]
    end
end

# Z = X + Y
@inline function XpY!(
            idxs,#::Vector{Int64},
            Z::AbstractVecOrMat{Float64},
            X::AbstractVecOrMat{Float64},
            Y::AbstractVecOrMat{Float64}
            )
    @inbounds for i in idxs
        Z[i] = X[i] + Y[i]
    end
end

# Z = X + a*Y
@inline function XpaY!(
            idxs,#::Vector{Int64},
            Z::AbstractVecOrMat{Float64},
            X::AbstractVecOrMat{Float64},
            a::Float64,
            Y::AbstractVecOrMat{Float64}
            )
    @inbounds for i in idxs
        Z[i] = X[i] + a*Y[i]
    end
end

# Z = X + a*Y
@inline function XdotY(
            idxs,#::Vector{Int64},
            X::AbstractVecOrMat{Float64},
            Y::AbstractVecOrMat{Float64}
            )
    d = 0.0
    @inbounds for i in idxs
        d += X[i] * Y[i]
    end
    return d
end

# convert v into a UnitRange if possible
function consec_range(v)
    if isempty(v)
        return v
    else
        sort!(v)
        if v[end] - v[1] + 1 == length(v)
            return v[1]:v[end]
        else
            return v
        end
    end
end
