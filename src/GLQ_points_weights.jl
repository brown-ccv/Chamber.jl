
include("./Data/GLQ_points_weights.jl")

function GLQ_points_weights_hard(n::Int64, GLQ_points_weights::Dict{Int64, Matrix{Float64}}=GLQ_points_weights)
    if haskey(GLQ_points_weights, n)
        A = GLQ_points_weights[n]
        quadpts = A[:,3]
        weights = A[:,2]
        return [quadpts, weights]
    else
        @warn("Weight and points not defined for this number.")
        return [[0], [0]]
    end
end