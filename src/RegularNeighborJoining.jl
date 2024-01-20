###
# Naive neighbor joining
###

module RegularNeighborJoining

using ArgCheck
using ..NeighborJoining: NJClust

"""
    regNJ(d::AbstractMatrix{<:Number})

regNJ algorithm is the traditional NeighborJoining algorithm from 

> Saitou, N. & Nei, M. The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4, 406-425 (1987).
 
This algorithm is guarenteed to infer the tree for additive distance matrices, but it does have an algorithmic complexity of `O(n^3)`, so it can be slow for distance matrices on the order of >10³.

args:
* d is an n by n square symetric distance matrix

returns:
* NJClust struct with fields merges and heights
"""
function regNJ(d::AbstractMatrix{<:Number})
    @argcheck allequal(size(d))
    
    n = size(d, 1)
    working_d = zeros(eltype(d), size(d))
    working_d .= d
    currentindices = trues(n)
    merges = zeros(Int, n-1, 2)
    heights = zeros(Float64, n-1, 2)
    regNJ!(merges, heights, working_d, currentindices)
    return NJClust(merges, heights)
end


function regNJ!(merges::AbstractMatrix, heights::AbstractMatrix, d::AbstractMatrix{<:Number}, currentindices::AbstractVector{<:Bool})
    n=sum(currentindices)
    nodelabels = collect(1:n)

    for mergestep in 1:n-1
        @debug "nodelabels = $(_mergeidx.(nodelabels, n))"
        @debug "d" d

        # find min dist to merge using Q metric
        min_i, min_j, min_Q = _Q(d)
        @debug min_i, min_j, min_Q

        # save merge to merge list
        merges[mergestep, 1] = _mergeidx(nodelabels[min_i], n) 
        merges[mergestep, 2] = _mergeidx(nodelabels[min_j], n)
        heights[mergestep, 1], heights[mergestep, 2] = _distance_to_parent(d, min_i, min_j)
        # update distance matrix with distances to new node
        k = min_j # put new node in place of node i
        for c in axes(d, 2)
            c in (min_i, min_j) && continue
            d[c, k] = d[k, c] = _distance_to_new_node(d, min_i, min_j, c)
        end
        # delete remaining node j
        nodelabels[k] = n+mergestep
        currentindices[min_i] = false # prep for deletion
        d = @view d[currentindices, currentindices]
        nodelabels = @view nodelabels[currentindices]
        deleteat!(currentindices, min_i)
    end
end

# """ return -x if less then or equal to n otherwise return x-n """
_mergeidx(i, n) = i ≤ n ? -i : i-n

function _Q(d::AbstractMatrix{<:Number})
    n = size(d,1)::Int64
    min_Q = Inf
    min_i, min_j = 0, 0
    marginsums = vec(sum(d, dims=1))
    for j in axes(d, 2), i in (j+1):lastindex(d, 1)
        @inbounds cur_Q = (n-2) * d[i,j] - marginsums[i] - marginsums[j]
        if cur_Q < min_Q
            min_i, min_j, min_Q = i, j, cur_Q
        end
    end
    return min_i, min_j, min_Q
end

# """ take distance matrix and current indices and calculates position of parent linking nodes i and j """
function _distance_to_parent(d::AbstractMatrix{<:Number}, i::Integer, j::Integer)
    n = size(d, 2)
    sum_i = sum(@view d[:,i])
    sum_j = sum(@view d[:,j])
    β = n > 2 ? (1/(2*(n-2)))*(sum_i - sum_j) : 0.0
    δa = ((0.5) * d[i, j]) + β
    δb = d[i, j] - δa
    # if branch length is negative 
    # exchange length to other branch
    # does not affect topology 
    # (Kuhner and Felsenstein, 1994)
    if δb < 0 
        δa -= δb
        δb = 0.
    elseif δa < 0
        δb -= δa
        δa = 0.
    end

    return δa,  δb
end

# """ distance from existing nodes to new nodes """
function _distance_to_new_node(d::AbstractMatrix{<:Number}, a::Integer, b::Integer, c::Integer)
    return (d[a,c] + d[b,c] - d[a, b]) / 2
end

end