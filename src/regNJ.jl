###
# Naive neighbor joining
###
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
    n = size(d, 1)
    n == size(d,2) || ArgumentError("d must be a square matrix")
    merges = zeros(Int, n-1, 2)
    heights = zeros(Float64, n-1, 2)
    sd = zeros(Float64, 2n-1, 2n-1)
    sd[1:n, 1:n] .= d
    currentindices = zeros(Bool, n+n-1)
    currentindices[1:n] .= true
    regNJ!(merges, heights, sd, currentindices)
    return NJClust(merges, heights)
end


function regNJ!(merges::AbstractMatrix, heights::AbstractMatrix, d::AbstractMatrix{<:Number}, currentindices::AbstractVector{<:Bool})
    n=sum(currentindices)
    orig_idx = collect(1:(2n-1))

    for mergestep in 1:n-1
        tmpd = @view d[currentindices, currentindices]
        tmp_orig_idx = @view orig_idx[currentindices]
        idxs = _Q(tmpd, tmp_orig_idx)
        merges[mergestep, 1] = _mergeidx(idxs[1], n) 
        merges[mergestep, 2] = _mergeidx(idxs[2], n)
        heights[mergestep, 1], heights[mergestep, 2] = _distance_to_parent(d, currentindices, idxs[1], idxs[2])
        
        currentindices[idxs[1]] = false
        currentindices[idxs[2]] = false
        for c in findall(currentindices)
            d[c, n+mergestep] = _distance_to_new_node(d, idxs[1], idxs[2], c)
            d[n+mergestep, c] = d[c, n+mergestep]
        end
        currentindices[n+mergestep] = true
    end
end


function _Q(d::AbstractMatrix{<:Number}, orig_idx::AbstractVector)
    n = size(d, 1)
    min_Q = Inf
    min_i, min_j = 0, 0
    marginsums = vec(mapslices(sum, d, dims=1))
    for j in axes(d, 2), i in (j+1):n
        cur_Q = (n-2) * d[i,j] - marginsums[i] - marginsums[j]
        if cur_Q < min_Q
            min_Q = cur_Q
            min_i, min_j = orig_idx[i], orig_idx[j]
        end
    end
    return min_i, min_j, min_Q
end
