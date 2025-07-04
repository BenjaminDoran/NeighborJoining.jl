module NeighborJoining
using Printf
using Base: hash, isequal, ==

"""
    NJClust(merges::Matrix{M}, heights::Matrix{H})

fields:
* merges is a n-1 x 2 matrix of integers: absolute values of negative integers indicate index into the distance matrix (i.e., leaves). positive integers are the index into the merge list (i.e., the kth internal node)
* heights is an n-1 x 2 matrix where each value is the distance from the left (1) or right (2) chield from its parent. Specifically `heights[i,j]` is the `j`th childs distance to the parent node, row `i`.
"""
struct NJClust{M, H}
    merges::Matrix{M}
    heights::Matrix{H}
end

"""
    merges(t::NJClust)
extracts the merge list from NJClust object. rowindex is the internal node id,
and each column represents the children nodes.
Leaf nodes are represented by negative integers corresponding to the index in the original distance matrix
"""
merges(t::NJClust) = t.merges

"""
    heights(t::NJClust)
extracts the heights list from NJClust object. rowindex is the internal node id,
and each column represents the distance from the internal node to it's left and right child respectively.
"""
heights(t::NJClust) = t.heights


"""
    order(clust::NJClust)

Returns the left-to-right order of leaf nodes (tips) in the phylogenetic tree.

Performs a depth-first traversal of the tree structure to determine the order in which
leaf nodes appear when reading the tree from left to right. This is useful for plotting
or arranging data according to the tree topology.

# Arguments
* `clust::NJClust`: A neighbor-joining clustering result containing the tree structure
  with merges and heights matrices

# Returns
* `Vector{Int}`: A vector of integers representing the indices of leaf nodes in their
  left-to-right order in the tree. The indices correspond to the original positions
  in the distance matrix used to construct the tree.

# Example
```jldoctest
julia> d = [
           0  5  9  9 8
           5  0 10 10 9
           9 10  0  8 7
           9 10  8  0 3
           8  9  7  3 0
       ];

julia> njclusts = regNJ(d)
NJClust{Int64, Float64}([-2 -1; -3 1; -4 2; -5 3], [3.0 2.0; 4.0 3.0; 2.0 2.0; 0.5 0.5])

julia> leaf_order = order(njclusts)
5-element Vector{Int64}:
 5
 4
 3
 2
 1
```
"""
function order(clust::NJClust)
	# Initialize
    leaf_order = Int[]

    # Stack to store nodes to process: (node_index, is_processed)
    stack = [(size(clust.merges, 1), false)]
    while !isempty(stack)
        node_idx, processed = pop!(stack)
        if node_idx < 0
            # Leaf node - add to result
            push!(leaf_order, abs(node_idx))
        elseif processed
            # Internal node already processed - skip
            continue
        else
            # Internal node - add children to stack (right first, then left)
            # This ensures left-to-right traversal
            push!(stack, (node_idx, true))  # Mark as processed
            push!(stack, (clust.merges[node_idx, 2], false))  # Right child
            push!(stack, (clust.merges[node_idx, 1], false))  # Left child
        end
    end

    return leaf_order
end

export merges, heights, order, NJClust

include("RegularNeighborJoining.jl")
using .RegularNeighborJoining: regNJ
export regNJ

include("FastNeighborJoining.jl")
using .FastNeighborJoining: fastNJ
export fastNJ

include("newickstring.jl")
export newickstring


Base.hash(a::NJClust, h::UInt) = hash(a.merges, hash(a.heights, hash(:NJClust, h)))
Base.isequal(a::NJClust, b::NJClust) = isequal(hash(a), hash(b))
Base.:(==)(a::NJClust, b::NJClust) = a.merges == b.merges && a.heights == b.heights

end
