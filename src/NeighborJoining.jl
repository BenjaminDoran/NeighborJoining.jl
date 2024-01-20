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

merges(t::NJClust) = t.merges
heights(t::NJClust) = t.heights

export merges, heights, NJClust

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
