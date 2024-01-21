
###
# Convert to newick string
###
"""
    newickstring(njc::NJClust, tiplabels=AbstractVector{<:String}; labelinternalnodes=false)
    newickstring(merges::AbstractArray{<:Integer}, heights::AbstractArray{<:AbstractFloat}, tiplabels::AbstractVector{<:String}; labelinternalnodes=false)

Converts a list of merges into a newicktree formatted string.

# args:
* njc: is a struct that has merges and heights
* merges and heights from a NJClust struct:
    * merges is a n-1 x 2 matrix of integers: absolute values of negative integers indicate index into the distance matrix (i.e., leaves). positive integers are the index into the merge list (i.e., the kth internal node)
    * heights is an n-1 x 2 matrix where each value is the distance from the left (1) or right (2) chield from its parent. Specifically `heights[i,j]` is the `j`th childs distance to the parent node, row `i`.
* tiplabels: vector of string labels corresponding to the order of leaves in the distance matrix
* labelinternalnodes: whether to generate node labels for the internal nodes. defaults to `false`.

# returns:
* newicktree formatted string

# example:

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

julia> nwstring = newickstring(njclusts)
"(5:5.000000e-01,(4:2.000000e+00,(3:4.000000e+00,(2:3.000000e+00,1:2.000000e+00):3.000000e+00):2.000000e+00):5.000000e-01):0.000000e+00;"
```
"""
function newickstring(njc::NJClust, tiplabels=1:(size(heights(njc), 1)+1); labelinternalnodes=false) 
    newickstring(merges(njc), heights(njc), string.(tiplabels); labelinternalnodes)
end
function newickstring(merges::A, heights::B, tiplabels::AbstractVector{<:AbstractString}; labelinternalnodes=false) where {
    A<:AbstractArray{<:Integer}, B<:AbstractArray{<:AbstractFloat}
}
    r = size(heights, 1)
    dist = 0.0
    _newickstring(view(merges, :, :), view(heights, :, :), r, r, dist, view(tiplabels, :); labelinternalnodes) * ";"
end
function _newickstring(merges::A, heights::B, i::C, p::C, dist::AbstractFloat, tiplabels::D; labelinternalnodes=false)::String where {
    A<:AbstractArray{<:Integer}, B<:AbstractArray{<:AbstractFloat},
    C<:Integer, D<:AbstractVector{<:AbstractString}
}
    j::Int64 = merges[i,1] # left subtree pointer
    k::Int64 = merges[i,2] # right subtree pointer
    a::String = if j < 0 # if tip format tip
            tiplabels[abs(j)] * ':' * @sprintf("%e", heights[i, 1])
        else # recurse and format internal node
            _newickstring(view(merges, :, :), view(heights, :, :), j, i, heights[i, 1], view(tiplabels, :); labelinternalnodes)
        end
    b::String = if k < 0 # if tip format tip
            tiplabels[abs(k)] * ':' * @sprintf("%e", heights[i, 2])
        else # recurse and format internal node
            _newickstring(view(merges, :, :), view(heights, :, :), k, i, heights[i, 2], view(tiplabels, :); labelinternalnodes)
        end
    nid = labelinternalnodes ? "node" * string(size(heights, 1) + i + 1) : ""
    stringdist = @sprintf("%e", dist)
    _newick_merge_strings(a,b,nid,stringdist)
end
function _newick_merge_strings(a::S, b::S, n::S, d::S) where S<:String
    '(' *  a *  ',' * b * ')' * n * ':' * d
end