module NeighborJoining
using Printf
using Base: hash, isequal, ==


struct NJClust{M, H}
    merges::Matrix{M}
    heights::Matrix{H}
end

merges(t::NJClust) = t.merges
heights(t::NJClust) = t.heights

Base.hash(a::NJClust, h::UInt) = hash(a.merges, hash(a.heights, hash(:NJClust, h)))
Base.isequal(a::NJClust, b::NJClust) = isequal(hash(a), hash(b))
Base.:(==)(a::NJClust, b::NJClust) = a.merges == b.merges && a.heights == b.heights

export fastNJ, regNJ, newickstring, merges, heights, NJClust

###
# Naive neighbor joining
###

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

###
# Fast Neighbor Joining
###

function fastNJ(d::AbstractMatrix{<:Number}) 
    n = size(d, 1)
    n == size(d,2) || ArgumentError("d must be a square matrix")
    merges = zeros(Int, n-1, 2)
    heights = zeros(Float64, n-1, 2)
    sd = zeros(eltype(d), 2n-1, 2n-1)
    sd[1:n, 1:n] .= d
    currentindices = zeros(Bool, n+n-1)
    currentindices[1:n] .= true
    Qk = Vector{Tuple{Int, Int, Float64}}(undef, binomial(n, 2))
    fastNJ!(merges, heights, sd, currentindices, Qk)
    return NJClust(merges, heights)
end

function fastNJ!(merges::AbstractMatrix, heights::AbstractMatrix, d::AbstractMatrix{<:Number}, currentindices::AbstractVector{<:Bool}, Qk::AbstractVector{Tuple{Int64, Int64, Float64}})
    n = sum(currentindices)
    mergestep = 1
    marginsums = zeros(size(d, 1))
    originalindices = collect(1:(2n-1))
    while (sum(currentindices) > 1) && (mergestep < n)
        wd = @view d[currentindices, currentindices]
        wmarginsums = @view marginsums[currentindices]
        woriginalindices = @view originalindices[currentindices]
        newjoins = _Qlist!(Qk, wd, currentindices, wmarginsums, woriginalindices)
        for (i, j, qij) in newjoins
            @debug "i=$i, j=$j, q=$qij"
            # add nodes to merge list
            merges[mergestep, 1] = _mergeidx(i, n) 
            merges[mergestep, 2] = _mergeidx(j, n)
            # collect each nodes distance to new node
            heights[mergestep, 1], heights[mergestep, 2] = _distance_to_parent(d, currentindices, i, j)
            
            # remove nodes from additional merges
            currentindices[i] = false
            currentindices[j] = false

            # find all remaining nodes distance to new node
            for c in findall(currentindices)
                d[c, n+mergestep] = _distance_to_new_node(d, i, j, c)
                d[n+mergestep, c] = d[c, n+mergestep]
            end

            # add new node
            currentindices[n+mergestep] = true
            # increment how many merges have been performed 
            mergestep += 1
            # print("\r$(lpad(mergestep, 10)) / $n")
        end
    end
end

function _Qlist!(
    Qk::AbstractVector{Tuple{Int64, Int64, Float64}},
    d::AbstractMatrix{<:Number},
    currentindices::AbstractVector{<:Bool},
    marginsums::AbstractVector{<:Number},
    indxs::AbstractVector{<:Number}
    )
    n = sum(currentindices)
    # indxs = findall(currentindices)
    marginsums .= vec(mapslices(sum, d, dims=1))
    binomial(n,2) <= length(Qk) || 
        ArgumentError("Qk is too short!")
    wQk = @view Qk[begin:binomial(n,2)]
    k = 1
    for j in axes(d, 2), i in (j+1):n
        wQk[k] = (indxs[i], indxs[j], (n-2) * d[i,j] - marginsums[i] - marginsums[j])
        k += 1
    end
    # sort by adjusted distance between nodes (Q value)
    sort!(wQk, by=x->x[3])
    @debug display(wQk)
    return _get_independent_merges(wQk)
end

function _get_independent_merges(q)
    s = Set{Int64}()
    isnew(x) = length(s) != length(push!(s, x))
    return Iterators.filter((x)->isnew(x[1]) & isnew(x[2]), q)
end

# function _get_independent_merges(q)
#     s = Set{Int64}()
#     isnew(x) = length(s) != length(push!(s, x))
#     return Iterators.takewhile((x)->isnew(x[1]) & isnew(x[2]), q)
# end

# function _get_independent_merges(q)
#     s = Set{Int64}()
#     isnew(x) = length(s) != length(push!(s, x))
#     return (
#         (i, j) for (i,j,_) in q if 
#             isnew(i) && isnew(j)
#     )
# end


""" return -x if less then or equal to n otherwise return x-n """
_mergeidx(i, n) = i ≤ n ? -i : i-n

""" take distance matrix and current indices and calculates position of parent linking nodes i and j """
function _distance_to_parent(d::AbstractMatrix{<:Number}, currentindices::AbstractVector{<:Bool}, i::Integer, j::Integer)
    n = sum(currentindices)
    sum_i = sum(@view d[currentindices,i])
    sum_j = sum(@view d[currentindices,j])
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

""" distance from existing nodes to new nodes """
function _distance_to_new_node(d::AbstractMatrix{<:Number}, a::Integer, b::Integer, c::Integer)
    return (d[a,c] + d[b,c] - d[a, b]) / 2
end



###
# Convert to newick string
###
function newickstring(njc::NJClust, tiplabels=1:(size(heights(njc), 1)+1); labelinternalnodes=false) 
    newickstring(merges(njc), heights(njc), string.(tiplabels); labelinternalnodes)
end
function newickstring(merges::A, heights::B, tiplabels::AbstractVector{<:String}; labelinternalnodes=false) where {
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

# """    
#     patristic_distances(t::Node)
# shortest branch length path between all leafs with `t` as ancester
# sibling nodes `i`, and `j` of parent `p` would have patristic `distance(i, p) + distance(j, p)`
# """
# function patristic_distances(tree::Node)
#     leaves = getleaves(tree) |> x->sort(x; by=name)
#     dists = zeros(length(leaves), length(leaves))
#     for j in axes(dists, 2), i in (j+1):lastindex(dists, 1)
#         dists[i, j] = NewickTree.getdistance(leaves[i], leaves[j])
#     end
#     return Symmetric(dists, :L)
# end

end
