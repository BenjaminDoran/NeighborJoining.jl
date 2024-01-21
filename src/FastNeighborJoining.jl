module FastNeighborJoining

###
# Fast Neighbor Joining
###
using ..NeighborJoining: NJClust
"""
    fastNJ(d::AbstractMatrix{<:Number})

fastNJ algorithm finds k independent pairs to merge and merges them for each iteration. 
This is significantly faster because it is not recalculating the full pairwise Q for each pair joined.

This algorithm is nearly additive, but there are instances where the topology 
'(a,(b,(c,d)))' is inferred as '((a,b), (c,d))'. This happens when 
'b' is not as close to 'c' or 'd' as it is to 'a', 
but 'b' is closer to the common ancestor '(c,d)' than it is to 'a'.


args:
* d is an n by n square symetric distance matrix

returns:
* NJClust struct with fields merges and heights

```@jldoctest
julia> d = [
           0  5  9  9 8
           5  0 10 10 9
           9 10  0  8 7
           9 10  8  0 3
           8  9  7  3 0
       ];

julia> njclusts = fastNJ(d)
NJClust{Int64, Float64}([-2 -1; -5 -4; 1 -3; 3 2], [3.0 2.0; 1.0 2.0; 3.0 4.0; 1.0 1.0])

julia> nwstring = newickstring(njclusts)
"(5:5.000000e-01,(4:2.000000e+00,(3:4.000000e+00,(2:3.000000e+00,1:2.000000e+00):3.000000e+00):2.000000e+00):5.000000e-01):0.000000e+00;"
```
"""
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
    wQk = @view Qk[firstindex(Qk):binomial(n,2)]
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

# """ return -x if less then or equal to n otherwise return x-n """
_mergeidx(i, n) = i ≤ n ? -i : i-n

# """ take distance matrix and current indices and calculates position of parent linking nodes i and j """
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

# """ distance from existing nodes to new nodes """
function _distance_to_new_node(d::AbstractMatrix{<:Number}, a::Integer, b::Integer, c::Integer)
    return (d[a,c] + d[b,c] - d[a, b]) / 2
end

end