module DynamicNeighborJoining

using ArgCheck
using Base.Threads: @threads, nthreads
using ..NeighborJoining: NJClust

"""
    dynamicNJ(d::AbstractMatrix{<:Number}; parallel::Bool=true)

`dynamicNJ` is the exact neighbor-joining algorithm accelerated with the dynamic
algorithm from

> Clausen, "Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining", Bioinformatics (2023). <https://doi.org/10.1093/bioinformatics/btac774>

Like `regNJ` it is an *exact* implementation of the canonical Saitou & Nei
neighbor-joining criterion, so it recovers the same tree, but it is substantially
faster on large inputs.

args:
* `d` is an n by n square symmetric distance matrix

keyword args:
* `parallel` toggles multithreaded search/update (default `true`)

returns:
* `NJClust` struct with fields `merges` and `heights`

examples:

```jldoctest
julia> d = [
           0  5  9  9 8
           5  0 10 10 9
           9 10  0  8 7
           9 10  8  0 3
           8  9  7  3 0
       ];

julia> njclusts = dynamicNJ(d)
NJClust{Int64, Float64}([-1 -2; 1 -3; 2 -4; 3 -5], [2.0 3.0; 3.0 4.0; 2.0 2.0; 0.5 0.5])

julia> nwstring = newickstring(njclusts)
"((((1:2.000000e+00,2:3.000000e+00):3.000000e+00,3:4.000000e+00):2.000000e+00,4:2.000000e+00):5.000000e-01,5:5.000000e-01):0.000000e+00;"
```

## Extended help

The speed-up comes from the observation that the neighbor-joining criterion `Q`
is monotonically weakened with each iteration:
the minimum value of `Q` found within a row in one iteration is a valid lower
bound for that row in every subsequent iteration. Rows whose lower bound already
exceeds the current global optimum can therefore be skipped entirely, so most of
the distance matrix does not need to be re-scanned each iteration.

The initial full scan and the per-iteration row search are parallelized with
Julia's `Base.Threads`; start Julia with multiple threads (e.g. `julia -t auto`)
to take advantage of this. Threading is engaged adaptively only for iterations
whose estimated work is large enough to amortise the thread-spawn cost, so
well-structured inputs (which skip most rows) never run slower than serial. The
result is bit-for-bit deterministic regardless of the number of threads or the
value of `parallel`.
"""
function dynamicNJ(d::AbstractMatrix{<:Number}; parallel::Bool = true)
    @argcheck allequal(size(d))

    n = size(d, 1)
    merges = zeros(Int, max(n - 1, 0), 2)
    heights = zeros(Float64, max(n - 1, 0), 2)
    n <= 1 && return NJClust(merges, heights)

    # Working type: promote integer matrices to Float64 so branch lengths match regNJ.
    Tf = float(eltype(d))

    # Full symmetric working distance matrix. We keep both triangles because the
    # update step reads distances to the merged nodes in arbitrary row order.
    D = Matrix{Tf}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        D[i, j] = d[i, j]
    end

    # Divergence: R[i] is the sum of distances from node i to every other active node.
    R = vec(sum(D, dims = 1))

    # Bookkeeping for the fixed-slot scheme (rows/columns are never compacted).
    obsolete = falses(n)          # slots retired after their node was merged
    idmap = collect(1:n)          # slot -> node id (leaves 1:n, internal nodes n+k)
    Q = fill(typemax(Tf), n)      # per-row minimum join criterion (lower bounds)

    # Only spawn threads once there is enough work to outweigh the thread-spawn
    # barrier. The dynamic algorithm skips most rows in later iterations, so blindly
    # threading every iteration is slower than serial. The serial and threaded paths
    # are bit-identical, so choosing per-iteration is safe.
    threads_ok = parallel && nthreads() > 1

    # --- First iteration: full scan of the matrix. ---
    # The initial scan is a single dense O(n^2) pass and parallelises cleanly, so it
    # uses a lower threshold than the per-iteration search.
    coeff = Tf(n - 2)
    x, y, d_xy = _full_scan!(Q, D, R, n, coeff, threads_ok && n >= _SCAN_PAR_MIN)
    _record_merge!(merges, heights, 1, x, y, d_xy, R, idmap, n, n)

    prev_z = 0
    if n > 2
        prev_z = _update!(D, R, obsolete, idmap, x, y, d_xy, n, n + 1, false)
    end

    # --- Remaining iterations. ---
    # Estimate the search cost of the next iteration from the number of rows that
    # actually had to be scanned in the previous one (which the row-skipping bound
    # makes highly variable). The initial full scan touched every row.
    last_searched = n
    for it in 2:(n - 1)
        n_rem = n - (it - 1)   # number of active nodes at the start of this iteration
        if n_rem > 2
            coeff = Tf(n_rem - 2)
            # Predicted work ~ (rows to scan) * (row length). Thread only when this is
            # large enough to amortise the spawn barrier.
            dothreads = threads_ok && last_searched * n_rem >= _WORK_PAR_MIN
            x, y, d_xy, last_searched = dothreads ?
                _dynamic_search_parallel(Q, D, R, obsolete, prev_z, n, coeff) :
                _dynamic_search(Q, D, R, obsolete, prev_z, n, coeff)
            _record_merge!(merges, heights, it, x, y, d_xy, R, idmap, n, n_rem)
            prev_z = _update!(D, R, obsolete, idmap, x, y, d_xy, n, n + it, false)
        else
            # Termination: join the two remaining nodes at their midpoint.
            a, b = _two_active(obsolete, n)
            _record_merge!(merges, heights, it, a, b, D[a, b], R, idmap, n, 2)
        end
    end

    return NJClust(merges, heights)
end

# Minimum predicted search work (rows-to-scan * row-length, i.e. distance reads)
# before a per-iteration search is run with threads. Below this the work is too
# small to amortise the thread-spawn barrier and serial execution wins. Tuned so
# that well-structured inputs (which skip most rows) never regress while large,
# search-heavy inputs still parallelise.
const _WORK_PAR_MIN = 500_000

# The one-off initial full scan is dense O(n^2) work, so it is worth threading at a
# much smaller size than the per-iteration search.
const _SCAN_PAR_MIN = 1024

# --- Row search helpers -----------------------------------------------------

# Scan row `i` (all columns active) computing its minimum join criterion, storing
# it in Q[i] and recording the winning column in rowbestj[i]. Used only for the
# initial full scan, so it omits the obsolete check.
@inline function _scan_row_full!(Q, rowbestj, D, R, n, coeff, i)
    Tf = eltype(Q)
    q_min = typemax(Tf)
    j_min = 0
    r_i = @inbounds R[i]
    @inbounds for k in (i + 1):n
        q = coeff * D[k, i] - r_i - R[k]
        if q < q_min
            q_min = q
            j_min = k
        end
    end
    @inbounds Q[i] = q_min
    @inbounds rowbestj[i] = j_min
    return nothing
end

# Search active columns `> i` of row `i`, returning (argmin column, min Q, distance
# at the min) and caching the row minimum in Q[i]. Reads D[k, i] (column i) which
# is contiguous in Julia's column-major layout.
@inline function _search_row!(Q, D, R, obsolete, n, coeff, i)
    Tf = eltype(Q)
    q_min = typemax(Tf)
    d_min = typemax(eltype(D))
    j_min = 0
    r_i = @inbounds R[i]
    @inbounds for k in (i + 1):n
        obsolete[k] && continue
        d = D[k, i]
        q = coeff * d - r_i - R[k]
        if q < q_min
            q_min = q
            d_min = d
            j_min = k
        end
    end
    @inbounds Q[i] = q_min
    return j_min, q_min, d_min
end

# --- Full scan (first iteration) --------------------------------------------

function _full_scan!(Q, D, R, n, coeff, parallel)
    rowbestj = Vector{Int}(undef, n)
    if parallel
        @threads for i in 1:n
            _scan_row_full!(Q, rowbestj, D, R, n, coeff, i)
        end
    else
        for i in 1:n
            _scan_row_full!(Q, rowbestj, D, R, n, coeff, i)
        end
    end

    # Reduce to the global minimum. Scanning ascending with a strict `<` picks the
    # smallest row index among ties, and _scan_row_full! already picked the smallest
    # column, so the result is deterministic and independent of threading.
    x = 0
    y = 0
    bestq = typemax(eltype(Q))
    @inbounds for i in 1:n
        if Q[i] < bestq
            bestq = Q[i]
            x = i
            y = rowbestj[i]
        end
    end
    return x, y, @inbounds(D[x, y])
end

# --- Dynamic search ---------------------------------------------------------

function _dynamic_search(Q, D, R, obsolete, z, n, coeff)
    # The row of the node created last iteration was never bounded, so scan it in
    # full and use it as the starting point.
    y, qbest, dbest = _search_row!(Q, D, R, obsolete, n, coeff, z)
    x = z
    nsearched = 1  # rows fully scanned (drives the next iteration's threading choice)

    @inbounds for i in 1:n
        (obsolete[i] || i == z) && continue

        if i < z
            # The pair (i, z) is not covered by row z's scan nor by row i's stale
            # bound, so evaluate it explicitly and tighten Q[i] if needed.
            q_iz = coeff * D[z, i] - R[i] - R[z]
            if q_iz < Q[i]
                Q[i] = q_iz
            end
        end

        # Core optimisation: skip rows whose lower bound already exceeds the best.
        Q[i] > qbest && continue

        j, q_ij, d_ij = _search_row!(Q, D, R, obsolete, n, coeff, i)
        nsearched += 1
        if q_ij < qbest || (q_ij == qbest && i < x)
            qbest = q_ij
            dbest = d_ij
            x = i
            y = j
        end
    end

    return x, y, dbest, nsearched
end

function _dynamic_search_parallel(Q, D, R, obsolete, z, n, coeff)
    # Scan the new node's row serially to seed the search.
    y0, q0, d0 = _search_row!(Q, D, R, obsolete, n, coeff, z)

    nt = nthreads()
    res_q = fill(q0, nt)
    res_d = fill(d0, nt)
    res_x = fill(z, nt)
    res_y = fill(y0, nt)
    res_ns = zeros(Int, nt)

    # Striped work distribution: thread t handles rows t, t+nt, t+2nt, ...
    # Each row is owned by exactly one thread, so writes to Q[i] never race.
    @threads for t in 1:nt
        lq = q0
        ld = d0
        lx = z
        ly = y0
        ns = 0
        i = t
        @inbounds while i <= n
            if !(obsolete[i] || i == z)
                if i < z
                    q_iz = coeff * D[z, i] - R[i] - R[z]
                    if q_iz < Q[i]
                        Q[i] = q_iz
                    end
                end
                if Q[i] <= lq
                    j, q_ij, d_ij = _search_row!(Q, D, R, obsolete, n, coeff, i)
                    ns += 1
                    if q_ij < lq || (q_ij == lq && i < lx)
                        lq = q_ij
                        ld = d_ij
                        lx = i
                        ly = j
                    end
                end
            end
            i += nt
        end
        res_q[t] = lq
        res_d[t] = ld
        res_x[t] = lx
        res_y[t] = ly
        res_ns[t] = ns
    end

    # Deterministic reduction: smallest Q, ties broken by smallest row index.
    bestq = res_q[1]
    bx = res_x[1]
    by = res_y[1]
    bd = res_d[1]
    @inbounds for t in 2:nt
        if res_q[t] < bestq || (res_q[t] == bestq && res_x[t] < bx)
            bestq = res_q[t]
            bx = res_x[t]
            by = res_y[t]
            bd = res_d[t]
        end
    end
    return bx, by, bd, sum(res_ns) + 1  # + 1 for the seed row z scanned above
end

# --- Distance / divergence update -------------------------------------------

# Retire node `y`, reuse slot `x` for the new node `z`, and refresh distances and
# divergences for all remaining active nodes. Returns the reused slot `z == x`.
function _update!(D, R, obsolete, idmap, x, y, d_xy, n, parent, parallel)
    @inbounds obsolete[y] = true
    z = x
    @inbounds idmap[z] = parent
    Tr = eltype(R)

    if parallel
        nt = nthreads()
        @threads for t in 1:nt
            k = t
            @inbounds while k <= n
                if !(obsolete[k] || k == x)
                    d_kx = D[k, x]
                    d_ky = D[k, y]
                    d_kz = (d_kx + d_ky - d_xy) / 2
                    D[k, z] = d_kz
                    D[z, k] = d_kz
                    R[k] = R[k] - d_kx - d_ky + d_kz
                end
                k += nt
            end
        end
        # Accumulate the new node's divergence in a fixed ascending order so the
        # result is bit-identical to the serial path (thread-order summation would
        # round differently and could flip a near-tie downstream).
        r_z = zero(Tr)
        @inbounds for k in 1:n
            (obsolete[k] || k == x) && continue
            r_z += D[k, z]
        end
        @inbounds R[z] = r_z
    else
        r_z = zero(Tr)
        @inbounds for k in 1:n
            (obsolete[k] || k == x) && continue
            d_kx = D[k, x]
            d_ky = D[k, y]
            d_kz = (d_kx + d_ky - d_xy) / 2
            D[k, z] = d_kz
            D[z, k] = d_kz
            R[k] = R[k] - d_kx - d_ky + d_kz
            r_z += d_kz
        end
        @inbounds R[z] = r_z
    end

    return z
end

# --- Merge bookkeeping ------------------------------------------------------

# Record the join of slots `x` and `y` as row `it` of the output. `n_rem` is the
# number of active nodes before this join (used for the branch-length formula).
@inline function _record_merge!(merges, heights, it, x, y, d_xy, R, idmap, n, n_rem)
    @inbounds merges[it, 1] = _mergeidx(idmap[x], n)
    @inbounds merges[it, 2] = _mergeidx(idmap[y], n)

    d = Float64(d_xy)
    β = n_rem > 2 ? (1 / (2 * (n_rem - 2))) * (Float64(R[x]) - Float64(R[y])) : 0.0
    δa = 0.5 * d + β
    δb = d - δa

    # If a branch length is negative, shift it onto the sibling branch. This does
    # not change the topology (Kuhner & Felsenstein, 1994) and matches regNJ.
    if δb < 0
        δa -= δb
        δb = 0.0
    elseif δa < 0
        δb -= δa
        δa = 0.0
    end

    @inbounds heights[it, 1] = δa
    @inbounds heights[it, 2] = δb
    return nothing
end

# Return -x for leaves (slot id <= n) or the internal-node row index otherwise.
_mergeidx(i, n) = i <= n ? -i : i - n

# Find the two remaining active slots (used at termination).
@inline function _two_active(obsolete, n)
    a = 0
    b = 0
    @inbounds for k in 1:n
        if !obsolete[k]
            if a == 0
                a = k
            else
                b = k
                break
            end
        end
    end
    return a, b
end

end
