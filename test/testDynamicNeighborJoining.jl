using NeighborJoining
using NeighborJoining: NJClust
using NewickTree
using NewickTree: Node
using Test


"""
    patristic_distances(t::Node)
shortest branch length path between all leafs with `t` as ancester
sibling nodes `i`, and `j` of parent `p` would have patristic `distance(i, p) + distance(j, p)`
"""
function patristic_distances(tree::Node)
    leaves = getleaves(tree) |> x->sort(x; by=name)
    dists = zeros(length(leaves), length(leaves))
    for j in axes(dists, 2), i in (j+1):lastindex(dists, 1)
        dists[i, j] = dists[j, i] = NewickTree.getdistance(leaves[i], leaves[j])
    end
    return dists
end

# Compare two NJClust results at the tree level (topology + branch lengths) rather
# than the raw merge/height representation, which may legitimately differ between
# implementations even when the inferred tree is identical.
function same_tree(a::NJClust, b::NJClust, names; atol=1e-4)
    ta = readnw(NeighborJoining.newickstring(a, names))
    tb = readnw(NeighborJoining.newickstring(b, names))
    return isapprox(patristic_distances(ta), patristic_distances(tb); atol)
end

## Load test data

d = [
    0  5  9  9 8
    5  0 10 10 9
    9 10  0  8 7
    9 10  8  0 3
    8  9  7  3 0
]

Dnt_leafnames = [
    "Azoto","Ecoli_K12","Mtuberculosis_H37Rv","Xylella_fastidiosa","Xac306","Pseudo","Mbovis_AF2122_97","Xfus4834",
]
Dnt = [
    0.0       0.995214  0.985736   0.998027  0.972612  0.970401  0.985776   0.973262
    0.995214  0.0       0.998149   0.998489  0.995557  0.994368  0.998137   0.995447
    0.985736  0.998149  0.0        0.998713  0.988001  0.994798  0.0208279  0.987918
    0.998027  0.998489  0.998713   0.0       0.993199  0.998325  0.998698   0.993084
    0.972612  0.995557  0.988001   0.993199  0.0       0.986404  0.988083   0.281812
    0.970401  0.994368  0.994798   0.998325  0.986404  0.0       0.994762   0.987006
    0.985776  0.998137  0.0208279  0.998698  0.988083  0.994762  0.0        0.987954
    0.973262  0.995447  0.987918   0.993084  0.281812  0.987006  0.987954   0.0
]

## input validation
@test_throws ArgumentError dynamicNJ(rand(3, 4))

## dynamicNJ is exact: an explicit representation anchor for the 5x5 example
@test dynamicNJ(d) == NJClust{Int64, Float64}([-1 -2; 1 -3; 2 -4; 3 -5], [2.0 3.0; 3.0 4.0; 2.0 2.0; 0.5 0.5])

## the result is deterministic and independent of threading
@test dynamicNJ(d) == dynamicNJ(d; parallel = false)

t_d = readnw(NeighborJoining.newickstring(dynamicNJ(d)))
@test patristic_distances(t_d) ≈ d atol=1e-4
## exact algorithm: identical tree to regNJ
@test same_tree(dynamicNJ(d), regNJ(d), string.(1:size(d, 1)))

njc = dynamicNJ(Dnt)
@test njc == dynamicNJ(Dnt; parallel = false)
t_Dnt = readnw(NeighborJoining.newickstring(njc, Dnt_leafnames))
pdistorder = sortperm(Dnt_leafnames)
@test patristic_distances(t_Dnt) ≈ Dnt[pdistorder,pdistorder] atol=1e-1
## exact algorithm: identical tree to regNJ
@test same_tree(njc, regNJ(Dnt), Dnt_leafnames)


truetrees = readnw.(readlines(joinpath(@__DIR__, "testtrees", "yuletree_l10_n10_s10_scaled_rounded.txt")))
truedists = patristic_distances.(truetrees)
mrghgts = dynamicNJ.(truedists)
leafnames = ["Tip$i" for i in 0:9]
predtrees = NeighborJoining.newickstring.(merges.(mrghgts), heights.(mrghgts), Ref(leafnames))
preddists = patristic_distances.(readnw.(predtrees))
for i in 1:10
    @test truedists[i] ≈ preddists[i] atol=1e-4
    @test mrghgts[i] == dynamicNJ(truedists[i]; parallel = false)
end

truetrees = readnw.(readlines(joinpath(@__DIR__, "testtrees", "yuletree_l10_n10_s10.txt")))
truedists = patristic_distances.(truetrees)
mrghgts = dynamicNJ.(truedists)
leafnames = ["Tip$i" for i in 0:9]
predtrees = NeighborJoining.newickstring.(merges.(mrghgts), heights.(mrghgts), Ref(leafnames))
preddists = patristic_distances.(readnw.(predtrees))
for i in 1:10
    @test truedists[i] ≈ preddists[i] atol=1e-4
    @test mrghgts[i] == dynamicNJ(truedists[i]; parallel = false)
end

truetrees = readnw.(readlines(joinpath(@__DIR__, "testtrees", "tree_varieties.nw")))
truedists = patristic_distances.(truetrees)
mrghgts = dynamicNJ.(truedists)
leafnames = map(x->name.(x), getleaves.(truetrees))
predtrees = NeighborJoining.newickstring.(merges.(mrghgts), heights.(mrghgts), sort.(leafnames))
preddists = patristic_distances.(readnw.(predtrees))
for (i, (td, pd)) in enumerate(zip(truedists, preddists))
    @test td ≈ pd atol=1e-4
    ## exact algorithm: identical tree to regNJ
    @test same_tree(mrghgts[i], regNJ(truedists[i]), sort(leafnames[i]))
end

@test order(dynamicNJ(reshape([abs(i-j) for i in 1:6 for j in 1:6], 6, 6))) == collect(1:6)
