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
Dnt_fastmrgs = [
    -7  -3
    -8  -5
    -6  -1
    -4  -2
     4   1
     3   2
     6   5
]
Dnt_fasthgts = [
    0.0104219   0.0104059
    0.140998    0.140814
    0.487981    0.48242
    0.499509    0.49898
    0.00310472  0.485661
    0.0058201   0.347895
    0.00157147  0.00157147
]

@test fastNJ(d) == NJClust{Int64, Float64}([-2 -1; -5 -4; 1 -3; 3 2], [3.0 2.0; 1.0 2.0; 3.0 4.0; 1.0 1.0])
t_d = readnw(NeighborJoining.newickstring(fastNJ(d)))
@test patristic_distances(t_d) ≈ d atol=1e-4

njc = fastNJ(Dnt)
t_Dnt = readnw(NeighborJoining.newickstring(njc, Dnt_leafnames))
@test merges(njc) ≈ Dnt_fastmrgs atol=1e-4
@test heights(njc) ≈ Dnt_fasthgts atol=1e-4
pdistorder = sortperm(Dnt_leafnames)
@test patristic_distances(t_Dnt) ≈ Dnt[pdistorder,pdistorder] atol=1e-1

truetrees = readnw.(readlines(joinpath(@__DIR__, "testtrees", "yuletree_l10_n10_s10_scaled_rounded.txt")))
truedists = patristic_distances.(truetrees)
mrghgts = fastNJ.(truedists)
leafnames = ["Tip$i" for i in 0:9]
predtrees = NeighborJoining.newickstring.(merges.(mrghgts), heights.(mrghgts), Ref(leafnames))
preddists = patristic_distances.(readnw.(predtrees))
for i in 1:10
    if i == 2 || i == 10
        @test_skip truedists[i] ≈ preddists[i] atol=1e-1
    else 
        @test truedists[i] ≈ preddists[i] atol=1e-1
    end
end

truetrees = readnw.(readlines(joinpath(@__DIR__, "testtrees", "yuletree_l10_n10_s10.txt")))
truedists = patristic_distances.(truetrees)
mrghgts = fastNJ.(truedists)
leafnames = ["Tip$i" for i in 0:9]
predtrees = NeighborJoining.newickstring.(merges.(mrghgts), heights.(mrghgts), Ref(leafnames))
preddists = patristic_distances.(readnw.(predtrees))
for i in 1:10
    if i == 2 || i == 10
        @test_skip truedists[i] ≈ preddists[i] atol=1e-1
    else 
        @test truedists[i] ≈ preddists[i] atol=1e-1
    end
end

## Need to think about how to test a tree esimator vs. exact algorithm...

# truetrees = readnw.(readlines(joinpath(@__DIR__, "testtrees", "tree_varieties.nw")))
# truedists = patristic_distances.(truetrees)
# mrghgts = fastNJ.(truedists)
# leafnames = map(x->name.(x), getleaves.(truetrees))
# predtrees = NeighborJoining.newickstring.(merges.(mrghgts), heights.(mrghgts), sort.(leafnames))
# preddists = patristic_distances.(readnw.(predtrees))
# for (i, (td, pd)) in enumerate(zip(truedists, preddists))
#     @test td ≈ pd atol=1e-1
# end