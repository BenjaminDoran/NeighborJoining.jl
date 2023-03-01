using NeighborJoining
using NeighborJoining: NJClust
using NewickTree
using NewickTree: Node
using LinearAlgebra: Symmetric
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
        dists[i, j] = NewickTree.getdistance(leaves[i], leaves[j])
    end
    return Symmetric(dists, :L)
end

@testset "NeighborJoining.jl" begin
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
        0.0                0.995213839457987  0.98573559031375    0.998027348864072   0.972611680703211  0.970401190581917  0.985776054066666  0.973262192179632
        0.995213839457987  0.0                0.998148938539192   0.998488816108287   0.995557411172427  0.994367546028565  0.998137116752552  0.995446830617351
        0.98573559031375   0.998148938539192  0.0                 0.99871270422122    0.988000892034872  0.9947982723887    0.0208278520191335 0.987917785278169
        0.998027348864072  0.998488816108287  0.99871270422122    0.0                 0.993199212089598  0.998324711025772  0.998697800631595  0.993083884624155
        0.972611680703211  0.995557411172427  0.988000892034872   0.993199212089598   0.0                0.986404377162655  0.988083283435624  0.28181153134629
        0.970401190581917  0.994367546028565  0.9947982723887     0.998324711025772   0.986404377162655  0.0                0.994762168983286  0.987006080464416
        0.985776054066666  0.998137116752552  0.0208278520191335  0.998697800631595   0.988083283435624  0.994762168983286  0.0                0.987953833947853
        0.973262192179632  0.995446830617351  0.987917785278169   0.993083884624155   0.28181153134629   0.987006080464416  0.987953833947853  0.0  
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
    Dnt_regmrgs = [
        -7  -3
        -8  -5
        -6  -1
         3   2
         4   1
        -4  -2
         6   5
    ]
    Dnt_reghgts = [
        0.0104219   0.0104059
        0.140998    0.140814
        0.487981    0.48242
        0.0058103   0.347904
        0.00314293  0.485661
        0.499513    0.498976
        0.00155236  0.00155236
    ]


    @testset "RegNJ" begin 
        @test regNJ(d) == NJClust([-2 -1; 1 -3; -5 -4; 3 2], [3.0 2.0; 3.0 4.0; 1.0 2.0; 1.0 1.0])
        njc = regNJ(Dnt)
        @test merges(njc) ≈ Dnt_regmrgs atol=1e-4
        @test heights(njc) ≈ Dnt_reghgts atol=1e-4

        truetrees = readnw.(open(readlines, "./testtrees/yuletree_l10_n10_s10_scaled_rounded.txt"))
        truedists = patristic_distances.(truetrees)
        mrghgts = regNJ.(truedists)
        leafnames = ["Tip$i" for i in 0:9]
        predtrees = newickstring.(merges.(mrghgts), heights.(mrghgts), Ref(leafnames))
        preddists = patristic_distances.(readnw.(predtrees))
        for i in 1:10
            @test truedists[i] ≈ preddists[i] atol=1e-4
        end

        truetrees = readnw.(open(readlines, "./testtrees/yuletree_l10_n10_s10.txt"))
        truedists = patristic_distances.(truetrees)
        mrghgts = regNJ.(truedists)
        leafnames = ["Tip$i" for i in 0:9]
        predtrees = newickstring.(merges.(mrghgts), heights.(mrghgts), Ref(leafnames))
        preddists = patristic_distances.(readnw.(predtrees))
        for i in 1:10
            @test truedists[i] ≈ preddists[i] atol=1e-4
        end

    end
    @testset "fastNJ" begin 
        @test fastNJ(d) == NJClust([-2 -1; -5 -4; 1 -3; 3 2], [3.0 2.0; 1.0 2.0; 3.0 4.0; 1.0 1.0])
        njc = fastNJ(Dnt)
        @test merges(njc) ≈ Dnt_fastmrgs atol=1e-4
        @test heights(njc) ≈ Dnt_fasthgts atol=1e-4

        truetrees = readnw.(open(readlines, "./testtrees/yuletree_l10_n10_s10_scaled_rounded.txt"))
        truedists = patristic_distances.(truetrees)
        mrghgts = fastNJ.(truedists)
        leafnames = ["Tip$i" for i in 0:9]
        predtrees = newickstring.(merges.(mrghgts), heights.(mrghgts), Ref(leafnames))
        preddists = patristic_distances.(readnw.(predtrees))
        for i in 1:10
            if i == 2 || i == 10
                @test_skip truedists[i] ≈ preddists[i] atol=1e-1
            else 
                @test truedists[i] ≈ preddists[i] atol=1e-1
            end
        end

        truetrees = readnw.(open(readlines, "./testtrees/yuletree_l10_n10_s10.txt"))
        truedists = patristic_distances.(truetrees)
        mrghgts = fastNJ.(truedists)
        leafnames = ["Tip$i" for i in 0:9]
        predtrees = newickstring.(merges.(mrghgts), heights.(mrghgts), Ref(leafnames))
        preddists = patristic_distances.(readnw.(predtrees))
        for i in 1:10
            if i == 2 || i == 10
                @test_skip truedists[i] ≈ preddists[i] atol=1e-1
            else 
                @test truedists[i] ≈ preddists[i] atol=1e-1
            end
        end
    end
    @testset "newickstring" begin 
        njc = fastNJ(d)
        ptree = readnw(newickstring(merges(njc), heights(njc), ["a", "b", "c", "d", "e"]))
        truetree = nw"(((b:3.0,a:2.0):3.0,c:4.0):1.0,(e:1.0,d:2.0):1.0);"
        @test patristic_distances(ptree) ≈ patristic_distances(truetree) atol=1e-6
    end
end
