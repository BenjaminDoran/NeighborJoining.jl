```@meta
CurrentModule = NeighborJoining
```

# NeighborJoining

Documentation for [NeighborJoining](https://github.com/BenjaminDoran/NeighborJoining.jl).

This package contains algorithms for [neighbor joining](https://en.wikipedia.org/wiki/Neighbor_joining)

### What is currently implemented?
* regular Neighborjoining `regNJ()`: Saitou, N. & Nei, M. The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4, 406-425 (1987).
* fast NeighborJoining `fastNJ()`: Li, J. F. A fast neighbor joining method. Genet Mol Res 14, 8733–8743 (2015).
    * This is an estimation algorithm and uses heuristics. Thus, it will not always find the exact additive tree, but it should be close 



## Installation

```
add NeighborJoining
```

## Examples

```jldoctest
julia> using NeighborJoining

julia> Dnt_leafnames = [
           "Azoto","Ecoli_K12","Mtuberculosis_H37Rv","Xylella_fastidiosa","Xac306","Pseudo","Mbovis_AF2122_97","Xfus4834",
       ];

julia> Dnt = [
           0.0       0.995214  0.985736   0.998027  0.972612  0.970401  0.985776   0.973262
           0.995214  0.0       0.998149   0.998489  0.995557  0.994368  0.998137   0.995447
           0.985736  0.998149  0.0        0.998713  0.988001  0.994798  0.0208279  0.987918
           0.998027  0.998489  0.998713   0.0       0.993199  0.998325  0.998698   0.993084
           0.972612  0.995557  0.988001   0.993199  0.0       0.986404  0.988083   0.281812
           0.970401  0.994368  0.994798   0.998325  0.986404  0.0       0.994762   0.987006
           0.985776  0.998137  0.0208279  0.998698  0.988083  0.994762  0.0        0.987954
           0.973262  0.995447  0.987918   0.993084  0.281812  0.987006  0.987954   0.0
       ];

julia> njc = regNJ(Dnt) # or fastNJ(Dnt) if regNJ is too slow
NJClust{Int64, Float64}([-7 -3; -8 -5; … ; -2 5; -4 6], [0.010421866666666738 0.010406033333333262; 0.1409980999999999 0.1408139000000001; … ; 0.4989758125 0.003104687499999925; 0.24975659374999998 0.24975659374999998])

julia> merges(njc)
7×2 Matrix{Int64}:
 -7  -3
 -8  -5
 -6  -1
  2   3
  1   4
 -2   5
 -4   6

julia> heights(njc)
7×2 Matrix{Float64}:
 0.0104219  0.010406
 0.140998   0.140814
 0.487981   0.48242
 0.347904   0.00581042
 0.485661   0.00314294
 0.498976   0.00310469
 0.249757   0.249757

julia> newickstring(njc, Dnt_leafnames)
"(Xylella_fastidiosa:2.497566e-01,(Ecoli_K12:4.989758e-01,((Mbovis_AF2122_97:1.042187e-02,Mtuberculosis_H37Rv:1.040603e-02):4.856611e-01,((Xfus4834:1.409981e-01,Xac306:1.408139e-01):3.479041e-01,(Pseudo:4.879810e-01,Azoto:4.824200e-01):5.810417e-03):3.142937e-03):3.104687e-03):2.497566e-01):0.000000e+00;"
```

