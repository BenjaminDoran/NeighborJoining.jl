var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = NeighborJoining","category":"page"},{"location":"#NeighborJoining","page":"Home","title":"NeighborJoining","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for NeighborJoining.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [NeighborJoining]","category":"page"},{"location":"#NeighborJoining.NJClust","page":"Home","title":"NeighborJoining.NJClust","text":"NJClust(merges::Matrix{M}, heights::Matrix{H})\n\nfields:\n\nmerges is a n-1 x 2 matrix of integers: absolute values of negative integers indicate index into the distance matrix (i.e., leaves). positive integers are the index into the merge list (i.e., the kth internal node)\nheights is an n-1 x 2 matrix where each value is the distance from the left (1) or right (2) chield from its parent. Specifically heights[i,j] is the jth childs distance to the parent node, row i.\n\n\n\n\n\n","category":"type"},{"location":"#NeighborJoining.fastNJ-Tuple{AbstractMatrix{<:Number}}","page":"Home","title":"NeighborJoining.fastNJ","text":"fastNJ(d::AbstractMatrix{<:Number})\n\nfastNJ algorithm finds k independent pairs to merge and merges them for each iteration.  This is significantly faster because it is not recalculating the full pairwise Q for each pair joined.\n\nThis algorithm is nearly additive, but there are instances where the topology  '(a,(b,(c,d)))' is inferred as '((a,b), (c,d))'. This happens when  'b' is not as close to 'c' or 'd' as it is to 'a',  but 'b' is closer to the common ancestor '(c,d)' than it is to 'a'.\n\nargs:\n\nd is an n by n square symetric distance matrix\n\nreturns:\n\nNJClust struct with fields merges and heights\n\n\n\n\n\n","category":"method"},{"location":"#NeighborJoining.newickstring","page":"Home","title":"NeighborJoining.newickstring","text":"newickstring(njc::NJClust, tiplabels=AbstractVector{<:String}; labelinternalnodes=false)\nnewickstring(merges::AbstractArray{<:Integer}, heights::AbstractArray{<:AbstractFloat}, tiplabels::AbstractVector{<:String}; labelinternalnodes=false)\n\nConverts a list of merges into a newicktree formatted string.\n\nargs:\n\nnjc: is a struct that has merges and heights\nmerges and heights from a NJClust struct:\nmerges is a n-1 x 2 matrix of integers: absolute values of negative integers indicate index into the distance matrix (i.e., leaves). positive integers are the index into the merge list (i.e., the kth internal node)\nheights is an n-1 x 2 matrix where each value is the distance from the left (1) or right (2) chield from its parent. Specifically heights[i,j] is the jth childs distance to the parent node, row i.\ntiplabels: vector of string labels corresponding to the order of leaves in the distance matrix\nlabelinternalnodes: whether to generate node labels for the internal nodes. defaults to false.\n\nreturns:\n\nnewicktree formatted string\n\nexample:\n\nd = rand(10, 10)^2\nfor i in 1:size(d,1) d[i,i]=0 end;\nnjclusts = regNJ(d)\nnwstring = newicktstring(njclusts)\n\n\n\n\n\n","category":"function"},{"location":"#NeighborJoining.regNJ-Tuple{AbstractMatrix{<:Number}}","page":"Home","title":"NeighborJoining.regNJ","text":"regNJ(d::AbstractMatrix{<:Number})\n\nregNJ algorithm is the traditional NeighborJoining algorithm from \n\nSaitou, N. & Nei, M. The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4, 406-425 (1987).\n\nThis algorithm is guarenteed to infer the tree for additive distance matrices, but it does have an algorithmic complexity of O(n^3), so it can be slow for distance matrices on the order of >10³.\n\nargs:\n\nd is an n by n square symetric distance matrix\n\nreturns:\n\nNJClust struct with fields merges and heights\n\n\n\n\n\n","category":"method"}]
}