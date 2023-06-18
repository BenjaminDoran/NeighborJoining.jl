# NeighborJoining

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://BenjaminDoran.github.io/NeighborJoining.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://BenjaminDoran.github.io/NeighborJoining.jl/dev/)
[![Build Status](https://github.com/BenjaminDoran/NeighborJoining.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/BenjaminDoran/NeighborJoining.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/BenjaminDoran/NeighborJoining.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BenjaminDoran/NeighborJoining.jl)

# NeighborJoining

This package contains algorithms for [neighbor joining](https://en.wikipedia.org/wiki/Neighbor_joining)

### What is currently implemented?
* regular Neighborjoining `regNJ()`: Saitou, N. & Nei, M. The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution 4, 406-425 (1987).
* fast NeighborJoining `fastNJ()`: Li, J. F. A fast neighbor joining method. Genet Mol Res 14, 8733–8743 (2015).
    * uses heuristics, will not always find the additive tree


## Installation

```
add git@github.com:BenjaminDoran/NeighborJoining.jl.git
```

## Examples

```
# make distance matrix
d = rand(10, 10)^2
for i in 1:size(d,1) d[i,i]=0 end;

# regular neighbor joining
njclusts = regNJ(d)

# or fast neighbor joining
njclusts = fastNJ(d)

# convert to newicktree string for export to other packages
nwstring = newickstring(njclusts)

# adding labels to the leaves
labels = ["leaf $i" for i in 1:10]
nwstring = newickstring(njclusts, labels)

# adding labels to the internal nodes
nwstring = newickstring(njclusts, labels; labelinternalnodes=true)
```
