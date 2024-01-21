using NeighborJoining, BenchmarkTools, Random

load_phylip_matrix(file) = open(load_phylip_matrix, file, "r")
function load_phylip_matrix(io::IO)
    header = readline(io)
    n = parse(Int, header)
    @debug n
    taxanames = String[]
    reldists = []
    for line in readlines(io)
        taxaname, reldistline... = split(line, r"\s")
        @debug  taxaname, reldistline
        push!(taxanames, taxaname)
        push!(reldists, parse.(Float64, reldistline))
    end
    matrix = stack(reldists)
    @assert n == length(taxanames) "Incorrect number of samples, should be $n found $(length(taxanames)) lines)"
    @assert n == size(matrix, 2) "Incorrect number of columns, should be $n found $(size(matrix, 2)) lines)"
    (;n, taxanames, matrix)
end

const SUITE = BenchmarkGroup()
SUITE["RegularNeighborJoining"] = BenchmarkGroup(["RNJ", "RegularNeighborJoining", "ExactAlgorithm"])
SUITE["FastNeighborJoining"] = BenchmarkGroup(["FNJ", "FastNeighborJoining", "InexactAlgorithm"])

for i in 2 .^ (2:10)
    n, labels, Dij = load_phylip_matrix("benchmark/benchmarkdata/dist_mtx_$i.txt")
    SUITE["RegularNeighborJoining"]["regNJ_size$i"] = @benchmarkable regNJ($Dij)
    SUITE["FastNeighborJoining"]["fastNJ_size$i"] = @benchmarkable fastNJ($Dij)
end
