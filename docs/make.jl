using NeighborJoining
using Documenter

DocMeta.setdocmeta!(NeighborJoining, :DocTestSetup, :(using NeighborJoining); recursive=true)

makedocs(;
    modules=[NeighborJoining],
    authors="Benjamin Doran and collaborators",
    repo="https://github.com/BenjaminDoran/NeighborJoining.jl/blob/{commit}{path}#{line}",
    sitename="NeighborJoining.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://BenjaminDoran.github.io/NeighborJoining.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/BenjaminDoran/NeighborJoining.jl",
    devbranch="main",
)
