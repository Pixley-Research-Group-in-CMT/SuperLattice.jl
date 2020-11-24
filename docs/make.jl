using SuperLattice
using Documenter

makedocs(;
    modules=[SuperLattice],
    authors="Yixing Fu",
    repo="https://github.com/yixingfu/SuperLattice.jl/blob/{commit}{path}#L{line}",
    sitename="SuperLattice.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yixingfu.github.io/SuperLattice.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yixingfu/SuperLattice.jl",
)
