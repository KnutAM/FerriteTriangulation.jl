using FerriteTriangulation
using Documenter

DocMeta.setdocmeta!(FerriteTriangulation, :DocTestSetup, :(using FerriteTriangulation); recursive=true)

makedocs(;
    modules=[FerriteTriangulation],
    authors="Knut Andreas Meyer and contributors",
    sitename="FerriteTriangulation.jl",
    format=Documenter.HTML(;
        canonical="https://KnutAM.github.io/FerriteTriangulation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/KnutAM/FerriteTriangulation.jl",
    devbranch="main",
)
