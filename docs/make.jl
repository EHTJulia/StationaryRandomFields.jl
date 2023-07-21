using StationaryRandomFields
using Documenter

DocMeta.setdocmeta!(StationaryRandomFields, :DocTestSetup, :(using StationaryRandomFields); recursive=true)

makedocs(;
    modules=[StationaryRandomFields],
    authors="Anna Tartaglia, Kazunori Akiyama",
    repo="https://github.com/EHTJulia/StationaryRandomFields.jl/blob/{commit}{path}#{line}",
    sitename="StationaryRandomFields.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EHTJulia.github.io/StationaryRandomFields.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
    ],
)

deploydocs(;
    repo="github.com/EHTJulia/StationaryRandomFields.jl",
    devbranch="main",
)
