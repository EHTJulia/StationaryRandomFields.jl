using StationaryRandomFields
using Documenter
using DocumenterVitepress
using FFTW

DocMeta.setdocmeta!(StationaryRandomFields, :DocTestSetup, :(using StationaryRandomFields); recursive=true)

makedocs(;
    modules=[StationaryRandomFields],
    authors="Anna Tartaglia, Kazunori Akiyama",
    repo="https://github.com/EHTJulia/StationaryRandomFields.jl/blob/{commit}{path}#{line}",
    sitename="StationaryRandomFields.jl",
#    format=Documenter.HTML(;
#        prettyurls=get(ENV, "CI", "false") == "true",
#        canonical="https://EHTJulia.github.io/StationaryRandomFields.jl",
#        edit_link="main",
#        assets=String[],
#    ),
    format = MarkdownVitepress(
        repo = "https://github.com/EHTJulia/StationaryRandomFields.jl",
        devurl = "dev",
        deploy_url = "https://EHTJulia.github.io/StationaryRandomFields.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "Mathematical Basis" => "math.md",
        "Tutorial" => "tutorial.md",
        "StationaryRandomFields.jl API" => "autodocs.md",
    ],
)

deploydocs(;
    repo="github.com/EHTJulia/StationaryRandomFields.jl",
    devbranch="main",
    push_preview=true,
)
