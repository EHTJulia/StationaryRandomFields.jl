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
    format = MarkdownVitepress(
        repo = "https://github.com/EHTJulia/StationaryRandomFields.jl",
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
    target = "build",
    devbranch="main",
    branch = "gh-pages",
    push_preview=true,
)
