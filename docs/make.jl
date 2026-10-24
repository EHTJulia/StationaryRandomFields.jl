using StationaryRandomFields
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(
    StationaryRandomFields, :DocTestSetup, :(using StationaryRandomFields); recursive = true
)

makedocs(;
    modules = [StationaryRandomFields],
    authors = "Anna Tartaglia, Kazunori Akiyama",
    sitename = "StationaryRandomFields.jl",
    repo = "https://github.com/EHTJulia/StationaryRandomFields.jl/blob/{commit}{path}#{line}",
    format = if get(ENV, "SRF_DOCS_FORMAT", "vitepress") == "html"
        Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true", edit_link = "main")
    else
        MarkdownVitepress(; repo = "https://github.com/EHTJulia/StationaryRandomFields.jl")
    end,
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "Mathematical Basis" => "math.md",
        "Tutorial" => "tutorial.md",
        "Covariance Kernels" => "kernels.md",
        "API" => "autodocs.md",
    ],
)

if get(ENV, "SRF_DEPLOY_DOCS", "false") == "true"
    deploydocs(; repo = "github.com/EHTJulia/StationaryRandomFields.jl", devbranch = "main")
end
