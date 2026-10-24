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
    format = MarkdownVitepress(;
        repo = "https://github.com/EHTJulia/StationaryRandomFields.jl",
        devbranch = "main",
        description = "Correlated noise and latent-field priors for astronomical signals and images",
    ),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "Installation" => "installation.md",
        "Mathematical Basis" => "math.md",
        "Tutorial" => "tutorial.md",
        "Covariance Kernels" => "kernels.md",
        "API" => "autodocs.md",
    ],
)

if get(ENV, "SRF_DEPLOY_DOCS", "false") == "true"
    deploydocs(; repo = "github.com/EHTJulia/StationaryRandomFields.jl", devbranch = "main")
end
