using StationaryRandomFields
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(
    StationaryRandomFields, :DocTestSetup, :(using StationaryRandomFields); recursive = true
)

makedocs(;
    modules = [StationaryRandomFields],
    authors = "Kazunori Akiyama, Anna Tartaglia, Paul Tiede",
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
        "Mathematical Basis" => [
            "Overview" => "math.md",
            "Latent fields and inference" => "math-latent.md",
            "Spectral transforms" => "math-spectral.md",
            "Wavelet transforms" => "math-wavelet.md",
        ],
        "Tutorials" => [
            "Overview" => "tutorial.md",
            "Gaussian time series" => "tutorial-gaussian-1d.md",
            "Gaussian image textures" => "tutorial-gaussian-2d.md",
            "Daubechies wavelet signals" => "tutorial-wavelet-1d.md",
            "Daubechies wavelet images" => "tutorial-wavelet-2d.md",
            "Custom fields" => "tutorial-custom-field.md",
            "Custom transforms" => "tutorial-custom-transform.md",
        ],
        "KernelFunctions.jl extension" => "kernels.md",
        "Performance" => "performance.md",
        "API" => "autodocs.md",
    ],
)

# Preserve the downloadable measurement table in the rendered site.
cp(
    joinpath(@__DIR__, "src", "assets", "dmog-dimensions-20261024.csv"),
    joinpath(@__DIR__, "build", "final_site", "assets", "dmog-dimensions-20261024.csv");
    force = true
)

if get(ENV, "SRF_DEPLOY_DOCS", "false") == "true"
    deploydocs(; repo = "github.com/EHTJulia/StationaryRandomFields.jl", devbranch = "main")
end
