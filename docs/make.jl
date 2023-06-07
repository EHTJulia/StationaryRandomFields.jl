using PowerSpectrumNoises
using Documenter

DocMeta.setdocmeta!(PowerSpectrumNoises, :DocTestSetup, :(using PowerSpectrumNoises); recursive=true)

makedocs(;
    modules=[PowerSpectrumNoises],
    authors="Anna Tartaglia, Kazunori Akiyama",
    repo="https://github.com/EHTJulia/PowerSpectrumNoises.jl/blob/{commit}{path}#{line}",
    sitename="PowerSpectrumNoises.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EHTJulia.github.io/PowerSpectrumNoises.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EHTJulia/PowerSpectrumNoises.jl",
    devbranch="main",
)
