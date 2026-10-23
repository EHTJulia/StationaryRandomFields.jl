using Pkg

# Use the validated companion implementations until matching releases are available.
const COMPANIONS = (
    ("AbstractFastHartleyTransforms", "bd08b5244d829d5d2dac2acaf9e260f5ff0da7e6"),
    ("FastHartleyTransform", "a1484b1df9b90aabe95e9c412cd631afeb6f3b9f"),
)

function setup(project=dirname(@__DIR__))
    Pkg.activate(project)
    checkouts = map(COMPANIONS) do (name, _)
        joinpath(dirname(@__DIR__), ".ci", name)
    end
    if all(isdir, checkouts)
        Pkg.develop([PackageSpec(; path=path) for path in checkouts])
    else
        Pkg.add([
            PackageSpec(; url="https://github.com/EHTJulia/$(name).jl", rev=revision) for
            (name, revision) in COMPANIONS
        ])
    end
    if abspath(project) != dirname(@__DIR__)
        Pkg.develop(PackageSpec(; path=dirname(@__DIR__)))
    end
    Pkg.instantiate()
    return nothing
end

setup(isempty(ARGS) ? dirname(@__DIR__) : abspath(first(ARGS)))
