using Pkg

function setup(project = dirname(@__DIR__))
    root = dirname(@__DIR__)

    Pkg.activate(project)

    if abspath(project) != root
        Pkg.develop(PackageSpec(; path = root))
    end

    Pkg.instantiate()
    return nothing
end

setup(isempty(ARGS) ? dirname(@__DIR__) : abspath(first(ARGS)))
