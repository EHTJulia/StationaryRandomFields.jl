using Runic

function main(args)
    all(arg -> arg == "--write", args) || error("Usage: format.jl [--write]")
    root = dirname(@__DIR__)
    # Include new files while excluding ignored builds and companion checkouts.
    files = cd(root) do
        split(read(`git ls-files -z --cached --others --exclude-standard -- '*.jl' '*.md'`, String), '\0'; keepempty = false)
    end
    paths = [joinpath(root, file) for file in files if isfile(joinpath(root, file))]
    isempty(paths) && error("No Julia or Markdown files found.")
    mode = "--write" in args ? ["--inplace"] : ["--check", "--diff"]
    return Runic.main(vcat(mode, ["--docstrings", "--languages=julia,julia-repl,jldoctest,@example,@repl,@setup,@eval"], paths))
end

exit(main(ARGS))
