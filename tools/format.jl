using JuliaFormatter

function main(args)
    all(arg -> arg == "--write", args) || error("Usage: format.jl [--write]")
    root = dirname(@__DIR__)
    formatted = format(root; overwrite="--write" in args, format_markdown=false)
    docs_formatted = format(
        joinpath(root, "docs", "src"); overwrite="--write" in args, format_markdown=true
    )
    if !(formatted && docs_formatted) && !("--write" in args)
        error("Julia files do not match the root BlueStyle configuration.")
    end
    return nothing
end

main(ARGS)
