using Documenter, RemoteS

makedocs(
    modules = [RemoteS],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Joaquim Luis",
    sitename = "RemoteS.jl",
    pages = Any["Index"                    => "index.md"],

    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/GenericMappingTools/RemoteS.jl.git",
    push_preview = true
)