using Documenter, RemoteS, GMT

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
    #repo = "github.com/GenericMappingTools/RemoteS.jl.git",
    repo = "https://www.generic-mapping-tools.org/RemoteS.jl",
	target  = "build",
    versions = ["stable" => "v^", "v#.#", devurl => devurl],
    push_preview = true
)