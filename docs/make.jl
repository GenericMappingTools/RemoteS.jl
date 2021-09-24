using Documenter, RemoteS

makedocs(
    modules = [RemoteS],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true", assets = ["assets/custom.css"]),
    authors = "Joaquim Luis",
    sitename = "RemoteS.jl",
    pages = Any[
        "Gallery"                  => [
            "Landsat 8 imgs"       => "gallery/L8cube_img/remotes_L8_cube_img.md",
        ],
        "Index"                    => "index.md"],

    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/GenericMappingTools/RemoteS.jl.git",
	target  = "build",
    devbranch = "main",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", "devurl" => devurl],
    push_preview = true
)