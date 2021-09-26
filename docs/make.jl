using Documenter, RemoteS

makedocs(
    modules = [RemoteS],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true", assets = ["assets/custom.css"]),
    authors = "Joaquim Luis",
    sitename = "RemoteS.jl",
    pages = Any[
        "Gallery"                  => [
            "Aqua orbits"          => "gallery/Aqua_orbits/remotes_sat_tracks.md",
            "Landsat 8 imgs"       => "gallery/L8cube_img/remotes_L8_cube_img.md",
            "Landsat 8 NDVI"       => "gallery/L8cube_ndvi/remotes_L8_NDVI.md",
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
    #versions = ["stable" => "v^", "v#.#", devurl => devurl],
    push_preview = true
)