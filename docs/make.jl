using Documenter, RemoteS

makedocs(
    modules = [RemoteS],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true", assets = ["assets/custom.css"]),
    authors = "Joaquim Luis",
    sitename = "RemoteS.jl",
    pages = Any[
        "Gallery"                  => [
            "Aqua orbits"          => "gallery/Aqua_orbits/remotes_sat_tracks.md",
            "Aqua SST"             => "gallery/Aqua_sst/remotes_L2_SST.md",
            "Landsat 8 images"     => "gallery/L8cube_img/remotes_L8_cube_img.md",
            "Landsat 8 NDVI"       => "gallery/L8cube_ndvi/remotes_L8_NDVI.md",
            "Cloud-Native HLS dat" => "gallery/HLS/cloud-native-hls-data.md",
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