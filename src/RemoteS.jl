module RemoteS

using GMT, Printf, Statistics, Dates
using PrecompileTools
using DecisionTree

const SCENE_HALFW = Dict("AQUA" => 1163479, "TERRA" => 1163479, "LANDSAT8" => 92500)	# half widths

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
	@eval Base.Experimental.@optlevel 1
end

export
	cutcube, subcube, dn2temperature, dn2radiance, dn2reflectance, reflectance_surf, grid_at_sensor, truecolor,
	clg, clre, evi, evi2, gndvi, mndwi, mtci, mcari, msavi, nbri, ndvi, ndwi, ndwi2, ndrei1,
	ndrei2, satvi, savi, slavi,
	classify, train_raster,
	clip_orbits, findscenes, sat_scenes, sat_tracks, reportbands

include("grid_at_sensor.jl")
include("spectral_indices.jl")
include("utils.jl")
include("sat_tracks.jl")

@setup_workload begin
	truecolor(mat2img(rand(UInt16,128,128)), mat2img(rand(UInt16,128,128)), mat2img(rand(UInt16,128,128)));
	get_MODIS_scene_name(datetime2julian(DateTime("2020-09-20")), "A");
	reportbands(mat2img(rand(UInt16, 4,4,3), names=["Band 1", "Band 2", "Band 3"]), 3)[1];
end


"""
Package to perform operations with satellite data. Easy to use in computing true color images with automatic contrast stretch, many spectral indices and processing of MODIS L2 files.
"""
RemoteS

end # module
