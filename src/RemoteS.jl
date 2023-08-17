module RemoteS

using GMT, Printf, Statistics, Dates#, Requires
using SatelliteToolboxTle, SatelliteToolboxPropagators, SatelliteToolboxTransformations
using PrecompileTools

const SCENE_HALFW = Dict("AQUA" => 1163479, "TERRA" => 1163479, "LANDSAT8" => 92500)	# half widths

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
	@eval Base.Experimental.@optlevel 1
end

export
	cutcube, subcube, dn2temperature, dn2radiance, dn2reflectance, reflectance_surf, grid_at_sensor, truecolor,
	clg, clre, evi, evi2, gndvi, mndwi, mtci, mcari, msavi, mbri, ndvi, ndwi, ndwi2, ndrei1,
	ndrei2, satvi, savi, slavi,
	clip_orbits, findscenes, sat_scenes, sat_tracks, reportbands

include("grid_at_sensor.jl")
include("spectral_indices.jl")
include("utils.jl")
include("sat_tracks.jl")

@setup_workload begin
	sat_tracks(tle=["1 27424U 02022A   21245.83760660  .00000135  00000-0  39999-4 0  9997"; "2 27424  98.2123 186.0654 0002229  67.6025 313.3829 14.57107527 28342"], duration=100, geocentric=true)
	sat_tracks(position=true, tle=["1 27424U 02022A   21245.83760660  .00000135  00000-0  39999-4 0  9997"; "2 27424  98.2123 186.0654 0002229  67.6025 313.3829 14.57107527 28342"]);
	truecolor(mat2img(rand(UInt16,128,128)), mat2img(rand(UInt16,128,128)), mat2img(rand(UInt16,128,128)));
	get_MODIS_scene_name(datetime2julian(DateTime("2020-09-20")), "A");
	reportbands(mat2img(rand(UInt16, 4,4,3), names=["Band 1", "Band 2", "Band 3"]), 3)[1];
end


"""
Package to perform operations with satellite data. Easy to use in computing true color images with automatic contrast stretch, many spectral indices and processing of MODIS L2 files.
"""
RemoteS

end # module
