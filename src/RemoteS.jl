module RemoteS

using GMT, Printf, Requires, Dates

const SCENE_HALFW = Dict("AQUA" => 1163479, "TERRA" => 1163479, "LANDSAT8" => 95236)	# half widths

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
	@eval Base.Experimental.@optlevel 1
end

export
	dn2temperature, dn2radiance, dn2reflectance, reflectance_surf, grid_at_sensor, truecolor,
	clg, clre, evi, evi2, gndvi, mndwi, mtci, mcari, msavi, mbri, ndvi, ndwi, ndwi2, ndrei1,
	ndrei2, satvi, savi, slavi,
	sat_tracks, findscenes

include("grid_at_sensor.jl")
include("utils.jl")

function __init__()
	@require SatelliteToolbox="6ac157d9-b43d-51bb-8fab-48bf53814f4a" include("sat_tracks.jl")
end

end # module
