module RemoteS

using GMT, Printf, Requires, Dates

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
	@eval Base.Experimental.@optlevel 1
end

export
	dn2temperature, dn2radiance, dn2reflectance, reflectance_surf, grid_at_sensor, truecolor,
	clg, clre, evi, evi2, gndvi, mndwi, mtci, mcari, msavi, mbri, ndvi, ndwi, ndwi2, ndrei1,
	ndrei2, satvi, savi, slavi,
	sat_tracks, findscene

include("grid_at_sensor.jl")
include("utils.jl")

function __init__()
	@require SatelliteToolbox="6ac157d9-b43d-51bb-8fab-48bf53814f4a" include("sat_tracks.jl")
end

end # module
