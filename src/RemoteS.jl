module RemoteS

using GMT, Printf

export
	dn2temperature, dn2radiance, dn2reflectance, reflectance_surf, grid_at_sensor, truecolor,
	clg, clre, evi, evi2, gndvi, mndwi, mtci, mcari, msavi, mbri, ndvi, ndwi, ndwi2, ndrei1,
	ndrei2, satvi, savi, slavi

include("grid_at_sensor.jl")
include("utils.jl")

end # module
