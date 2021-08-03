module RemoteS

export
	bright_T, grid_at_sensor, ndvi, truecolor, radiance_TOA, reflectance_TOA, reflectance_surf

include("grid_at_sensor.jl")
include("utils.jl")

end # module
