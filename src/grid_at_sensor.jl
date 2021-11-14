"""
    G = grid_at_sensor(fname::String, sds_name::String=""; V::Bool=false, kw...)

Read one of those netCDF files that are not regular grids but have instead the coordinates in the
LONGITUDE abd LATITUDE arrays. MODIS L2 files are a good example of this. Data in theses files are
not layed down on a regular grid and we must interpolate to get one. Normally the lon and lat arrays
are called 'longitude' and 'latitude' and these it's what is seek for by default. But files exist
that pretend to comply to CF but use other names. In this case, use the kwargs `xarray` & `yarray`
to pass in the variable names. For example: `xarray="XLONG"`, `yarray="XLAT"`
The other fundamental info to pass in is the name of the array to be read/interpolated. We do that
via the `sds_name` arg.

In simpler cases the variable to be interpolated lays down on a 2D array but it is also possible that
it is stored in a 3D array. If that is the case, use the keyword 'band' to select a band (ex: 'band=2')
Bands are numbered from 1.

The interpolation is done so far with 'nearneighbor'. Both the region (-R) and increment (-I) are estimated
from data but they can be set with `region` and `inc` kwargs as well.
For MODIS data we can select the quality flag to filter by data quality. By default the best quality (=0) is
used, but one can select another with the `quality=val` kwarg. Positive 'val' values select data of quality
<= quality, whilst negative 'val' values select only data with quality >= abs(val). This allows for example
to extract only the cloud coverage.

If instead of calculating a grid (returned as a GMTgrid type) user wants the x,y,z data intself, use the
keywords `dataset`, or `outxyz` and the output will be in a GMTdataset (i.e. use `dataset=true`).

To inquire just the list of available arrays use `list=true` or `gdalinfo=true` to get the full file info.

    Examples:

    G = grid_at_sensor("AQUA_MODIS.20020717T135006.L2.SST.nc", "sst", V=true);

    G = grid_at_sensor("TXx-narr-annual-timavg.nc", "T2MAX", xarray="XLONG", yarray="XLAT", V=true);
"""
function grid_at_sensor(fname::String, sds_name::String=""; quality::Int=0, V::Bool=false, inc=0.0, kw...)

	d = KW(kw)
	(inc >= 1) && error("Silly value $(inc) for the resolution of L2 MODIS grid")
	!isfile(fname) && error("This file $(fname) does not exist")
	info = gdalinfo(fname)
	if (info === nothing)		# Try again after calling resetGMT because still something screws time-to-time
		GMT.resetGMT()
		info = gdalinfo(fname)
	end
	(haskey(d, :gdalinfo)) && (return println(info))
	((ind = findfirst("Subdatasets:", info)) === nothing) && error("This file " * fame * " is not a MODS L2 file")
	is_MODIS = (findfirst("MODISA Level-2", info) !== nothing) ? true : false
	info = info[ind[1]+12:end]		# Chop up the long string into smaller chunk where all needed info lives
	ind = findlast("SUBDATASET_", info)
	info = info[1:ind[1]]			# Chop even last SUBDATASET_X_DESC string that we wouldn't use anyway
	ind_EOLs = findall("\n", info)

	if (haskey(d, :list))
		c = ((ind = findlast("/", info[ind_EOLs[1][1] : ind_EOLs[2][1]-1])) !== nothing) ? '/' : ':'
		println("List of bands in this file:")
		[println("\t",split(info[ind_EOLs[k-1][1] : ind_EOLs[k][1]-1], c)[end]) for k = 2:2:length(ind_EOLs)]
		return nothing
	end

	(sds_name == "") && error("Must provide the band name to process. Try grid_at_sensor(\"\", list=true) to print available bands")

	# Get the arrays  SUBDATASET names
	sds_z  = helper_find_sds(sds_name, info, ind_EOLs)		# Return the full SUBDATASET name (a string)
	x_name = ((val = find_in_dict(d, [:xarray])[1]) !== nothing) ? string(val) : "longitude"
	y_name = ((val = find_in_dict(d, [:yarray])[1]) !== nothing) ? string(val) : "latitude"
	sds_qual = (is_MODIS) ? helper_find_sds("qual_" * sds_name, info, ind_EOLs) : ""
	sds_lon = helper_find_sds(x_name, info, ind_EOLs)
	sds_lat = helper_find_sds(y_name, info, ind_EOLs)

	# Get the arrays with the data
	band = ((val = find_in_dict(d, [:band])[1]) !== nothing) ? Int(val) : 1
	lon, lat, z_vals, inc, proj4 = get_xyz_qual(sds_lon, sds_lat, sds_z, quality, sds_qual, inc, band, V)

	if ((opt_R = GMT.parse_R(d, "")[1]) == "")	# If != "" believe it makes sense as a -R option
		inc_txt = split("$(inc)", '.')[2]		# To count the number of decimal digits to use in rounding
		nd = length(inc_txt)
		min_lon, max_lon = extrema(lon)
		min_lat, max_lat = extrema(lat)
		west, east   = round(min_lon-inc; digits=nd), round(max_lon+inc; digits=nd)
		south, north = round(min_lat-inc; digits=nd), round(max_lat+inc; digits=nd)
		if (x_name == "longitude")			# The "inc" pad above may "overflow" geogs. Revet the padding if need.
			((east - west) > 360) && (west += inc;	east -= inc)
			(north >  90) && (north -= inc)
			(south < -90) && (south += inc)
		end
		opt_R = @sprintf("%.10g/%.10g/%.10g/%.10g", west, east, south, north)
	else
		opt_R = opt_R[4:end]		# Because it already came with " -R....." from parse_R()
	end

	if (haskey(d, :outxyz) || haskey(d, :dataset))
		O = mat2ds([lon lat z_vals])
		O[1].proj4 = proj4
	else
		O = nearneighbor([lon lat z_vals], I=inc, R=opt_R, S=2*inc, Vd=(V) ? 1 : 0)
		O.proj4 = proj4
	end
	return O
end

# ---------------------------------------------------------------------------------------------------
function helper_find_sds(sds::String, info::String, ind_EOLs::Vector{UnitRange{Int64}})::String
	if ((ind = findfirst("/" * sds, info)) === nothing)
		((ind = findfirst(":" * sds, info)) === nothing) && error("The band name -- " * sds * " -- does not exist")
	end
	k = 1;	while (ind_EOLs[k][1] < ind[1])  k += 1  end
	return info[ind_EOLs[k-1][1]+3:ind_EOLs[k][1]-1]	# +3 because 1 = \n and the other 2 are blanks
end

# ---------------------------------------------------------------------------------------------------
function get_xyz_qual(sds_lon::String, sds_lat::String, sds_z::String, quality::Int, sds_qual::String="",
	                  inc::Float64=0., band::Int=1, V::Bool=false)
	# Get a Mx3 matrix with data to feed interpolator. Filter with quality if that's the case
	# If INC != 0, also estimates a reasonable increment for interpolation
	(V) && println("Extract lon, lat, " * sds_z * " from file")
	G = gd2gmt(sds_z; band=band);		z_vals = G.z;		proj4 = G.proj4

	if (sds_qual != "")
		Gqual = gd2gmt(sds_qual)
		if (quality >= 0)  qual = (Gqual.image .< quality + 1)		# Best (0), Best+Intermediate (1) or all (2)
		else               qual = (Gqual.image .> -quality - 1)		# Pick only, Intermediate+Lousy (-1) or Lousy (-2)
		end
		qual = reshape(qual, size(qual,1), size(qual,2))
		z_vals = z_vals[qual]
		lon, lat, dx, dy = get_lon_lat_qual(sds_lon, sds_lat, qual, inc)
		(proj4 == "") && (proj4 = GMT.prj4WGS84)	# Almost for sure that's always the case
	else
		info = gdalinfo(GMT.trim_SUBDATASET_str(sds_z))
		fill_val, got_fill_val = GMT.get_FillValue(info)
		if (got_fill_val)
			qual = (z_vals .!= fill_val)
			z_vals = z_vals[qual]
			lon, lat, dx, dy = get_lon_lat_qual(sds_lon, sds_lat, qual, inc)
		else
			G = gd2gmt(sds_lon);	lon = G.z
			G = gd2gmt(sds_lat);	lat = G.z
			(inc == 0.0) && (dx = diff(lon[:, round(Int, size(lon, 1)/2)]))
			(inc == 0.0) && (dy = diff(lat[round(Int, size(lat, 2)/2), :]))
		end
		(proj4 == "") && (proj4 = GMT.seek_wkt_in_gdalinfo(info))
	end
	(inc == 0) && (inc = guess_increment_from_coordvecs(dx, dy))
	(V) && println("Finished extraction ($(length(z_vals)) points), now intepolate")
	return lon, lat, z_vals, inc, proj4
end

function get_lon_lat_qual(sds_lon::String, sds_lat::String, qual, inc)
	# Another helper function to get only the lon, lat values that pass the 'qual' criteria
	dx, dy = Vector{Float32}(), Vector{Float32}()
	G = gd2gmt(sds_lon);
	(inc == 0.0) && (dx = diff(G.z[:, round(Int, size(G.z, 1)/2)]))
	lon = G.z[qual]
	G = gd2gmt(sds_lat);
	(inc == 0.0) && (dy = diff(G.z[round(Int, size(G.z,2)/2), :]))
	lat = G.z[qual]
	return lon, lat, dx, dy
end

# ---------------------------------------------------------------------------------------------------
function guess_increment_from_coordvecs(dx, dy)
	# Guess a good -I<inc> from the spacings in the x (lon), y(lat) arrays
	x_mean = abs(Float64(mean(dx)));		y_mean = abs(Float64(mean(dy)));
	xy_std = max(Float64(std(dx)), Float64(std(dy)))
	return (xy_std == 0) ? (x_mean + y_mean) / 2 : round((x_mean + y_mean) / 2, digits=round(Int, abs(log10(xy_std))))
end