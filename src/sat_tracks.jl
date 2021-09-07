"""
sat_tracks(; tiles::Bool=false, position::Bool=false, kwargs...)

Compute satellite tracks using the TLE, or Two Line Elements set, a data format that contains
information about the orbit at a specific epoch of an Earth-orbiting object. It can also
calculate polygons arround the scene extents of AQUA and TERRA satellites as well as create
the scene names, which provides a mean to direct download that data.
  - `start`: A DateTime object or a string convertable to a DateTime with `DateTime(start)`
     specifying the start of the orbit calculation. If omited, current time in UTC will be used.
  - `duration`: Length of time for which the orbit is calculated. Accepts duration in days, hours,
     minutes or seconds. The default is minutes (100 minutes). To use other units use a string with
     the value appended with 'D', 'h', 'm' or 's'. _e.g._ `duration="55m"` to compute orbit 55 minutes from `start`
  - `step` or `inc` or `dt`: The time interval at which to compute locations along the orbit. The default
     unit here is seconds (30 sec) but minutes can be used as well by appending 'm'. _e.g._ `step="1m"`
  - `stop`: As alternative to `duration` provide the end date for the orbit. Same conditions as `start`
  - `position`: Computes only first location af the `start` time. Boolean, use `position=true`
  - `tle` or `TLE`: a file name with the TLE data for a specific satellite and period. It can also be a
    two elements string vector with the first and second lines of the TLE file.
  - `tiles`: Compute the scene limits and file names for some satellites. Currently AQUA only.

### Returns
A GMTdataset with the orbit or the scene polygons

## Example: 
Compute ~one orbit of the AQUA satellite starting at current local time. Note, this
will be accurate for the month of September 2021. For other dates it needs an updated TLE.

    tle1 = "1 27424U 02022A   21245.83760660  .00000135  00000-0  39999-4 0  9997";
    tle2 = "2 27424  98.2123 186.0654 0002229  67.6025 313.3829 14.57107527 28342";
    orb = sat_tracks(tle=[tle1; tle2], duration=100);

and the orbit track can be visualized with

    imshow(orb,  proj=:Robinson, region=:global, coast=true)
"""
function sat_tracks(; geocentric::Bool=false, tiles::Bool=false, position::Bool=false, kwargs...)
	# ...
	(position && tiles) && error("Cannot require tiles and a single position. Makes no sense.")
	d = KW(kwargs)

	function getitDTime()
		if     isa(val, String) || isa(val, Tuple)  ret::DateTime = DateTime(val)
		elseif (isa(val, DateTime))  ret = val
		else   error("Bad input type ($(typeof(val)). Must be a DateTime, a String or a Tuple(Int)")
		end
		return ret
	end

	start = ((val = find_in_dict(d, [:start])[1]) === nothing) ? now(Dates.UTC) : getitDTime()
	
	(tiles) && (start = round(start, Dates.Minute(5)))	# This is for MODIS only

	if ((val = find_in_dict(d, [:duration])[1]) !== nothing)
		if (isa(val, String))		# Accept duration in D(ays), h(ours), m(inutes) or s(econds)
			if     (endswith(val,"D"))  dur = Day(parse(Int, val[1:end-1]))
			elseif (endswith(val,"h"))  dur = Hour(parse(Int, val[1:end-1]))
			elseif (endswith(val,"m"))  dur = Minute(parse(Int, val[1:end-1]))
			elseif (endswith(val,"s"))  dur = Second(parse(Int, val[1:end-1]))
			else    error("Only 'D', 'h', 'm' or 's' are accepted in duration")
			end
		elseif (isa(val, Real))		# Assume duration was given in hours
			dur = Hour(trunc(Int, val))
		end
		stop = start + dur
	else
		if ((val = find_in_dict(d, [:stop :end])[1]) !== nothing)
			stop = getitDTime()
		else
			stop = start + Minute(100)		# Default is ~Terra period
		end
	end
	if ((val = find_in_dict(d, [:step :inc :dt])[1]) !== nothing)	# Steps are in seconds
		if (isa(val, String))
			if     (endswith(val,"m"))  dt = parse(Int, val[1:end-1]) * 60
			elseif (endswith(val,"s"))  dt = parse(Int, val[1:end-1])
			else    error("Only 's' or 'm' are accepted in increment")
			end
		else
			dt = trunc(Int, val)
		end
	else
		dt = 30
	end

	if (tiles)
		#if ((val = find_in_dict(d, [:AQUA])[1]) !== nothing)
			dt, halfwidth = 60, 2326958 / 2	# This is if for AQUA (and got by measuring over a L2 grid
			sat_name = "AQUA"
		#end
	end

	if ((val = find_in_dict(d, [:tle :TLE])[1]) !== nothing)
		if (isa(val, String))  tle = SatelliteToolbox.read_tle(val)
		elseif (isa(val, Vector{String}) && length(val) == 2)
			tle = SatelliteToolbox.read_tle_from_string(val[1], val[2])
		else
			error("BAD input TLE data")
		end
	else
		#tle = SatelliteToolbox.read_tle("C:\\v\\Landsat8.tle")
		tle = SatelliteToolbox.read_tle("C:\\v\\AQUA.tle")
	end

	orbp = SatelliteToolbox.init_orbit_propagator(Val(:sgp4), tle[1])

	startmfe = (datetime2julian(DateTime(start)) - tle[1].epoch) * 24 * 3600
	stopmfe  = (datetime2julian(DateTime(stop))  - tle[1].epoch) * 24 * 3600
	t = startmfe:dt:stopmfe
	
	(position) && (t = [t[1]])		# Single position. Doing it here wastes work above but code is way cleaner

	out = Matrix{Float64}(undef, length(t), 4)
	r, = SatelliteToolbox.propagate!(orbp, t)

	for n = 1:length(t)
		jd = tle[1].epoch + t[n] / (24 * 3600)
		#jd = datetime2julian(DateTime(start)) + (n-1)*dt / (24 * 3600)
		tt = SatelliteToolbox.r_eci_to_ecef(SatelliteToolbox.TEME(), SatelliteToolbox.PEF(), jd) * r[n]
		out[n,1], out[n,2], out[n,3], out[n, 4] = tt[1], tt[2], tt[3], jd
	end

	if (tiles)
		out = mapproject(out, E=true, I=true)
		return make_sat_tiles(out[1].data, halfwidth, sat_name)
	end

	return (geocentric) ? out : mapproject(out, E=true, I=true)[1]
end

# --------------------------------------------------------------------------------------------
function sat_scenes(track, halfwidth, sat_name)
	# Compute polygons with the scene limits. 'halfwidth' is half the scene width.
	_, azim, = invgeod(track[1:end-1, 1:2], track[2:end, 1:2])	# distances and azimuths along the tracks
	append!(azim, azim[end])		# To make it same size of tracks
	D = Vector{GMTdataset}(undef, trunc(Int, (length(azim)-1)/5))
	n = 0
	for k = 1:5:length(azim)-5
		ll14 = geod(track[k,1:2],   [azim[k]+90, azim[k]-90], halfwidth)[1]
		ll23 = geod(track[k+5,1:2], [azim[k+5]+90, azim[k+5]-90], halfwidth)[1]
		sc = get_MODIS_scene_name(track[k,4], sat_name)
		D[n+=1] = GMTdataset([ll14[1:1,:]; ll23[1:1,:]; ll23[2:2,:]; ll14[2:2,:]; ll14[1:1,:]], String[], sc, String[], "", "", GMT.Gdal.wkbPolygon)
	end
	D[1].proj4 = geo_proj4
	D
end

# --------------------------------------------------------------------------------------------
function get_MODIS_scene_name(jd::Float64, prefix::String, sst::Bool=true)
	# Build the scene name after it's acquisition time coded in the JD argin.
	DT = julian2datetime(jd)
	if (sst)
		p = (prefix[1] == 'A') ? "AQUA_MODIS." : "TERRA_MODIS."
		@sprintf("%s%.4d%.2d%.2dT%.2d%.2d01.L2.SST.NRT.nc", p, year(DT), month(DT), day(DT), hour(DT), minute(DT))
	else
		@sprintf("%s%.4d%.3d%.2d%.2d00.L2_LAC_OC.nc", prefix[1], year(DT), dayofyear(DT), hour(DT), minute(DT))
	end
end

# --------------------------------------------------------------------------------------------
function within_BB(track, bb::Vector{<:Real})
	# Select the parts of the track that are contained inside the BB.
	segments = GMT.gmtselect(track, region=bb, f=:g)
	azim = invgeod(segments[1].data[1:end-1, 1:2], segments[1].data[2:end, 1:2])[2]
	inds = findall(abs.(diff(azim)) .> 10)			# Detect line breaks

	begin_seg = [1; inds[2:2:length(inds)].+1]		# Not easy to explain why those indices give seg boundaries
	end_seg   = [inds[2:2:length(inds)]; length(azim)+1]
	D = Vector{GMTdataset}(undef, length(begin_seg))
	for k = 1:length(begin_seg)
		D[k] = GMTdataset(segments[1].data[begin_seg[k]:end_seg[k], :], String[], "", String[], "", "", GMT.Gdal.wkbLineString)
	end
	D[1].proj4 = GMT.geo_proj4
	D
end

# --------------------------------------------------------------------------------------------
function findscene(track::Union{GMTdataset, Vector{GMTdataset}}, lon::Real, lat::Real, distmin::Real=Inf)
	# ...
	dists = mapproject(track, G=(lon,lat))
	for k = 1:length(dists)
		d, ind = findmin(view(dists[k].data, :,5))
		if (d <= distmin)
			jd = datetime2julian(floor(julian2datetime(dists[k].data[ind,4]), Dates.Minute(5)))
			sc = get_MODIS_scene_name(jd, "AQUA", false)
			println(sc)
		end
	end
end

findscene(date::String, lon::Real, lat::Real, distmin::Real=Inf) = findscene(DateTime(date), lon, lat, distmin)
function findscene(date::DateTime, lon::Real, lat::Real, distmin::Real=Inf)
	# ...
	_date = date - Dates.Day(1)
	tracks = sat_tracks(start=_date, duration="2D")
	BB = [lon-10, lon+10, lat-10, lat+10]
	D = within_BB(tracks, BB)
	findscene(D, lon, lat, 1163479)
end