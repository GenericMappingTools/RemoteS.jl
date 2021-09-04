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

	function get_MODIS_scene_name(jd)
		DT = julian2datetime(jd)
		@sprintf("A%.4d%.3d%.2d%.2d%.2d.L2", year(DT), dayofyear(DT), hour(DT), minute(DT), second(DT))
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
		#if ((val = find_in_dict(d, [:MODIS])[1]) !== nothing)
			dt, halfwidth = 60, 2326958 / 2	# This is if for AQUA (and got by measuring over a L2 grid
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
		#@show(get_MODIS_scene_name(jd))
	end

	if (tiles)
		out = mapproject(out, E=true, I=true)
		return make_sat_tiles(out[1].data, halfwidth)
	end

	return (geocentric) ? out : mapproject(out, E=true, I=true)[1]
end

# --------------------------------------------------------------------------------------------
function make_sat_tiles(track, halfwidth)
	_, azim, = invgeod(track[1:end-1, 1:2], track[2:end, 1:2])	# distances and azimuths along the tracks
	append!(azim, azim[end])		# To make it same size of tracks
	D = Vector{GMTdataset}(undef, trunc(Int, (length(azim)-1)/5))
	n = 0
	for k = 1:5:length(azim)-5
		ll14 = geod(track[k,1:2],   [azim[k]+90, azim[k]-90], halfwidth)[1]
		ll23 = geod(track[k+5,1:2], [azim[k+5]+90, azim[k+5]-90], halfwidth)[1]
		D[n+=1] = GMTdataset([ll14[1:1,:]; ll23[1:1,:]; ll23[2:2,:]; ll14[2:2,:]; ll14[1:1,:]])
	end
	D[1].proj4 = "geog"
	D
end