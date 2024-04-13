module RemoteSSatTbxExt
	using RemoteS, GMT, Dates
	using SatelliteToolboxTle, SatelliteToolboxPropagators, SatelliteToolboxTransformations
	using PrecompileTools

	function RemoteS.sat_tracks_ext(; geocentric::Bool=false, tiles::Bool=false, position::Bool=false, kwargs...)
		(position && tiles) && error("Cannot require tiles and a single position. Makes no sense.")
		d = GMT.KW(kwargs)
	
		start = ((val = find_in_dict(d, [:start])[1]) === nothing) ? now(Dates.UTC) : RemoteS.getitDTime(val)
		
		(tiles) && (sat_name = get_sat_name(d))
		(tiles) && (start = round(start, Dates.Minute(5)))	# This is for MODIS only
	
		if ((val = find_in_dict(d, [:duration])[1]) !== nothing)
			if (isa(val, String))		# Accept duration in D(ays), h(ours), m(inutes) or s(econds)
				if     (endswith(val,"D"))  dur = Day(parse(Int, val[1:end-1]))
				elseif (endswith(val,"h"))  dur = Hour(parse(Int, val[1:end-1]))
				elseif (endswith(val,"m"))  dur = Minute(parse(Int, val[1:end-1]))
				elseif (endswith(val,"s"))  dur = Second(parse(Int, val[1:end-1]))
				else    error("Only 'D', 'h', 'm' or 's' are accepted in duration")
				end
			elseif (isa(val, Real))		# Assume duration was given in minutes
				dur = Minute(trunc(Int, val))
			end
			stop = start + dur
		else
			if ((val = find_in_dict(d, [:stop])[1]) !== nothing)
				stop = RemoteS.getitDTime(val)
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
	
		if ((val_tle = find_in_dict(d, [:tle_obj])[1]) !== nothing)  tle = val_tle	# Some other fun already got it
		else                                                         tle = RemoteS.loadTLE_ext(d)
		end
		orbp = SatelliteToolboxPropagators.Propagators.init(Val(:SGP4), tle)
	
		epoch_jd = orbp.sgp4d.epoch
		startmfe = (datetime2julian(DateTime(start)) - epoch_jd) * 24 * 3600
		stopmfe  = (datetime2julian(DateTime(stop))  - epoch_jd) * 24 * 3600
		(tiles) && (dt = 60)			# Arbitrary choice that works well for MODIS but may need revision for others
		t = startmfe:dt:stopmfe
		
		(position) && (t = [t[1]])		# Single position. Doing it here wastes work above but code is way cleaner
	
		out = Matrix{Float64}(undef, length(t), 4)
		#r = SatelliteToolboxPropagators.Propagators.propagate!.(orbp, t)[1]	# Doesn't work. Why?
	
		for n = 1:GMT.numel(t)
			r = SatelliteToolboxPropagators.Propagators.propagate!(orbp, t[n])[1]
			jd = epoch_jd + t[n] / (24 * 3600)
			tt = SatelliteToolboxTransformations.r_eci_to_ecef(SatelliteToolboxTransformations.TEME(), SatelliteToolboxTransformations.PEF(), jd) * r
			out[n,1], out[n,2], out[n,3], out[n, 4] = tt[1], tt[2], tt[3], jd
		end
	
		if (tiles)
			out = mapproject(out, E=true, I=true)
			return make_sat_tiles(out[1].data, SCENE_HALFW[sat_name], sat_name)
		end
	
		if (geocentric)
			D = mat2ds(out, geom=UInt32(2), colnames=["X", "Y", "Z", "JulianDay"])
		else
			D = mapproject(out, E=true, I=true)
			D.colnames, D.proj4, D.geom = ["lon", "lat", "alt", "JulianDay"], GMT.prj4WGS84, UInt32(2)
		end
		D
	end
	
	# --------------------------------------------------------------------------------------------
	function RemoteS.loadTLE_ext(d::Dict)
		# Load a TLE or use a default one. In a function because it's used by two functions
		if ((val = find_in_dict(d, [:tle :TLE])[1]) !== nothing)
			if (isa(val, String))  tle = SatelliteToolboxTle.read_tle(val)
			elseif (isa(val, Vector{String}) && length(val) == 2)
				tle = SatelliteToolboxTle.read_tle(val[1], val[2])
			else
				error("BAD input TLE data")
			end
		else
			tle = SatelliteToolboxTle.read_tle("C:\\v\\AQUA.tle")
		end
		return tle
	end

	@setup_workload begin
		RemoteS.sat_tracks_ext(tle=["1 27424U 02022A   21245.83760660  .00000135  00000-0  39999-4 0  9997"; "2 27424  98.2123 186.0654 0002229  67.6025 313.3829 14.57107527 28342"], duration=100, geocentric=true)
		RemoteS.sat_tracks_ext(position=true, tle=["1 27424U 02022A   21245.83760660  .00000135  00000-0  39999-4 0  9997"; "2 27424  98.2123 186.0654 0002229  67.6025 313.3829 14.57107527 28342"]);
	end

end