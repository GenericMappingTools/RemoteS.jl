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

	#tle = SatelliteToolbox.read_tle("C:\\v\\Landsat8.tle")
	tle = SatelliteToolbox.read_tle("C:\\v\\AQUA.tle")

	#tle1 = "1 27424U 02022A   21245.83760660  .00000135  00000-0  39999-4 0  9997"
	#tle2 = "2 27424  98.2123 186.0654 0002229  67.6025 313.3829 14.57107527 28342"
	#tle = SatelliteToolbox.read_tle_from_string(tle1, tle2, true)

	orbp = SatelliteToolbox.init_orbit_propagator(Val(:sgp4), tle[1])

#=
	out, = propagate_to_epoch!(orbp, tle[1].epoch .+ collect(0.:60:60*100) ./ (3600*24))	# StaticArray
	out = collect(transpose(reshape(reinterpret(Float64, out), 3, size(out,1))))	# Now a Matrix
=#

	startmfe = (datetime2julian(DateTime(start)) - tle[1].epoch) * 24 * 3600
	stopmfe  = (datetime2julian(DateTime(stop))  - tle[1].epoch) * 24 * 3600
	t = startmfe:dt:stopmfe
	
	(position) && (t = [t[1]])		# Single position. Doing it here wastes work above but code is way cleaner

	out = Matrix{Float64}(undef, length(t), 4)
	r, = SatelliteToolbox.propagate!(orbp, t)
	#eop_IAU1980 = SatelliteToolbox.get_iers_eop()

	#convrt = pi / (3600 * 180)		# arc sec to rad
	#xp, yp = 0.0075 * convrt, 0.2925 * convrt		# xp and yp from Naval observatory, 1/26/2016

	for n = 1:length(t)
		jd = tle[1].epoch + t[n] / (24 * 3600)
		#jd = datetime2julian(DateTime(start)) + (n-1)*dt / (24 * 3600)
		##
		tt = SatelliteToolbox.r_eci_to_ecef(SatelliteToolbox.TEME(), SatelliteToolbox.PEF(), jd) * r[n]
		out[n,1], out[n,2], out[n,3], out[n, 4] = tt[1], tt[2], tt[3], jd
		@show(get_MODIS_scene_name(jd))
		#=
		DT = julian2datetime(jd)
		Y,M,D = yearmonthday(DT)
		jdut1, ttt = convtime(Y,M,D, Dates.hour(DT), Dates.minute(DT), Dates.second(DT), 0, 0, 0)
		recef = teme2ecef([r[n][1]; r[n][2]; r[n][3]], ttt, jdut1, 0, xp, yp)
		out[n,1], out[n,2], out[n,3], out[n, 4] = recef[1], recef[2], recef[3], jd
		=#
	end

	if (tiles)
		out = mapproject(out, E=true, I=true)
		return make_sat_tiles(out[1].data, halfwidth)
	end

	return (geocentric) ? out : mapproject(out, E=true, I=true)[1]
end

#=
# Dates.dayofyear(Dates.now())
# ju = julian2datetime(2459457.24647)
# 2021-08-30T17:54:55.008
# Dates.year(ju); Dates.second(ju); yearmonthday(ju)
=#

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
	D[1]
end


# --------------------------------------------------------------------------------------------
function teme2ecef(rteme,ttt,jdut1,lod,xp,yp)
	gmst = gstime(jdut1)

	st      = zeros(3,3);
	st[1,1] =  cos(gmst);
	st[1,2] = -sin(gmst);
	st[2,1] =  sin(gmst);
	st[2,2] =  cos(gmst);
	st[3,3] =  1.0;

	pm = polarm(xp, yp, ttt, "80");
	rpef  = st' * rteme;
	recef = (pm' * rpef)';

end

# --------------------------------------------------------------------------------------------
function convtime(year, mon, day, hr, min, sec, timezone, dut1, dat)
#  this function finds the time parameters and julian century values for inputs
#    of utc or ut1. numerous outputs are found as shown in the local variables.
#    because calucations are in utc, you must include timezone if ( you enter a
#    local time, otherwise it should be zero.
#
#  author        : david vallado                  719-573-2600    4 jun 2002
#
#  revisions
#    vallado     - add tcg, tcb, etc                              6 oct 2005
#    vallado     - fix documentation for dut1                     8 oct 2002
#
#  inputs          description                    range / units
#    year        - year                           1900 .. 2100
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28,29,30,31
#    hr          - universal time hour            0 .. 23
#    min         - universal time min             0 .. 59
#    sec         - universal time sec (utc)            0.0  .. 59.999
#    timezone    - offset to utc from local site  0 .. 23 hr
#    dut1        - delta of ut1 - utc             sec
#    dat         - delta of tai - utc             sec

	# ------------------ start if ( ut1 is known ------------------
	localhr = timezone + hr
	utc = localhr * 3600 + min * 60 + sec

	ut1 = utc + dut1
	hrtemp,mintemp,sectemp = sec2hms(ut1)
	#jdut1 = jday(year,mon,day, hrtemp, mintemp, sectemp)
	jdut1 = datetime2julian(DateTime(year,mon,day, hrtemp, mintemp, sectemp))
	#tut1= (jdut1 - 2451545.0) / 36525.0

	tai = utc + dat
	tt = tai + 32.184;   # sec
	hrtemp,mintemp,sectemp = sec2hms(tt)
	#jdtt = jday(year,mon,day, hrtemp, mintemp, sectemp)
	jdtt = datetime2julian(DateTime(year,mon,day, hrtemp, mintemp, sectemp))
	ttt= (jdtt - 2451545.0) / 36525.0

	return jdut1, ttt
end

# --------------------------------------------------------------------------------------------
function sec2hms(utsec)
#  this function converts seconds from the beginning of the day into hours, minutes and seconds.
	temp  = utsec / 3600.0
	hr    = trunc(Int,temp)
	min   = trunc(Int, (temp - hr)* 60.0 )
	sec   = trunc(Int, (temp - hr - min/60.0 ) * 3600.0)	# NOT RIGHT, SHOULD ROUND AND CHECK IF == 60
	return hr,min,sec
end

# --------------------------------------------------------------------------------------------
function polarm(xp, yp, ttt, opt)
	#  this function calulates the transformation matrix that accounts for polar
	#    motion. both the 1980 and 2000 theories are handled. note that the rotation 
	#    order is different between 1980 and 2000 .
	#
	#  author        : david vallado                  719-573-2600   25 jun 2002
	#
	#  revisions
	#    vallado     - consolidate with iau 2000                     14 feb 2005
	#
	#  inputs          description                    range / units
	#    xp          - polar motion coefficient       rad
	#    yp          - polar motion coefficient       rad
	#    ttt         - julian centuries of tt (00 theory only)
	#    opt         - method option                  '01', '02', '80'
	#
	#  outputs       :
	#    pm          - transformation matrix for ecef - pef

	cosxp = cos(xp);
	sinxp = sin(xp);
	cosyp = cos(yp);
	sinyp = sin(yp);

	pm = zeros(3,3)
	if (opt == "80")
		pm[1,1] =  cosxp;
		pm[1,2] =  0.0;
		pm[1,3] = -sinxp;
		pm[2,1] =  sinxp * sinyp;
		pm[2,2] =  cosyp;
		pm[2,3] =  cosxp * sinyp;
		pm[3,1] =  sinxp * cosyp;
		pm[3,2] = -sinyp;
		pm[3,3] =  cosxp * cosyp;
	else  
		convrt = pi / (3600.0*180.0);
		# approximate sp value in rad
		sp = -47.0e-6 * ttt * convrt;
		cossp = cos(sp);
		sinsp = sin(sp);

		# form the matrix
		pm[1,1] =  cosxp * cossp;
		pm[1,2] = -cosyp * sinsp + sinyp * sinxp * cossp;
		pm[1,3] = -sinyp * sinsp - cosyp * sinxp * cossp;
		pm[2,1] =  cosxp * sinsp;
		pm[2,2] =  cosyp * cossp + sinyp * sinxp * sinsp;
		pm[2,3] =  sinyp * cossp - cosyp * sinxp * sinsp;
		pm[3,1] =  sinxp;
		pm[3,2] = -sinyp * cosxp;
		pm[3,3] =  cosyp * cosxp;
	end
	pm
end

# --------------------------------------------------------------------------------------------
function gstime(jdut1)
	#  this function finds the greenwich sidereal time (iau-82).
	#
	#  author        : david vallado                  719-573-2600    7 jun 2002
	#
	#  inputs          description                    range / units
	#    jdut1       - julian date of ut1             days from 4713 bc
	#
	#  outputs       :
	#    gst         - greenwich sidereal time        0 to 2pi rad

	tut1 = (jdut1 - 2451545.0) / 36525.0;
	temp = - 6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;

	# 360/86400 = 1/240, to deg, to rad
	gst = rem(temp * pi / 180 / 240, 2*pi )

	# ------------------------ check quadrants --------------------
	(gst < 0) && (gst = gst + 2 * pi)
	gst
end