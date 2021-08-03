struct MTL_short
	band::Int
	rad_mul::Float64
	rad_add::Float64
	rad_max::Float64
	reflect_mul::Float64
	reflect_add::Float64
	reflect_max::Float64
	sun_azim::Float64
	sun_elev::Float64
	sun_dist::Float64
	K1::Float64
	K2::Float64
end

const KW = Dict{Symbol,Any}
find_in_dict = GMT.find_in_dict

"""
    Irgb = truecolor(bndR, bndG, bndB)

Take three Landsat8 UINT16 GMTimages or the file names of those bands and compose
an RGB true color image applying automatic histogram stretching.

Return an UInt8 RGB GMTimage
"""
function truecolor(bndR, bndG, bndB)
	I = isa(bndR, GMT.GMTimage) ? bndR : gmtread(bndR)
	img = Array{UInt8}(undef, size(I,1), size(I,2), 3)
	_ = mat2img(I.image, stretch=true, img8=view(img,:,:,1), scale_only=1)
	I = isa(bndG, GMT.GMTimage) ? bndG : gmtread(bndG)
	@assert size(I,1) == size(img,1) && size(I,2) == size(img,2)
	_ = mat2img(I.image, stretch=true, img8=view(img,:,:,2), scale_only=1)
	I = isa(bndB, GMT.GMTimage) ? bndB : gmtread(bndB)
	@assert size(I,1) == size(img,1) && size(I,2) == size(img,2)
	_ = mat2img(I.image, stretch=true, img8=view(img,:,:,3), scale_only=1)
	Io = mat2img(img, I)
	Io.layout = "TRBa"
	Io
end

# ----------------------------------------------------------------------------------------------------------
"""
read_mtl(band_name::String, mtl::String="")

Use the `band_name` of a Landsat8 band to find the MTL file with the parameters of the scene at which band
belongs and read the params needed to compute Brightness temperature, radiance at top of atmosphere, etc.
If the MTL file does not lieves next to the band file, send its name via the `mtl` argument.

Returns a tuple with:

(band=band, rad_mul=rad_mul, rad_add=rad_add, rad_max=rad_max, reflect_mul=reflect_mul, reflect_add=reflect_add, reflect_max=reflect_max, sun_azim=sun_azim, sun_elev=sun_elev, sun_dis=sun_azim, K1=K1, K2=K2)
"""
function read_mtl(fname::String, mtl::String="")#::NamedTuple{Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}
	_fname = splitext(fname)[1]
	((ind = findfirst("_B", fname)) === nothing) && error("This $(fname) is not a valid Landsat8 band file")
	if (mtl == "")
		mtl = fname[1:ind[1]] * "MTL.txt"
		if (!isfile(mtl))
			pato = splitdir(_fname)[1]
			lst = filter(x -> endswith(x, "MTL.txt"), readdir((pato == "") ? "." : pato))
			if (length(lst) == 1 && startswith(lst[1], _fname[1:16]))
				mtl = joinpath(pato, lst[1])
				!isfile(mtl) && error("MTL file was not transmitted in input and I couldn't find it next the band file.")
			end
		end
	end
	band = parse(Int, _fname[ind[1]+2:end])

	f = open(mtl);	lines = readlines(f);	close(f)

	function get_par(str)
		ind = findfirst(findfirst.(str, lines) .!== nothing)
		parse(Float64, split(lines[ind], "=")[2])
	end

	rad_mul = get_par("RADIANCE_MULT_BAND_$(band)")
	rad_add = get_par("RADIANCE_ADD_BAND_$(band)")
	rad_max = get_par("RADIANCE_MAXIMUM_BAND_$(band)")
	reflect_mul = (band < 10) ? get_par("REFLECTANCE_MULT_BAND_$(band)") : 1.0
	reflect_add = (band < 10) ? get_par("REFLECTANCE_ADD_BAND_$(band)") : 0.0
	reflect_max = (band < 10) ? get_par("REFLECTANCE_MAXIMUM_BAND_$(band)") : 0.0
	sun_azim = get_par("SUN_AZ")
	sun_elev = get_par("SUN_EL")
	sun_dist = get_par("EARTH_SUN")
	K1 = (band >= 10) ? get_par("K1_CONSTANT_BAND_$(band)") : 0.0
	K2 = (band >= 10) ? get_par("K2_CONSTANT_BAND_$(band)") : 0.0
	#(band=band, rad_mul=rad_mul, rad_add=rad_add, rad_max=rad_max, reflect_mul=reflect_mul, reflect_add=reflect_add, reflect_max=reflect_max, sun_azim=sun_azim, sun_elev=sun_elev, sun_dist=sun_azim, K1=K1, K2=K2)
	MTL_short(band, rad_mul, rad_add, rad_max, reflect_mul, reflect_add, reflect_max, sun_azim, sun_elev, sun_dist, K1, K2)
end

# ----------------------------------------------------------------------------------------------------------
function helper1_sats(fname::String)
	pars = read_mtl(fname)
	I::GMTimage{UInt16, 2} = gmtread(fname)
	indNaN = fill(false, size(I))
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)	# 5x faster than: indNaN = (I.image .== 0)
		(I.image[k] == 0) && (indNaN[k] = true)
	end
	o = Matrix{Float32}(undef, size(I))
	return I, pars, indNaN, o
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = bright_T(fname::String)

Returns a GMTgrid with the brigthness temperature of Landasat8 termal band (10 or 11)
"""
function bright_T(fname::String)
	I, pars, indNaN, o = helper1_sats(fname)
	(pars.band < 10) && error("Brightness temperature is only for bands 10 or 11. Not this: $(pars.band)")
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
		o[k] = pars.K2 / (log(pars.K1 / (I.image[k] * pars.rad_mul + pars.rad_add) + 1.0)) - 273.15
	end
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)  indNaN[k] && (o[k] = NaN)  end
	mat2grid(o, I)
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = radiance_TOA(fname::String)

Returns a GMTgrid with the radiance at TopOfAtmosphere for the Landsat8 band file `fname`
"""
function radiance_TOA(fname::String)
	I, pars, indNaN, o = helper1_sats(fname)
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
		o[k] = I.image[k] * pars.rad_mul + pars.rad_add
	end
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)  indNaN[k] && (o[k] = NaN)  end
	mat2grid(o, I)
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = reflectance_TOA(fname::String)

Returns a GMTgrid with the TopOfAtmosphere planetary reflectance for the Landsat8 band file `fname`
"""
function reflectance_TOA(fname::String)
	I, pars, indNaN, o = helper1_sats(fname)
	(pars.band >= 10) && error("Computing Reflectance for Thermal bands is not defined.")
	s_elev = sin(pars.sun_elev * pi/180)
	fact_x = pars.reflect_mul / s_elev
	fact_a = pars.reflect_add / s_elev
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
		o[k] = I.image[k] * fact_x + fact_a
	end
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)  indNaN[k] && (o[k] = NaN)  end
	mat2grid(o, I)
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = reflectance_surf(fname::String)

Compute the radiance-at-surface of Landsat8 band using the COST model.

Returns a Float32 GMTgrid type
"""
function reflectance_surf(fname::String)
	I, pars, indNaN, o = helper1_sats(fname)
	(pars.band >= 10) && error("Computing Surface Reflectance for Thermal bands is not defined.")

	s_elev = sin(pars.sun_elev * pi/180)
	Esun = (pi * pars.sun_dist ^2) * pars.rad_max / pars.reflect_max
	TAUv = 1.0;		TAUz = s_elev;		Esky = 0.0;		sun_prct=1;
	(pars.band == 6 || pars.band == 7 || pars.band == 9) && (TAUz = 1.0)
	Sun_Radiance = TAUv * (Esun * s_elev * TAUz + Esky) / (pi * pars.sun_dist ^2)

	tmp = I.image
	tmp = filter(tmp -> tmp != 0, tmp)		# This is 3 times faster then: tmp = I.image[I.image .> 0]
	darkDN = quantile!(tmp, 0.01)
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
		o[k] = I.image[k] * pars.rad_mul + pars.rad_add
	end
	radiance_dark = pars.rad_mul * darkDN + pars.rad_add	# 0.01%
	LHaze = radiance_dark - sun_prct * Sun_Radiance / 100
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
		o[k] = (o[k] - LHaze) / Sun_Radiance
	end

	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
		if     (o[k] < 0) o[k] = 0
		elseif (o[k] > 1) o[k] = 1
		end
		indNaN[k] && (o[k] = NaN)
	end
	mat2grid(o, I)
end

# ----------------------------------------------------------------------------------------------------------
"""
    NDVI = ndvi(bndR, bndNIR; threshold=0.4, mask::Bool=false)

Compute the NDVI vegetation index. Input can be either the bands file names, or GMTimage objects
with the band's data.

Returns either a Float32 GMTgrid or a UInt8 GMTimage if the `mask` option is set to true.
"""
function ndvi(bndR::String, bndNIR::String; radiance::Bool=false, threshold=0.4, mask::Bool=false)
	Ired = (radiance) ? radiance_TOA(bndR)   : gmtread(bndR)
	Inir = (radiance) ? radiance_TOA(bndNIR) : gmtread(bndNIR)
	ndvi(Ired, Inir; mask=mask, threshold=threshold)
end
function ndvi(Ired, Inir; threshold=0.4, mask::Bool=false)
	@assert size(Ired) == size(Inir)
	if (mask)
		img = fill(UInt8(0), size(Ired))
		@inbounds Threads.@threads for k = 1:size(Ired,1)*size(Ired,2)
			t = (Inir[k] - Ired[k]) / (Inir[k] + Ired[k])
			(t >= threshold && t <= 1) && (img[k] = 255)
		end
		I = mat2img(img, proj4=Ired.proj4, wkt=Ired.wkt, x=Ired.x, y=Ired.y)
		I.range[5], I.range[6] = 0, 255;	I.epsg = Ired.epsg;		I.layout = "BRPa"
		I
	else
		img = fill(NaN32, size(Ired))
		@inbounds Threads.@threads for k = 1:size(Ired,1)*size(Ired,2)
			t = (Inir[k] - Ired[k]) / (Inir[k] + Ired[k])
			(t >= threshold && t <= 1) && (img[k] = t)
		end
		G = mat2grid(img, Ired)
	end
end
