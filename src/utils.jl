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

Take three Landsat8/Sentinel2 UINT16 GMTimages or the file names of those bands and compose
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

# ----------------------------------------------------------------------------------------------------------
"""
    CLG = vi_clg(green, redEdge3; kw...)

Green cholorphyl index. Wu et al 2012.

CLG = (redEdge3)/(green)-1 
"""
vi_clg(green, redEdge3; kw...) = spectral_indices(green, redEdge3; index="CLG", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    CLRE = vi_clre(redEdge1, redEdge3; kw...)

RedEdge cholorphyl index. Clevers and Gitelson 2013.

CLRE = (redEdge3)/(redEdge1)-1
"""
vi_clre(redEdge1, redEdge3; kw...) = spectral_indices(redEdge1, redEdge3; index="CLRE", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    EVI = vi_evi(blue, red, nir; kw...)

Enhanced vegetation index. Huete et al 1990

EVI = G * ((nir - red) / (nir + C1 * red - C2 * blue + Levi));
C1, C2, G, Levi = 6.0, 7.5, 2.5, 1.
"""
vi_evi(blue, red, nir; kw...) = spectral_indices(blue, red, nir; index="EVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    EVI2 = vi_evi(red, nir; kw...)

Two-band Enhanced vegetation index. Jiang et al 2008

EVI2 = G * ((nir - red) / (nir + 2.4 * red ))
"""
vi_evi2(red, nir; kw...) = spectral_indices(red, nir; index="EVI2", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    GNDVI = vi_gndvi(green, nir; kw...)

green Normalized diff vegetation index: more sensitive to cholorphyll than ndvi. Gitelson, A., and M. Merzlyak

GNDVI = (nir - green) / (nir + green)
"""
vi_gndvi(green, nir; kw...) = spectral_indices(green, nir; index="GNDVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MNDWI = vi_mndwi(green, swir2; kw...)

Modified Normalised Difference Water Index. Xu2006

MNDWI = (green-swir2) / (green+swir2)
"""
vi_mndwi(green, swir2; kw...) = spectral_indices(green, swir2; index="MNDWI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MTCI = vi_mtci(red, redEdge1, redEdge2; kw...)

Meris Terrestrial Chlorophyll Index. Clevers and Gitelson 2013, Dash and Curran 2004

MTCI = (redEdge2-redEdge1) / (redEdge1-red)
"""
vi_mtci(red, redEdge1, redEdge2; kw...) = spectral_indices(red, redEdge1, redEdge2; index="MTCI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
"""
function spectral_indices(bnd1::String, bnd2::String, bnd3::String=""; index::String="", mask::Bool=false, kwargs...)
	do_radTOA = any(keys(kwargs) .== :radiance_TOA)
	do_refTOA = any(keys(kwargs) .== :reflectance_TOA)
	do_refSrf = any(keys(kwargs) .== :reflectance_surf)
	Bnd1 = (do_radTOA) ? radiance_TOA(bnd1) : (do_refTOA) ? reflectance_TOA(bnd1) : (do_refSrf) ? reflectance_surf(bnd1) : gmtread(bnd1)
	Bnd2 = (do_radTOA) ? radiance_TOA(bnd2) : (do_refTOA) ? reflectance_TOA(bnd2) : (do_refSrf) ? reflectance_surf(bnd2) : gmtread(bnd2)
	if (bnd3 != "")
		Bnd3 = (do_radTOA) ? radiance_TOA(bnd3) : (do_refTOA) ? reflectance_TOA(bnd3) : (do_refSrf) ? reflectance_surf(bnd3) : gmtread(bnd3)
	else
		Bnd3 = nothing
	end
	spectral_indices(bnd1, bnd2, bnd3; index=index, mask=mask, kwargs...)
end
function spectral_indices(bnd1, bnd2, bnd3=nothing; index::String="", mask::Bool=false, kwargs...)
	(index == "") && error("Must select which index to compute")
	@assert size(bnd1) == size(bnd2)
	img = (mask) ? fill(UInt8(0), size(bnd1)) : fill(NaN32, size(bnd1))
	mn = size(bnd1,1)*size(bnd1,2)
	C1, C2, G, L, Levi = 6.0, 7.5, 2.5, 0.5, 1.0
	if (index == "CLG" || index == "CLRE")
		# Green cholorphyl index. Wu et al 2012		CLG               (redEdge3)/(green)-1
		# RedEdge cholorphyl index. Clevers and Gitelson 2013	CLRE  (redEdge3)/(redEdge1)-1 
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = bnd2[k] / bnd1[k] - 1
			end
		end
	elseif (index == "EVI")			# Enhanced vegetation index. Huete et al 1990
		# G * ((nir - red) / (nir + C1 * red - C2 * blue + Levi))
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = G * (bnd3[k] - bnd2[k]) / (bnd3[k] + C1 * bnd2[k] - C2 * bnd1[k] + Levi)
			end
		end
	elseif (index == "EVI2")		# Two-band Enhanced vegetation index. Jiang et al 2008 
		# G * ((nir - red) / (nir + 2.4 * red ))
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				t = G * (bnd2[k] - bnd1[k]) / (bnd2[k] + 2.4 * bnd1[k])
				(t >= -1 && t <= 1) && (img[k] = t)
			end
		end
	elseif (index == "GNDVI" || index == "MNDWI")
		# green Normalized diff vegetation index: more sensitive to cholorphyll than ndvi. Gitelson, A., and M. Merzlyak. 
		# GNDVI		=> (nir - green)/( nir + green)
		# Modified Normalised Difference Water Index. Xu2006;  MNDWI	=> (green-swir2) / (green+swir2)
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				t = (bnd2[k] - bnd1[k]) / (bnd2[k] + bnd1[k])
				(t >= -1 && t <= 1) && (img[k] = t)
			end
		end
	elseif (index == "MTCI")		# Meris Terrestrial Chlorophyll Index. Clevers and Gitelson 2013, Dash and Curran 2004 
		# (redEdge2-redEdge1) / (redEdge1-red)
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = (bnd3[k] - bnd2[k]) / (bnd2[k] - bnd1[k])
			end
		end
	elseif (index == "MCARI")		# Modified Chlorophyll Absorption ratio index. Daughtery et al. 2000 
		# (redEdge1 - red - 0.2 * (redEdge1 + green)) * (redEdge1 / red)
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = (bnd3[k] - bnd2[k] - 0.2 * (bnd3[k] - bnd1[k])) * (bnd3[k] / bnd2[k])
			end
		end
	elseif (index == "MSAVI")		# Modified soil adjusted vegetation index.
		# nir + 0.5 - (0.5 * sqrt(pow(2.0 * nir + 1.0, 2) - 8.0 * (nir - (2.0 * red))))
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = bnd2[k] + 0.5 - (0.5 * sqrt((2 * bnd2[k] + 1) ^2) - 8 * (bnd2[k] - (2 * bnd1[k])))
			end
		end
	#elseif (index == "NDVIC")		# Normalized difference vegetation index.
	elseif (index == "NBRI" || index == "NDWI" || index == "NDWI2" || index == "NDREI1" || index == "NDREI2")
		# Normalised Burn Ratio Index. NRBI => (nir - swir3) / (nir + swir3)
		# Normalized difference water index. McFeeters 1996. NDWI => (green - nir)/(green + nir)
		# NDBI, LSWI. Normalized difference water index. Gao 1996, Chen 2005; NDWI2 => (nir - swir2)/(nir + swir2)
		# Normalized difference red edge index. Gitelson and Merzlyak 1994; (redEdge2 - redEdge1)/(redEdge2 + redEdge1)
		# Normalized difference red edge index 2. Barnes et al 2000; (redEdge3 - redEdge1)/(redEdge3 + redEdge1)
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				t = (bnd1[k] - bnd2[k]) / (bnd1[k] + bnd2[k])
				(t >= -1 && t <= 1) && (img[k] = t)
			end
		end
	elseif (index == "SATVI")		# Soil adjusted total vegetation index.
		# ((swir2 - red) / (swir2 + red + L)) * (1.0 + L) - (swir3 / 2.0)
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = ((bnd2[k] - bnd1[k]) / (bnd2[k] + bnd1[k] + L)) * (1.0 + L) - (bnd3[k] / 2.0) 
			end
		end
	elseif (index == "SAVI")		# Soil adjusted vegetation index. Huete1988
		# (nir - red) * (1.0 + L) / (nir + red + L);
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = (bnd2[k] - bnd1[k]) * (1.0 + L) / (bnd2[k] - bnd1[k] + L) 
			end
		end
	elseif (index == "SLAVI")		# 
		# nir / (red + swir2)
		if (mask)
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = bnd2[k] / (bnd1[k] + bnd3[k]) 
			end
		end
	end
	return (mask) ? nothing : mat2grid(img, bnd1)
end