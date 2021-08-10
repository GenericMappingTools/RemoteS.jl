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
	Io = mat2img(img, I);	Io.layout = "TRBa"
	Io
end
function truecolor(cube::GMT.GMTimage{UInt16, 3}, wavelength; kw...)
	img = Array{UInt8}(undef, size(cube,1), size(cube,2), 3)
	if (wavelength == "Sentinel")
		_ = mat2img(@view(cube[:,:,4]), stretch=true, img8=view(img,:,:,1), scale_only=1)
		_ = mat2img(@view(cube[:,:,3]), stretch=true, img8=view(img,:,:,2), scale_only=1)
		_ = mat2img(@view(cube[:,:,2]), stretch=true, img8=view(img,:,:,3), scale_only=1)
		Io = mat2img(img, cube);	Io.layout = "TRBa"
		Io
	end
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
	I::GMT.GMTimage{UInt16, 2} = gmtread(fname)
	indNaN = fill(false, size(I))
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)	# 5x faster than: indNaN = (I.image .== 0)
		(I.image[k] == 0) && (indNaN[k] = true)
	end
	o = Matrix{Float32}(undef, size(I))
	return I, pars, indNaN, o
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = dn2temperature(fname::String)

Returns a GMTgrid with the brigthness temperature of Landasat8 termal band (10 or 11)
"""
function dn2temperature(fname::String)
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
    R = dn2radiance(fname::String)

Returns a GMTgrid with the radiance at TopOfAtmosphere for the Landsat8 band file `fname`
"""
function dn2radiance(fname::String)
	I, pars, indNaN, o = helper1_sats(fname)
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
		o[k] = I.image[k] * pars.rad_mul + pars.rad_add
	end
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)  indNaN[k] && (o[k] = NaN)  end
	mat2grid(o, I)
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = dn2reflectance(fname::String)

Returns a GMTgrid with the TopOfAtmosphere planetary reflectance for the Landsat8 band file `fname`
"""
function dn2reflectance(fname::String)
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

#=
function ndvi(bndR::String, bndNIR::String; radiance::Bool=false, threshold=0.4, mask::Bool=false)
	Ired = (radiance) ? dn2radiance(bndR)   : gmtread(bndR)
	Inir = (radiance) ? dn2radiance(bndNIR) : gmtread(bndNIR)
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
=#

# ----------------------------------------------------------------------------------------------------------
"""
    CLG = clg(green, redEdge3; kw...)

Green cholorphyl index. Wu et al 2012.

CLG = (redEdge3)/(green)-1 
"""
clg(green, redEdge3; kw...) = spectral_indices(green, redEdge3; index="CLG", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    CLRE = clre(redEdge1, redEdge3; kw...)

RedEdge cholorphyl index. Clevers and Gitelson 2013.

CLRE = (redEdge3)/(redEdge1)-1
"""
clre(redEdge1, redEdge3; kw...) = spectral_indices(redEdge1, redEdge3; index="CLRE", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    EVI = evi(blue, red, nir; kw...)

Enhanced vegetation index. Huete et al 1990

EVI = G * ((nir - red) / (nir + C1 * red - C2 * blue + Levi));
C1, C2, G, Levi = 6.0, 7.5, 2.5, 1.
"""
evi(blue, red, nir; kw...) = spectral_indices(blue, red, nir; index="EVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    EVI2 = evi(red, nir; kw...)

Two-band Enhanced vegetation index. Jiang et al 2008

EVI2 = G * ((nir - red) / (nir + 2.4 * red ))
"""
evi2(red, nir; kw...) = spectral_indices(red, nir; index="EVI2", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    GNDVI = gndvi(green, nir; kw...)

green Normalized diff vegetation index: more sensitive to cholorphyll than ndvi. Gitelson, A., and M. Merzlyak

GNDVI = (nir - green) / (nir + green)
"""
gndvi(green, nir; kw...) = spectral_indices(green, nir; index="GNDVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MNDWI = mndwi(green, swir2; kw...)

Modified Normalised Difference Water Index. Xu2006

MNDWI = (green-swir2) / (green+swir2)
"""
mndwi(green, swir2; kw...) = spectral_indices(swir2, green; index="MNDWI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MTCI = mtci(red, redEdge1, redEdge2; kw...)

Meris Terrestrial Chlorophyll Index. Clevers and Gitelson 2013, Dash and Curran 2004

MTCI = (redEdge2-redEdge1) / (redEdge1-red)
"""
mtci(red, redEdge1, redEdge2; kw...) = spectral_indices(red, redEdge1, redEdge2; index="MTCI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MCARI = mcari(green, red, redEdge1; kw...)

Modified Chlorophyll Absorption ratio index. Daughtery et al. 2000

MCARI = (redEdge1 - red - 0.2 * (redEdge1 + green)) * (redEdge1 / red)
"""
mcari(green, red, redEdge1; kw...) = spectral_indices(green, red, redEdge1; index="MCARI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MSAVI = msavi(red, nir; kw...)

Modified soil adjusted vegetation index. Qi 1994

MSAVI = nir + 0.5 - (0.5 * sqrt(pow(2.0 * nir + 1.0, 2) - 8.0 * (nir - (2.0 * red))))
"""
msavi(red, nir; kw...) = spectral_indices(red, nir; index="MSAVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NBRI = nbri(nir, swir3; kw...)

Normalised Burn Ratio Index. Garcia 1991

NBRI = (nir - swir3) / (nir + swir3)
"""
nbri(nir, swir3; kw...) = spectral_indices(swir3, nir; index="NBRI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDVI = ndvi(red, nir; kw...)

Compute the NDVI vegetation index. Input can be either the bands file names, or GMTimage objects
with the band's data.

NDVI = (nir - red) / (nir + red)

Returns either a Float32 GMTgrid or a UInt8 GMTimage if the `mask` option is set to true.
"""
ndvi(red, nir; kw...) = spectral_indices(red, nir; index="NDVI", kw...)
function ndvi(cube::GMTimage{UInt16, 3}, wavelength; kw...)
	if (wavelength == "Sentinel")
		o = spectral_indices(@view(cube[:,:,4]), @view(cube[:,:,8]); index="NDVI", kw...)
		return isa(o, Float32) ? mat2grid(o, cube) : mat2img(o, cube)
	end
end

# ----------------------------------------------------------------------------------------------------------
"""
    NDWI = ndwi(green, nir; kw...)

Normalized difference water index. McFeeters 1996. NDWI => (green - nir)/(green + nir)

NDWI = (green - nir)/(green + nir)
"""
ndwi(green, nir; kw...) = spectral_indices(nir, green; index="NDWI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDWI2 = ndwi2(nir, swir2; kw...)

Normalized difference water index. Gao 1996, Chen 2005 (also known as Normalized Difference Moisture Index
NDBI and LSWI)

NDWI2 = (nir - swir2)/(nir + swir2)
"""
ndwi2(nir, swir2; kw...) = spectral_indices(swir2, nir; index="NDWI2", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDREI1 = ndrei1(redEdge1, redEdge2; kw...)

Normalized difference red edge index. Gitelson and Merzlyak 1994

NDREI1 = (redEdge2 - redEdge1) / (redEdge2 + redEdge1)
"""
ndrei1(redEdge1, redEdge2; kw...) = spectral_indices(redEdge1, redEdge2; index="NDREI1", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDREI2 = ndrei2(redEdge1, redEdge3; kw...)

Normalized difference red edge index 2. Barnes et al 2000

NDREI2 = (redEdge3 - redEdge1) / (redEdge3 + redEdge1)
"""
ndrei2(redEdge1, redEdge3; kw...) = spectral_indices(redEdge1, redEdge3; index="NDREI2", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    SATVI = satvi(red, swir2, swir3; kw...)

Soil adjusted total vegetation index. Marsett 2006

SATVI = ((swir2 - red) / (swir2 + red + L)) * (1.0 + L) - (swir3 / 2.0)
"""
satvi(red, swir2, swir3; kw...) = spectral_indices(red, swir2, swir3; index="SATVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    SAVI = savi(red, nir; kw...)

Soil adjusted vegetation index. Huete 1988

SAVI = (nir - red) * (1.0 + L) / (nir + red + L)
"""
savi(red, nir; kw...) = spectral_indices(red, nir; index="SAVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    SLAVI = slavi(red, nir, swir2; kw...)

Specific Leaf Area Vegetation Index. Lymburger 2000

SLAVI = nir / (red + swir2)
"""
slavi(red, nir, swir2; kw...) = spectral_indices(red, nir, swir2; index="SLAVI", kw...)
function slavi(cube::GMT.GMTimage{UInt16, 3}, wavelength; kw...)
	if (wavelength == "Sentinel")	
		o = spectral_indices(@view(cube[:,:,4]), @view(cube[:,:,8]), @view(cube[:,:,11]); index="SLAVI", kw...)
		return isa(o, Float32) ? mat2grid(o, cube) : mat2img(o, cube)
	end
end

# ----------------------------------------------------------------------------------------------------------
function spectral_indices(bnd1::String, bnd2::String, bnd3::String=""; index::String="", kwargs...)
	do_radTOA = any(keys(kwargs) .== :dn2radiance)
	do_refTOA = any(keys(kwargs) .== :dn2reflectance)
	do_refSrf = any(keys(kwargs) .== :reflectance_surf)
	Bnd1 = (do_radTOA) ? dn2radiance(bnd1) : (do_refTOA) ? dn2reflectance(bnd1) : (do_refSrf) ? reflectance_surf(bnd1) : gmtread(bnd1)
	Bnd2 = (do_radTOA) ? dn2radiance(bnd2) : (do_refTOA) ? dn2reflectance(bnd2) : (do_refSrf) ? reflectance_surf(bnd2) : gmtread(bnd2)
	if (bnd3 != "")
		Bnd3 = (do_radTOA) ? dn2radiance(bnd3) : (do_refTOA) ? dn2reflectance(bnd3) : (do_refSrf) ? reflectance_surf(bnd3) : gmtread(bnd3)
	else
		Bnd3 = nothing
	end
	spectral_indices(Bnd1, Bnd2, Bnd3; index=index, kwargs...)
end

function spectral_indices(bnd1, bnd2, bnd3=nothing; index::String="", kwargs...)
	(index == "") && error("Must select which index to compute")
	@assert size(bnd1) == size(bnd2)
	mask = any(keys(kwargs) .== :mask)

	d = KW(kwargs)
	mask = (haskey(d, :mask))
	threshold = find_in_dict(d, [:threshold])[1]
	classes = (threshold === nothing) ? find_in_dict(d, [:classes])[1] : nothing
	(mask && threshold === nothing) && error("The `mask` option requires the `threshold=x` option")
	(classes !== nothing && length(classes) > 3) && (classes = classes[1:3]; @warn("`classes` maximum elements is 3. Clipping the others"))

	img = (mask || classes !== nothing) ? fill(UInt8(0), size(bnd1)) : fill(NaN32, size(bnd1))

	mn = size(bnd1,1) * size(bnd1,2)
	C1, C2, G, L, Levi = 6.0, 7.5, 2.5, 0.5, 1.0
	if (index == "CLG" || index == "CLRE")
		# Green cholorphyl index. Wu et al 2012		CLG               (redEdge3)/(green)-1
		# RedEdge cholorphyl index. Clevers and Gitelson 2013	CLRE  (redEdge3)/(redEdge1)-1 
		if (mask)
			@inbounds Threads.@threads for k = 1:mn
				t = bnd2[k] / bnd1[k] - 1
				(t >= threshold) && (img[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = bnd2[k] / bnd1[k] - 1
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "EVI")			# Enhanced vegetation index. Huete et al 1990
		# G * ((nir - red) / (nir + C1 * red - C2 * blue + Levi))
		if (mask)
			@inbounds Threads.@threads for k = 1:mn
				t = G * (bnd3[k] - bnd2[k]) / (bnd3[k] + C1 * bnd2[k] - C2 * bnd1[k] + Levi)
				(t >= threshold) && (img[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = G * (bnd3[k] - bnd2[k]) / (bnd3[k] + C1 * bnd2[k] - C2 * bnd1[k] + Levi)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "EVI2")		# Two-band Enhanced vegetation index. Jiang et al 2008 
		# G * ((nir - red) / (nir + 2.4 * red ))
		if (mask)
			@inbounds Threads.@threads for k = 1:mn
				t = (bnd2[k] - bnd1[k]) / (bnd1[k] + 2.4 * bnd2[k])
				(t >= threshold) && (img[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t = G * (bnd2[k] - bnd1[k]) / (bnd2[k] + 2.4 * bnd1[k])
				(t >= -1 && t <= 1) && (img[k] = t)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "MTCI")		# Meris Terrestrial Chlorophyll Index. Clevers and Gitelson 2013, Dash and Curran 2004 
		# (redEdge2-redEdge1) / (redEdge1-red)
		if (mask)
			@inbounds Threads.@threads for k = 1:mn
				t = (bnd3[k] - bnd2[k]) / (bnd2[k] - bnd1[k])
				(t >= threshold) && (img[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = (bnd3[k] - bnd2[k]) / (bnd2[k] - bnd1[k])
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "MCARI")		# Modified Chlorophyll Absorption ratio index. Daughtery et al. 2000 
		# (redEdge1 - red - 0.2 * (redEdge1 + green)) * (redEdge1 / red)
		if (mask)
			@inbounds Threads.@threads for k = 1:mn
				t = (bnd3[k] - bnd2[k] - 0.2 * (bnd3[k] - bnd1[k])) * (bnd3[k] / bnd2[k])
				(t >= threshold) && (img[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = (bnd3[k] - bnd2[k] - 0.2 * (bnd3[k] - bnd1[k])) * (bnd3[k] / bnd2[k])
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "MSAVI")		# Modified soil adjusted vegetation index.
		# nir + 0.5 - (0.5 * sqrt(pow(2.0 * nir + 1.0, 2) - 8.0 * (nir - (2.0 * red))))
		if (mask)
			@inbounds Threads.@threads for k = 1:mn
				t = bnd2[k] + 0.5 - (0.5 * sqrt((2 * bnd2[k] + 1) ^2) - 8 * (bnd2[k] - (2 * bnd1[k])))
				(t >= threshold) && (img[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = bnd2[k] + 0.5 - (0.5 * sqrt((2 * bnd2[k] + 1) ^2) - 8 * (bnd2[k] - (2 * bnd1[k])))
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "GNDVI" || index == "MNDWI" || index == "NBRI" || index == "NDVI" || index == "NDWI" || index == "NDWI2" || index == "NDREI1" || index == "NDREI2")
		# green Normalized diff vegetation index: more sensitive to cholorphyll than ndvi. Gitelson, A., and M. Merzlyak. 
		# GNDVI		=> (nir - green) / (nir + green)
		# Modified Normalised Difference Water Index. Xu2006;  MNDWI	=> (green-swir2) / (green+swir2)
		# Normalised Burn Ratio Index. NBRI => (nir - swir3) / (nir + swir3)
		# Normalized difference water index. McFeeters 1996. NDWI => (green - nir)/(green + nir)
		# NDVI Normalised Difference Vegetation Index. Rouse 1974. (nir - red) / (nir + red)
		# NDBI, LSWI. Normalized difference water index. Gao 1996, Chen 2005; NDWI2 => (nir - swir2)/(nir + swir2)
		# Normalized difference red edge index. Gitelson and Merzlyak 1994; (redEdge2 - redEdge1)/(redEdge2 + redEdge1)
		# Normalized difference red edge index 2. Barnes et al 2000; (redEdge3 - redEdge1)/(redEdge3 + redEdge1)
		if (mask)
			@inbounds Threads.@threads for k = 1:mn
				t = (bnd2[k] - bnd1[k]) / (bnd1[k] + bnd2[k])
				(t >= threshold && t <= 1) && (img[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t = (bnd2[k] - bnd1[k]) / (bnd1[k] + bnd2[k])
				(t >= -1 && t <= 1) && (img[k] = t)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "SATVI")		# Soil adjusted total vegetation index.
		# ((swir2 - red) / (swir2 + red + L)) * (1.0 + L) - (swir3 / 2.0)
		if (mask)
			@inbounds Threads.@threads for k = 1:mn
				t = ((bnd2[k] - bnd1[k]) / (bnd2[k] + bnd1[k] + L)) * (1.0 + L) - (bnd3[k] / 2.0) 
				(t >= threshold) && (img[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = ((bnd2[k] - bnd1[k]) / (bnd2[k] + bnd1[k] + L)) * (1.0 + L) - (bnd3[k] / 2.0) 
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "SAVI")		# Soil adjusted vegetation index. Huete1988
		# (nir - red) * (1.0 + L) / (nir + red + L);
		if (mask)
			@inbounds Threads.@threads for k = 1:mn
				t = (bnd2[k] - bnd1[k]) * (1.0 + L) / (bnd2[k] - bnd1[k] + L) 
				(t >= threshold) && (img[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = (bnd2[k] - bnd1[k]) * (1.0 + L) / (bnd2[k] - bnd1[k] + L) 
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "SLAVI")		# Soil Adjusted Vegetation Index Huete 1988.
		# nir / (red + swir2)
		if (mask)
			@inbounds Threads.@threads for k = 1:mn
				t = bnd2[k] / (bnd1[k] + bnd3[k]) 
				(t >= threshold) && (img[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = bnd2[k] / (bnd1[k] + bnd3[k]) 
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	end
	if (isa(bnd1, GMTimage))
		if (mask)
			I = mat2img(img, proj4=bnd1.proj4, wkt=bnd1.wkt, x=bnd1.x, y=bnd1.y)
			I.range[5], I.range[6] = 0, 255;	I.epsg = bnd1.epsg;		I.layout = "BRPa"
		end
		#return (mask || classes !== nothing) ? I : mat2grid(img, bnd1)
		return (mask || classes !== nothing) ? helper_si2(img, bnd1, mask, classes) : mat2grid(img, bnd1)
	else
		return img
	end
end

function helper_si2(img, refimg, mask, classes)
	I = mat2img(img, proj4=refimg.proj4, wkt=refimg.wkt, x=refimg.x, y=refimg.y)
	I.epsg = refimg.epsg;		I.layout = "BRPa"
	I.range[5], I.range[6] = 0, (mask) ? 255 : length(classes)
	return I
end

function helper_si!(img, threshold, classes)
	# Helper function to apply a threshold to a grid OR create a classes image
	mn = size(img, 1) * size(img, 2)
	if (threshold !== nothing)		# In this case we return one GMTgrid where < thres = NaN
		@inbounds Threads.@threads for k = 1:mn
			(img[k] < threshold) && (img[k] = NaN32)
		end
	elseif (classes !== nothing)	# Here we return a GMTimage with up to 4 classes (maximum allowed)
		@inbounds Threads.@threads for k = 1:mn
			(img[k] >= classes[1]) && (img[k] = 1)
			(img[k] >= classes[2]) && (img[k] = 2)
		end
		if (length(classes) == 3)
			@inbounds Threads.@threads for k = 1:mn 
				(img[k] >= classes[3]) && (img[k] = 3)
			end
		end
	end
end
