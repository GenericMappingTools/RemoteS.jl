const generic_docs = "
- The first form, `evi(blue, red, nir; kw...)`, accepts inputs as matrices, or file names of the data bands.
- The second form is more versatile but also more complex to describe.
  - `cube`: Is the file name of a 'cube', a multi-layered file normally created with the [`cutcube`](@ref) function.
     If this file was created with band descriptions one can use the `bands` or the `bandnames` options.
  - `bands`: _cubes_ created with [`cutcube`](@ref) assign descriptions starting with \"Band 1 ...\" an so on forth
    the other bands. So when `bands` is used we search for bands named \"Band band[k]\", where band[k] loops
    over all elements of the `bands` vector. WARNING: the elements order in the vector must be sorted in increasing
    wavelength numbers, _i.e._ like the example for the first form.
  - `layers`: Use this option when you are certain of the bands order in the cube or the it doesn't have a bands
    description. The selection will be made with cube[:,:,layer[1]], etc... WARNING: same warn as above.
  - `bandnames`: When we know the common designation of a band, for example \"Green\", or any part of a band
    description, for example \"NIR\", we can use that info to create a `bandnames` string vector that will be
    matched against the cube's bands descriptions.
- `kwargs`:
  - `threshold`: When a threshold is provided we return a GMTgrid where `vals[ij] < threshold = NaN`
  - `classes`: is a vector with up to 3 elements (class separators) and we return a  UInt8 GMTimage with the
    indices categorized into vals[ij] > classes[1] = 1; vals[ij] > classes[2] = 2; vals[ij] > classes[3] = 3 and 0 otherwise.
  - `mask`: Used together with `threshold` outputs a UInt8 GMTimage mask with `vals[ij] > threshold = 255` and 0 otherwise

If none of `bands`, `layers` or `bandnames` is provided, we use the default band names shown in the first form.

Returns either a Float32 GMTgrid or a UInt8 GMTimage if the `mask` or `classes` options are used.
"

# ----------------------------------------------------------------------------------------------------------
function helper_si_method(cube::String, index::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[],
	                      bandnames::Vector{String}=String[], defbandnames::Vector{String}=String[], kw...)
	# Helper function to compute Spectral Indices from a 'cube' file and somehow band selection.
	(isempty(bands) && isempty(bandnames) && isempty(layers)) && (bandnames = defbandnames)
	sp_indices(subcube(cube, bands=bands, layers=layers, bandnames=bandnames), [1,2]; index=index, kw...)
end

# ----------------------------------------------------------------------------------------------------------
"""
    CLG = clg(green, redEdge3; kw...)

Green cholorphyl index. Wu et al 2012.

CLG = (redEdge3)/(green)-1 
"""
clg(green, redEdge3; kw...) = sp_indices(green, redEdge3; index="CLG", kw...)
clg(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="CLG", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    CLRE = clre(redEdge1, redEdge3; kw...)

RedEdge cholorphyl index. Clevers and Gitelson 2013.

CLRE = (redEdge3)/(redEdge1)-1
"""
clre(redEdge1, redEdge3; kw...) = sp_indices(redEdge1, redEdge3; index="CLRE", kw...)
clre(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="CLRE", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    EVI = evi(blue, red, nir; kw...)
or

    EVI = evi(cube::String; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Enhanced vegetation index. Huete et al 1990

EVI = G * ((nir - red) / (nir + C1 * red - C2 * blue + Levi));
C1, C2, G, Levi = 6.0, 7.5, 2.5, 1.

$(generic_docs)

"""
evi(blue, red, nir; kw...) = sp_indices(blue, red, nir; index="EVI", kw...)
evi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 3); index="EVI", kw...)
evi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "EVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["blue", "red", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    EVI2 = evi2(red, nir; kw...)
or

    EVI2 = evi2(cube::String; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Two-band Enhanced vegetation index. Jiang et al 2008

EVI2 = G * ((nir - red) / (nir + 2.4 * red ))

$(generic_docs)
"""
evi2(red, nir; kw...) = sp_indices(red, nir; index="EVI2", kw...)
evi2(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="EVI2", kw...)
evi2(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "EVI2"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    GNDVI = gndvi(green, nir; kw...)
or

    GNDVI = gndvi(cube::String; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

green Normalized diff vegetation index: more sensitive to cholorphyll than ndvi. Gitelson, A., and M. Merzlyak

GNDVI = (nir - green) / (nir + green)

$(generic_docs)
"""
gndvi(green, nir; kw...) = sp_indices(green, nir; index="GNDVI", kw...)
gndvi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="GNDVI", kw...)
gndvi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "GNDVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MNDWI = mndwi(green, swir2; kw...)
or

    MNDWI = mndwi(cube::String; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Modified Normalised Difference Water Index. Xu2006

MNDWI = (green-swir2) / (green+swir2)

$(generic_docs)
"""
mndwi(green, swir2; kw...) = sp_indices(swir2, green; index="MNDWI", kw...)
mndwi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="MNDWI", kw...)
mndwi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "MNDWI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "swir 2"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MTCI = mtci(red, redEdge1, redEdge2; kw...)

Meris Terrestrial Chlorophyll Index. Clevers and Gitelson 2013, Dash and Curran 2004

MTCI = (redEdge2-redEdge1) / (redEdge1-red)
"""
mtci(red, redEdge1, redEdge2; kw...) = sp_indices(red, redEdge1, redEdge2; index="MTCI", kw...)
mtci(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 3); index="MTCI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MCARI = mcari(green, red, redEdge1; kw...)

Modified Chlorophyll Absorption ratio index. Daughtery et al. 2000

MCARI = (redEdge1 - red - 0.2 * (redEdge1 + green)) * (redEdge1 / red)
"""
mcari(green, red, redEdge1; kw...) = sp_indices(green, red, redEdge1; index="MCARI", kw...)
mcari(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 3); index="MCARI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MSAVI = msavi(red, nir; kw...)
or

    MSAVI = msavi(cube::String; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Modified soil adjusted vegetation index. Qi 1994

MSAVI = nir + 0.5 - (0.5 * sqrt(pow(2.0 * nir + 1.0, 2) - 8.0 * (nir - (2.0 * red))))

$(generic_docs)
"""
msavi(red, nir; kw...) = sp_indices(red, nir; index="MSAVI", kw...)
msavi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="MSAVI", kw...)
msavi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "MSAVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NBRI = nbri(nir, swir3; kw...)

Normalised Burn Ratio Index. Garcia 1991

NBRI = (nir - swir3) / (nir + swir3)
"""
nbri(nir, swir3; kw...) = sp_indices(swir3, nir; index="NBRI", kw...)
nbri(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="NBRI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDVI = ndvi(red, nir; kw...)
or

    NDVI = ndvi(cube::String; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Compute the NDVI vegetation index. Input can be either the bands file names, or GMTimage objects
with the band's data.

NDVI = (nir - red) / (nir + red)

$(generic_docs)
"""
ndvi(red, nir; kw...) = sp_indices(red, nir; index="NDVI", kw...)
ndvi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="NDVI", kw...)
ndvi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDWI = ndwi(green, nir; kw...)
or

    NDWI = ndwi(cube::String; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Normalized difference water index. McFeeters 1996. NDWI => (green - nir)/(green + nir)

NDWI = (green - nir)/(green + nir)

$(generic_docs)
"""
ndwi(green, nir; kw...) = sp_indices(nir, green; index="NDWI", kw...)
ndwi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="NDWI", kw...)
ndwi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDWI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDWI2 = ndwi2(nir, swir2; kw...)
or

    NDWI2 = ndwi2(cube::String; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Normalized difference water index. Gao 1996, Chen 2005 (also known as Normalized Difference Moisture Index
NDBI and LSWI)

NDWI2 = (nir - swir2)/(nir + swir2)

$(generic_docs)
"""
ndwi2(nir, swir2; kw...) = sp_indices(swir2, nir; index="NDWI2", kw...)
ndwi2(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="NDWI2", kw...)
ndwi2(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDWI2"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["nir", "swir 2"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDREI1 = ndrei1(redEdge1, redEdge2; kw...)

Normalized difference red edge index. Gitelson and Merzlyak 1994

NDREI1 = (redEdge2 - redEdge1) / (redEdge2 + redEdge1)
"""
ndrei1(redEdge1, redEdge2; kw...) = sp_indices(redEdge1, redEdge2; index="NDREI1", kw...)
ndrei1(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="NDREI1", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDREI2 = ndrei2(redEdge1, redEdge3; kw...)

Normalized difference red edge index 2. Barnes et al 2000

NDREI2 = (redEdge3 - redEdge1) / (redEdge3 + redEdge1)
"""
ndrei2(redEdge1, redEdge3; kw...) = sp_indices(redEdge1, redEdge3; index="NDREI2", kw...)
ndrei2(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="NDREI2", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    SATVI = satvi(red, swir2, swir3; kw...)

Soil adjusted total vegetation index. Marsett 2006

SATVI = ((swir2 - red) / (swir2 + red + L)) * (1.0 + L) - (swir3 / 2.0)
"""
satvi(red, swir2, swir3; kw...) = sp_indices(red, swir2, swir3; index="SATVI", kw...)
satvi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 3); index="SATVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    SAVI = savi(red, nir; kw...)
or

    SAVI = savi(cube::String; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Soil adjusted vegetation index. Huete 1988

SAVI = (nir - red) * (1.0 + L) / (nir + red + L)

$(generic_docs)
"""
savi(red, nir; kw...) = sp_indices(red, nir; index="SAVI", kw...)
savi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="SAVI", kw...)
savi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "SAVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    SLAVI = slavi(red, nir, swir2; kw...)
or

    SLAVI = slavi(cube::String; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Specific Leaf Area Vegetation Index. Lymburger 2000

SLAVI = nir / (red + swir2)

$(generic_docs)
"""
slavi(red, nir, swir2; kw...) = sp_indices(red, nir, swir2; index="SLAVI", kw...)
slavi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 3); index="SLAVI", kw...)
slavi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "SLAVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir", "swir 2"], kw...)

# ----------------------------------------------------------------------------------------------------------
function helper_sp_indices(kwargs...)
	# Helper function that is called by two different sp_indices methods
	d = KW(kwargs)
	mask = (haskey(d, :mask))
	threshold = find_in_dict(d, [:threshold])[1]
	classes = (threshold === nothing) ? find_in_dict(d, [:classes])[1] : nothing
	(mask && threshold === nothing) && error("The `mask` option requires the `threshold=x` option")
	(classes !== nothing && length(classes) > 3) && (classes = classes[1:3]; @warn("`classes` maximum elements is 3. Clipping the others"))
	return mask, classes, threshold
end

# ----------------------------------------------------------------------------------------------------------
function sp_indices(bnd1::String, bnd2::String, bnd3::String=""; index::String="", kwargs...)
	# Compute spectral indices
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
	sp_indices(Bnd1, Bnd2, Bnd3; index=index, kwargs...)
end

function sp_indices(cube::GMT.GMTimage{UInt16, 3}, bands::Vector{Int}; index::String="", kw...)
	# This method recieves the cube and a vector with the bands list and calls the worker with @view
	mask, classes, = helper_sp_indices(kw...)	# Do this first because if it errors no point in continuing
	if (length(bands) == 2)
		o = sp_indices(@view(cube[:,:,bands[1]]), @view(cube[:,:,bands[2]]); index=index, kw...)
	else
		o = sp_indices(@view(cube[:,:,bands[1]]), @view(cube[:,:,bands[2]]), @view(cube[:,:,bands[3]]); index=index, kw...)
	end
	if (mask || classes !== nothing)
		O = mat2img(o, cube)
		O.layout, O.range[5], O.range[6] = "BRPa", 0, (mask) ? 255 : length(classes)
	else
		O = mat2grid(o, cube)
	end
	O.names = [index * " index"]
	O
end

function sp_indices(bnd1, bnd2, bnd3=nothing; index::String="", kwargs...)
	# This is the method who does the real work.
	(index == "") && error("Must select which index to compute")
	@assert size(bnd1) == size(bnd2)
	mask, classes, threshold = helper_sp_indices(kwargs...)
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
	if (isa(bnd1, GMT.GMTimage))
		if (mask || classes !== nothing)
			I = mat2img(img, proj4=bnd1.proj4, wkt=bnd1.wkt, x=bnd1.x, y=bnd1.y)
			I.epsg = bnd1.epsg;		I.layout = "BRPa"
			I.range[5], I.range[6] = 0, (mask) ? 255 : length(classes)
			return I
		else
			return mat2grid(img, bnd1)
		end
	else
		return img
	end
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