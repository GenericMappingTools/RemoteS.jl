const generic_docs = "
- The first form accepts inputs as matrices, or file names of the data bands.
- The last form is more versatile but also more complex to describe.
  - `cube`: Is the file name of a 'cube', a multi-layered file normally created with the [`cutcube`](@ref) function.
     If this file was created with band descriptions one can use the `bands` or the `bandnames` options.
  - `bands`: _cubes_ created with [`cutcube`](@ref) assign descriptions starting with \"Band1 ...\" an so on
    the other bands. So when `bands` is used we search for bands named \"Band'band[k]'\", where band[k] loops
    over all elements of the `bands` vector. WARNING: the elements order in the vector must be sorted in increasing
    wavelength numbers, _i.e._ like the example for the first form.
  - `layers`: Use this option when you are certain of the bands order in the cube or the it doesn't have a bands
    description. The selection will be made with cube[:,:,layer[1]], etc... WARNING: same warn as above.
  - `bandnames`: When we know the common designation of a band, for example \"Green\", or any part of a band
    description, for example \"NIR\", we can use that info to create a `bandnames` string vector that will be
    matched against the cube's bands descriptions.

### Kwargs
  - `threshold`: When a threshold is provided we return a GMTgrid where `vals[ij] < threshold = NaN`
  - `classes`: is a vector with up to 3 elements (class separators) and we return a  UInt8 GMTimage with the
    indices categorized into vals[ij] > classes[1] = 1; vals[ij] > classes[2] = 2; vals[ij] > classes[3] = 3 and 0 otherwise.
  - `mask`: Used together with `threshold` outputs a UInt8 GMTimage mask with `vals[ij] >= threshold = 255` and 0 otherwise
     If `mask=-1` (or any other negative number) we compute instead a mask where `vals[ij] < threshold = 255` and 0 otherwise
  - `save`: Use `save=\"file_name.ext\"` to save the result in a disk file. File format is picked from file extension.
  - `order` | `bands_order` | `rgb`: For the ``GLI``, ``TGI`` and ``VARI`` (RGB) indices, we allow to reorder the bands
    and change the expected RGB order. Pass in a string, or symbol, with the color order. For example, `order=:rbg`
	will swap the green and blue components making the result index identify the _reds_ instead of the _greens_.
	Not good for vegetation indices, but potentially useful for other purposes.

If none of `bands`, `layers` or `bandnames` is provided, we use the default band names shown in the first form.

See also https://www.indexdatabase.de/ for a list of indices and the appropriate band names per sensor.

Returns either a Float32 GMTgrid or a UInt8 GMTimage if the `mask` or `classes` options are used.
"

# ----------------------------------------------------------------------------------------------------------
# Helper function to compute Spectral Indices from a 'cube' file and somehow band selection.
function helper_si_method(cube::String, index::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[],
	                      bandnames::Vector{String}=String[], defbandnames::Vector{String}=String[], kw...)
	if (index == "GLI" || index == "TGI" || index == "VARI") && ((ext = lowercase(splitext(cube)[2])) == ".png" || ext == ".jpg")
		return sp_indices(gdalread(cube); index=index, kw...)	# Try to read it as an image
	end
	(isempty(bands) && isempty(bandnames) && isempty(layers)) && (bandnames = defbandnames)
	sc = subcube(cube, bands=bands, layers=layers, bandnames=bandnames)
	sp_indices(sc, collect(1:size(sc,3)); index=index, kw...)	# Here we know that the layers are all of those in cube
end

# Method for in memory cubes
function helper_si_method(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}}, index::String;
                          bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], defbandnames::Vector{String}=String[], kw...)
	(isempty(bands) && isempty(bandnames) && isempty(layers)) && (bandnames = defbandnames)
	lay = !isempty(layers) ? layers : find_layers(cube, bandnames, bands)
	sp_indices(cube, lay; index=index, kw...)
end

# ----------------------------------------------------------------------------------------------------------
"""
    CLG = clg(green, redEdge3; kw...)
or

    CLG = clg(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Green cholorphyl index. Wu et al 2012.

CLG = (redEdge3)/(green)-1 
"""
clg(green, redEdge3; kw...) = sp_indices(green, redEdge3; index="CLG", kw...)
clg(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "CLG"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "redEdge3"], kw...)
clg(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
    bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "CLG"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "redEdge3"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    CLRE = clre(redEdge1, redEdge3; kw...)
or

    CLRE = clre(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

RedEdge cholorphyl index. Clevers and Gitelson 2013.

CLRE = (redEdge3)/(redEdge1)-1
"""
clre(redEdge1, redEdge3; kw...) = sp_indices(redEdge1, redEdge3; index="CLRE", kw...)
clre(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "CLRE"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["redEdge1", "redEdge3"], kw...)
clre(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "CLRE"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["redEdge1", "redEdge3"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    EVI = evi(blue, red, nir; kw...)
or

    EVI = evi(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Enhanced vegetation index. Huete et al 1990

EVI = G * ((nir - red) / (nir + C1 * red - C2 * blue + Levi));
C1, C2, G, Levi = 6.0, 7.5, 2.5, 1.

$(generic_docs)

"""
evi(blue, red, nir; kw...) = sp_indices(blue, red, nir; index="EVI", kw...)
evi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "EVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["blue", "red", "nir"], kw...)
evi(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "EVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["blue", "red", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    EVI2 = evi2(red, nir; kw...)
or

    EVI2 = evi2(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Two-band Enhanced vegetation index. Jiang et al 2008

EVI2 = G * ((nir - red) / (nir + 2.4 * red ))

$(generic_docs)
"""
evi2(red, nir; kw...) = sp_indices(red, nir; index="EVI2", kw...)
evi2(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "EVI2"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)
evi2(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "EVI2"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    GNDVI = gndvi(green, nir; kw...)
or

    GNDVI = gndvi(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

green Normalized diff vegetation index: more sensitive to cholorphyll than ndvi. Gitelson, A., and M. Merzlyak

GNDVI = (nir - green) / (nir + green)

$(generic_docs)
"""
gndvi(green, nir; kw...) = sp_indices(green, nir; index="GNDVI", kw...)
gndvi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "GNDVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "nir"], kw...)
gndvi(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "GNDVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    GLI = gli(red, green, blue; kw...)
or (here fname is a .png or .jpg file name)

    GLI = gli(fname::String; kw...)
or

    GLI = gli(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Green Leaf Index. Louhaichi, M., Borman, M.M., Johnson, D.E., 2001. 

GLI = (2green - red - blue) / (2green + red + blue)

$(generic_docs)
"""
gli(rgb::GMTimage{UInt8, 3}; kw...) = sp_indices(rgb; index="GLI", kw...)
gli(red, green, blue; kw...) = sp_indices(red, green, blue; index="GLI", kw...)
gli(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "GLI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "green", "blue"], kw...)
gli(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
    bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "GLI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "green", "blue"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MNDWI = mndwi(green, swir2; kw...)
or

    MNDWI = mndwi(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Modified Normalised Difference Water Index. Xu2006

MNDWI = (green-swir2) / (green+swir2)

$(generic_docs)
"""
mndwi(green, swir2; kw...) = sp_indices(swir2, green; index="MNDWI", kw...)
mndwi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "MNDWI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "swir2"], kw...)
mndwi(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "MNDWI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "swir2"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MTCI = mtci(red, redEdge1, redEdge2; kw...)

Meris Terrestrial Chlorophyll Index. Clevers and Gitelson 2013, Dash and Curran 2004

MTCI = (redEdge2-redEdge1) / (redEdge1-red)
"""
mtci(red, redEdge1, redEdge2; kw...) = sp_indices(red, redEdge1, redEdge2; index="MTCI", kw...)
mtci(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 3); index="MTCI", kw...)
mtci(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "MTCI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "redEdge1", "redEdge2"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MCARI = mcari(green, red, redEdge1; kw...)
or

	MCARI = mcari(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Modified Chlorophyll Absorption ratio index. Daughtery et al. 2000

MCARI = (redEdge1 - red - 0.2 * (redEdge1 - green)) * (redEdge1 / red)

(Sentinel-2 Band 5 (VNIR), Band 4 (Red) and Band 3 (Green)).
"""
mcari(green, red, redEdge1; kw...) = sp_indices(green, red, redEdge1; index="MCARI", kw...)
mcari(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "MCARI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "red", "redEdge1"], kw...)
mcari(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
      bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "MCARI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "red", "redEdge1"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MSAVI = msavi(red, nir; kw...)
or

    MSAVI = msavi(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Modified soil adjusted vegetation index. Qi 1994

MSAVI = nir + 0.5 - (0.5 * sqrt(pow(2.0 * nir + 1.0, 2) - 8.0 * (nir - (2.0 * red))))

$(generic_docs)
"""
msavi(red, nir; kw...) = sp_indices(red, nir; index="MSAVI", kw...)
msavi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "MSAVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)
msavi(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "MSAVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NBRI = nbri(nir, swir3; kw...)
or

	NBRI = nbri(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Normalised Burn Ratio Index. Garcia 1991

NBRI = (nir - swir2) / (nir + swir2)
"""
nbri(nir, swir2; kw...) = sp_indices(swir2, nir; index="NBRI", kw...)
nbri(cube::String;
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NBRI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["nir", "swir2"], kw...)
nbri(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NBRI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["nir", "swir2"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDVI = ndvi(red, nir; kw...)
or

    NDVI = ndvi(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Compute the NDVI vegetation index. Input can be either the bands file names, or GMTimage objects
with the band's data.

NDVI = (nir - red) / (nir + red)

$(generic_docs)
"""
ndvi(red, nir; kw...) = sp_indices(red, nir; index="NDVI", kw...)
ndvi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)
ndvi(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDWI = ndwi(green, nir; kw...)
or

    NDWI = ndwi(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Normalized difference water index. McFeeters 1996. NDWI => (green - nir)/(green + nir)

NDWI = (green - nir)/(green + nir)

$(generic_docs)
"""
ndwi(green, nir; kw...) = sp_indices(nir, green; index="NDWI", kw...)
ndwi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDWI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "nir"], kw...)
ndwi(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDWI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["green", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDWI2 = ndwi2(nir, swir2; kw...)
or

    NDWI2 = ndwi2(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Normalized difference water index. Gao 1996, Chen 2005 (also known as Normalized Difference Moisture Index
NDBI and LSWI)

NDWI2 = (nir - swir2)/(nir + swir2)

$(generic_docs)
"""
ndwi2(nir, swir2; kw...) = sp_indices(swir2, nir; index="NDWI2", kw...)
ndwi2(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDWI2"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["nir", "swir2"], kw...)
ndwi2(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDWI2"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["nir", "swir2"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDREI1 = ndrei1(redEdge1, redEdge2; kw...)

Normalized difference red edge index. Gitelson and Merzlyak 1994

NDREI1 = (redEdge2 - redEdge1) / (redEdge2 + redEdge1)
"""
ndrei1(redEdge1, redEdge2; kw...) = sp_indices(redEdge1, redEdge2; index="NDREI1", kw...)
ndrei1(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
       bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDREI1"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["redEdge1", "redEdge2"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDREI2 = ndrei2(redEdge1, redEdge3; kw...)
or

	NDREI2 = ndrei2(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Normalized difference red edge index 2. Barnes et al 2000

NDREI2 = (redEdge3 - redEdge1) / (redEdge3 + redEdge1)
"""
ndrei2(redEdge1, redEdge3; kw...) = sp_indices(redEdge1, redEdge3; index="NDREI2", kw...)
ndrei2(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDREI2"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["redEdge1", "redEdge3"], kw...)
ndrei2(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
       bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "NDREI2"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["redEdge1", "redEdge3"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    SATVI = satvi(red, swir2, swir3; kw...)
or

	SATVI = satvi(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Soil adjusted total vegetation index. Marsett 2006

SATVI = ((swir1 - red) / (swir1 + red + L)) * (1.0 + L) - (swir2 / 2.0)
"""
satvi(red, swir2, swir3; kw...) = sp_indices(red, swir2, swir3; index="SATVI", kw...)
satvi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "SATVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "swir1", "swir2"], kw...)
satvi(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "SATVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "swir2", "swir3"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    SAVI = savi(red, nir; kw...)
or

    SAVI = savi(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Soil adjusted vegetation index. Huete 1988

SAVI = (nir - red) * (1.0 + L) / (nir + red + L)

$(generic_docs)
"""
savi(red, nir; kw...) = sp_indices(red, nir; index="SAVI", kw...)
savi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "SAVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)
savi(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "SAVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    SLAVI = slavi(red, nir, swir2; kw...)
or

    SLAVI = slavi(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Specific Leaf Area Vegetation Index. Lymburger 2000

SLAVI = nir / (red + swir2)

$(generic_docs)
"""
slavi(red, nir, swir2; kw...) = sp_indices(red, nir, swir2; index="SLAVI", kw...)
slavi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "SLAVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir", "swir2"], kw...)
slavi(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
     bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "SLAVI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "nir", "swir2"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    TGI = tgi(red, green, blue; kw...)
or (here fname is a .png or .jpg file name)

    TGI = tgi(fname::String; kw...)
or

    TGI = tgi(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Triangular Greenness Index. Hunt et al. 2013

TGI = green - 0.39 * red - 0.61 * blue

$(generic_docs)
"""
tgi(rgb::GMTimage{UInt8, 3}; kw...) = sp_indices(rgb; index="TGI", kw...)
tgi(red, green, blue; kw...) = sp_indices(red, green, blue; index="TGI", kw...)
tgi(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "TGI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "green", "blue"], kw...)
tgi(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
    bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "TGI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "green", "blue"], kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    VARI = vari(red, green, blue; kw...)
or (here fname is a .png or .jpg file name)

    VARI = vari(fname::String; kw...)
or

    VARI = vari(cube::Union{String, GMTgrid}; [bands=Int[], bandnames=String[], layers=Int[]], kwargs...)

Visible Atmospherically Resistant Index. Gitelson, A.A., Kaufman, Y.J., Stark, R., Rundquist, D., 2002

VARI = (green - red) / (green + red - blue)

$(generic_docs)
"""
vari(rgb::GMTimage{UInt8, 3}; kw...) = sp_indices(rgb; index="VARI", kw...)
vari(red, green, blue; kw...) = sp_indices(red, green, blue; index="VARI", kw...)
vari(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "VARI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "green", "blue"], kw...)
vari(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
    bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], kw...) =
	helper_si_method(cube, "VARI"; bands=bands, layers=layers, bandnames=bandnames, defbandnames=["red", "green", "blue"], kw...)

# ----------------------------------------------------------------------------------------------------------
function helper_sp_indices(kwargs...)
	# Helper function that is called by two different sp_indices methods
	d = KW(kwargs)
	mask::Bool = ((val = find_in_dict(d, [:mask], false)[1]) !== nothing)
	rev_mask = (mask && val < 0) ? true : false		# See if we want to reverse the mask
	threshold = find_in_dict(d, [:threshold], false)[1]
	classes = (threshold === nothing) ? find_in_dict(d, [:classes])[1] : nothing
	(mask && threshold === nothing) && error("The `mask` option requires the `threshold=x` option")
	(classes !== nothing && length(classes) > 3) && (classes = classes[1:3]; @warn("`classes` maximum elements is 3. Clipping the others"))
	save_name::String = ((val = find_in_dict(d, [:save], false)[1]) !== nothing) ? string(val) : ""
	dbg = (find_in_dict(d, [:Vd :dbg])[1] !== nothing)
	dd = d		# Copy to report unused keys (errors) but can't do it in 'd' because this function may be called twice
	(haskey(dd, :mask))      && delete!(dd, :mask)
	(haskey(dd, :classes))   && delete!(dd, :classes)
	(haskey(dd, :threshold)) && delete!(dd, :threshold)
	(haskey(dd, :rgb))       && delete!(dd, :rgb)
	(haskey(dd, :order))     && delete!(dd, :order)
	(haskey(dd, :bands_order)) && delete!(dd, :bands_order)
	(length(d) > 0) && println("Warning: the following options were not consumed in sp_indices => ", keys(dd))
	return mask, rev_mask, classes, threshold, save_name, dbg
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

# ----------------------------------------------------------------------------------------------------------
function sp_indices(rgb::GMT.GMTimage{UInt8, 3}; index::String="", kw...)
	# This method applyies only in the case of the RGB vegetation indices (GLI, TGI, VARI)
	(index != "GLI" && index != "TGI" && index != "VARI") && error("With RGB images input, only `GLI`, `TGI` and `VARI` indices are supported, not $index")
	(rgb.layout[3] != 'B') && error("For now, only band interleavedRGB composition is supported and not $(rgb.layout)")
	# Here we are allowing cheating the indices by altering the bands order. These indices expect (were deffined)
	# the bands in RGB order but nothing stops us to to change that and convert a green index into a blue index.
	# For that pass string with the R,G,B in the wished order to the 'order' option. E.g. 'order="rbg"'
	bds::Vector{Int} = [1,2,3]			# The default RGB order
	if ((val = find_in_kwargs(kw, [:order :bands_order :rgb])[1]) !== nothing)
		o = lowercase(string(val))
		bds[1] = (o[1] == 'r') ? 1 : (o[1] == 'g') ? 2 : (o[1] == 'b') ? 3 : error("Non 'r', 'g' or 'b' in the 'order' option")
		bds[2] = (o[2] == 'r') ? 1 : (o[2] == 'g') ? 2 : (o[2] == 'b') ? 3 : error("Non 'r', 'g' or 'b' in the 'order' option")
		bds[3] = (o[3] == 'r') ? 1 : (o[3] == 'g') ? 2 : (o[3] == 'b') ? 3 : error("Non 'r', 'g' or 'b' in the 'order' option")
	end
	img = sp_indices(view(rgb, :, :, bds[1]), view(rgb, :, :, bds[2]), view(rgb, :, :, bds[3]); index=index, kw...)
	return mat2grid(img, rgb)
end

# ----------------------------------------------------------------------------------------------------------
function sp_indices(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}}, bands::Vector{Int}; index::String="", kw...)
	# This method recieves the cube and a vector with the bands list and calls the worker with @view
	mask, _, classes, _, save_name, dbg = helper_sp_indices(kw...)	# Do this first because if it errors no point in continuing
	(dbg) && println(cube.names)
	if (length(bands) == 2)
		o = sp_indices(@view(cube[:,:,bands[1]]), @view(cube[:,:,bands[2]]); index=index, kw...)
	else
		o = sp_indices(@view(cube[:,:,bands[1]]), @view(cube[:,:,bands[2]]), @view(cube[:,:,bands[3]]); index=index, kw...)
	end
	if (mask || classes !== nothing)
		if (isa(o, Matrix{<:Integer}))
			# Don't know if this case still happens. It was for old code. Now, for masks, 'o' is already a GMTimage
			O = mat2img(o, cube)
			O.layout, O.range[5], O.range[6] = "BRPa", 0, (mask) ? 255 : length(classes)
		else
			O = o
		end
	else
		O = isa(o, GMT.GMTgrid) ? o : mat2grid(o, cube)
	end
	O.names = [index * " index"]
	(save_name != "") && gmtwrite(save_name, O)
	return (save_name == "") ? O : nothing
end

# ----------------------------------------------------------------------------------------------------------
function sp_indices(bnd1, bnd2, bnd3=nothing; index::String="", kwargs...)
	# This is the method who does the real work.
	(index == "") && error("Must select which index to compute")
	@assert size(bnd1) == size(bnd2)
	ismask, rev_mask, classes, threshold, = helper_sp_indices(kwargs...)
	if ((ismask || classes !== nothing))
		mask = fill(UInt8(0), size(bnd1))
	else
		img = fill(0.0f0, size(bnd1))		# Changed mind. Make 0.0f0 the neutral instead of NaN. Classification algos don't work with NaN
	end

	if (rev_mask) fcomp = <	else  fcomp = >  end 	# Which one to use when computing masks
	mn = size(bnd1,1) * size(bnd1,2)
	i_tmax = (eltype(bnd1) <: Integer) ? 1. / typemax(eltype(bnd1)) : 1.0	# Floats go unchanged
	C1, C2, G, L, Levi = 6.0, 7.5, 2.5, 0.5, 1.0
	if (index == "CLG" || index == "CLRE")
		# Green cholorphyl index. Wu et al 2012		CLG               (redEdge3)/(green)-1
		# RedEdge cholorphyl index. Clevers and Gitelson 2013	CLRE  (redEdge3)/(redEdge1)-1 
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t = bnd2[k] / bnd1[k] - 1
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = bnd2[k] / bnd1[k] - 1
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "EVI")			# Enhanced vegetation index. Huete et al 1990
		# G * ((nir - red) / (nir + C1 * red - C2 * blue + Levi))
		(bnd3 === nothing) && error("NIR component cannot be empty")
		@assert size(bnd3) == size(bnd2)
		blue = bnd1;	red = bnd2;		nir = bnd3
		_C2 = C2 * i_tmax;
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t_red = red[k]*i_tmax;	t_nir = nir[k]*i_tmax
				t = G * ((t_nir - t_red) / (t_nir + C1 * t_red - _C2 * blue[k] + Levi))
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t_red = red[k]*i_tmax;	t_nir = nir[k]*i_tmax
				img[k] = G * ((t_nir - t_red) / (t_nir + C1 * t_red - _C2 * blue[k] + Levi))
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "EVI2")		# Two-band Enhanced vegetation index. Jiang et al 2008 
		# G * ((nir - red) / (nir + 2.4 * red ))
		red = bnd1;		nir = bnd2
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t_red = red[k]*i_tmax;	t_nir = nir[k]*i_tmax
				t = G * (t_nir - t_red) / (t_red + 2.4 * t_nir)
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t_red = red[k]*i_tmax;	t_nir = nir[k]*i_tmax
				t = G * (t_nir - t_red) / (t_red + 2.4 * t_nir)
				(t >= -1 && t <= 1) && (img[k] = t)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "GLI")			# Green Leaf Index 
		# (2green - red - blue) / (2green + red + blue)
		red = bnd1;		green = bnd2;		blue = bnd3
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t_red  = red[k]*i_tmax;	t_green = 2 * green[k]*i_tmax;	t_blue = blue[k]*i_tmax
				t_rb   = t_red + t_blue
				t =  (t_green - t_rb) / (t_green + t_rb)
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t_red  = red[k]*i_tmax;	t_green = 2 * green[k]*i_tmax;	t_blue = blue[k]*i_tmax
				t_rb   = t_red + t_blue
				img[k] =  (t_green - t_rb) / (t_green + t_rb)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "MTCI")		# Meris Terrestrial Chlorophyll Index. Clevers and Gitelson 2013, Dash and Curran 2004 
		# (redEdge2-redEdge1) / (redEdge1-red)
		red = bnd1;		redEdge1 = bnd2;	redEdge2 = bnd3
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t_redEdge1 = redEdge1[k]*i_tmax
				t = (redEdge2[k]*i_tmax - t_redEdge1) / (t_redEdge1 - red[k]*i_tmax)
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t_redEdge1 = redEdge1[k]*i_tmax
				img[k] = (redEdge2[k]*i_tmax - t_redEdge1) / (t_redEdge1 - red[k]*i_tmax)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "MCARI")		# Modified Chlorophyll Absorption ratio index. Daughtery et al. 2000 
		# (redEdge1 - red - 0.2 * (redEdge1 - green)) * (redEdge1 / red)
		# Sentinel-2 Band 5 (VNIR), Band 4 (Red) and Band 3 (Green).
		green = bnd1;	red = bnd2;		redEdge1 = bnd3
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t_redEdge1 = redEdge1[k]*i_tmax;	t_red = red[k]*i_tmax
				t = (t_redEdge1 - t_red - 0.2 * (t_redEdge1 - green[k]*i_tmax)) * (t_redEdge1 / t_red)
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t_redEdge1 = redEdge1[k]*i_tmax;	t_red = red[k]*i_tmax
				img[k] = (t_redEdge1 - t_red - 0.2 * (t_redEdge1 - green[k]*i_tmax)) * (t_redEdge1 / t_red)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "MSAVI")		# Modified soil adjusted vegetation index.
		# nir + 0.5 - (0.5 * sqrt(pow(2.0 * nir + 1.0, 2) - 8.0 * (nir - (2.0 * red))))
		red = bnd1;		nir = bnd2;
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t_nir = nir[k]*i_tmax
				t = t_nir + 0.5 - (0.5 * sqrt((2 * t_nir + 1) ^2) - 8 * (t_nir - (2 * red[k]*i_tmax)))
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t_nir = nir[k]*i_tmax
				img[k] = t_nir + 0.5 - (0.5 * sqrt((2 * t_nir + 1) ^2) - 8 * (t_nir - (2 * red[k]*i_tmax)))
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

		# Need to swap bnd1 and bnd2 for indices that are defined as (shorter_cdo - longer_cdo)
		_b1, _b2 = (index == "MNDWI" || index == "NBRI" || index == "NDWI" || index == "NDWI2" || index == "NDWI2") ? (bnd2, bnd1) : (bnd1, bnd2)
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t1 = _b1[k]*i_tmax;	t2 = _b2[k]*i_tmax
				t = (t2 - t1) / (t1 + t2)
				(fcomp(t, threshold) && t <= 1) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t1 = _b1[k]*i_tmax;	t2 = _b2[k]*i_tmax
				t = (t2 - t1) / (t1 + t2)
				(t >= -1 && t <= 1) && (img[k] = t)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "SATVI")		# Soil adjusted total vegetation index.
		# ((swir1 - red) / (swir1 + red + L)) * (1.0 + L) - (swir2 / 2.0)
		red = bnd1;		swir1 = bnd2;	swir2 = bnd3;
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t_red = red[k]*i_tmax;	t_swir1 = swir1[k]*i_tmax
				t = ((t_swir1 - t_red) / (t_swir1 + t_red + L)) * (1.0 + L) - (swir2[k]*i_tmax / 2.0)
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t_red = red[k]*i_tmax;	t_swir1 = swir1[k]*i_tmax
				img[k] = ((t_swir1 - t_red) / (t_swir1 + t_red + L)) * (1.0 + L) - (swir2[k]*i_tmax / 2.0)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "SAVI")		# Soil adjusted vegetation index. Huete1988
		# (nir - red) * (1.0 + L) / (nir + red + L);
		red = bnd1;		nir = bnd2;
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t_red = red[k]*i_tmax;	t_nir = nir[k]*i_tmax
				t = (t_nir - t_red) * (1.0 + L) / (t_nir - t_red + L)
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t_red = red[k]*i_tmax;	t_nir = nir[k]*i_tmax
				img[k] = (t_nir - t_red) * (1.0 + L) / (t_nir - t_red + L)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "SLAVI")		# Soil Adjusted Vegetation Index Huete 1988.
		# nir / (red + swir2)
		red = bnd1;		nir = bnd2;		swir2 = bnd3
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t = nir[k]*i_tmax / (red[k]*i_tmax + swir2[k]*i_tmax)
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = nir[k]*i_tmax / (red[k]*i_tmax + swir2[k]*i_tmax)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "TGI")			# Triangular Greenness Index. Hunt et al. 2013.
		# green - 0.39 * red - 0.61 * blue
		red = bnd1;		green = bnd2;		blue = bnd3
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t = green[k]*i_tmax - 0.39 * red[k]*i_tmax - 0.61 * blue[k]*i_tmax
				(fcomp(t, threshold)) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				img[k] = green[k]*i_tmax - 0.39 * red[k]*i_tmax - 0.61 * blue[k]*i_tmax
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	elseif (index == "VARI")			# Visible Atmospherically Resistant Index.
		#  (green - red) / (green + red - blue)
		red = bnd1;		green = bnd2;		blue = bnd3
		if (ismask)
			@inbounds Threads.@threads for k = 1:mn
				t_green = green[k]*i_tmax;	t_red = red[k]*i_tmax
				t = (t_green - t_red) / (t_green + t_red - blue[k]*i_tmax)
				(fcomp(t, threshold) && t <= 1) && (mask[k] = 255)
			end
		else
			@inbounds Threads.@threads for k = 1:mn
				t_green = green[k]*i_tmax;	t_red = red[k]*i_tmax
				t = (t_green - t_red) / (t_green + t_red - blue[k]*i_tmax)
				(t >= -1 && t <= 1) && (img[k] = t)
			end
			helper_si!(img, threshold, classes)		# Threshold or Classes if one of them is != nothing
		end
	end

	if (isa(bnd1, GMT.GMTimage) || isa(bnd1, GMT.GMTgrid) || isa(parent(bnd1), GItype))	# Last one is when slice is a SubArray{...
		if (ismask || classes !== nothing)
			I = mat2img(mask, isa(bnd1, GMT.GMTimage) ? bnd1 : parent(bnd1))
			I.layout="BRPa"
			I.range[5], I.range[6] = 0, (ismask) ? 255 : length(classes)
			return I
		else
			return mat2grid(img, isa(bnd1, GMT.GItype) ? bnd1 : parent(bnd1))
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
