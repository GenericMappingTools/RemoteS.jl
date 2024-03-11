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

const Lsat8_desc = Dict(
	1  => "Band1 - Coastal aerosol [0.43-0.45]",
	2  => "Band2 - Blue [0.45-0.51]",
	3  => "Band3 - Green [0.53-0.59]",
	4  => "Band4 - Red [0.64-0.67]",
	5  => "Band5 - NIR [0.85-0.88]",
	6  => "Band6 - SWIR1 [1.57-1.65]",
	7  => "Band7 - SWIR2 [2.11-2.29]",
	9  => "Band9 - Cirrus [1.36-1.38]",
	10 => "Band10 - Thermal IR1 [10.6-11.19]",
	11 => "Band11 - Thermal IR2 [11.50-12.51]")

const Sentinel2_10m_desc = Dict(
	2  => "Band2 - Blue [0.490]",
	3  => "Band3 - Green [0.560]",
	4  => "Band4 - Red [0.665]",
	8  => "Band8 - NIR [0.842]")

const Sentinel2_20m_desc = Dict(				# Common to both 20 & 60 m
	1  => "Band1 - Coastal aerosol [0.443]",	# Actually only available in 60 m.
	2  => "Band2 - Blue [0.490]",
	3  => "Band3 - Green [0.560]",
	4  => "Band4 - Red [0.665]",
	5  => "Band5 - Red Edge1 [0.705]",
	6  => "Band6 - Red Edge2 [0.740]",
	7  => "Band7 - Red Edge3 [0.783]",
	8  => "Band8A - Red Edge4 (NIR) [0.865]",	# ~Landsat8 Band 5
	9  => "Band9 - Water vapour [0.945]",
	10 => "Band10 - Cirrus [1.375]",			# ~Landsat8 Band 9
	11 => "Band11 - SWIR1 [1.610]",				# ~Landsat8 Band 6
	12 => "Band12 - SWIR2 [2.190]")				# ~Landsat8 Band 7

# ----------------------------------------------------------------------------------------------------------
"""
    Irgb = truecolor(bndR, bndG, bndB)

Take three Landsat8/Sentinel2 UINT16 GMTimages or the file names of those bands and compose
an RGB true color image applying automatic histogram stretching.

Return an UInt8 RGB GMTimage

    Irgb = truecolor(cube::GMTImage, bands::Vector{Int})

Make an RGB composition of the 3 bands passed in the vector 'bands' from the layers in the multi-layered GMTimage `cube`

Return an auto-stretched UInt8 RGB GMTimage

    Irgb = truecolor(cube::String, [bands::Vector{Int}], [bandnames::Vector{String}], [raw=false])

Make an RGB composition of 3 bands from the `cube` file holding a UInt16 multi-layered array (often created with `cutcube`)
The band selection can be made with `bands` vector, case in which we will search for bands named "Band band[k]"
or where the bands description contain the contents of `bandnames`. If none of `bands` or `bandnames` is used
we search for a made up `bandnames=["red", "green", "blue"]`.

Return an auto-stretched UInt8 RGB GMTimage OR a GMTimage{UInt16,3} if the `raw` option is set to `true`.

    Irgb = truecolor(cube::GMTgrid, [bands|layers::Vector{Int}], [bandnames::Vector{String}], [type=UInt8])

Make an RGB composition of 3 bands from the `cube` file holding a Float32 multi-layered array.
The band selection can be made with `bands` vector, case in which we will search for bands named "Band band[k]"
or where the bands description contain the contents of `bandnames`. If none of `bands` or `bandnames` is used
we search for a made up `bandnames=["red", "green", "blue"]`.

By default we scale the bands to 0-255. Use `type=UInt16` to scale the bands to 0-65535`. Note that this
will matter only for the guessing of the good limits to perform the histogram stretching.


### Example:
Make an RGB composite from data in the cube file "LC08__cube.tiff"
```julia
I = truecolor("LC08__cube.tiff");
```
"""
function truecolor(bndR, bndG, bndB)
	I = isa(bndR, GMT.GMTimage) ? bndR : gmtread(bndR)
	img = Array{UInt8}(undef, size(I,1), size(I,2), 3)
	if (eltype(I) == UInt8)
		img[:,:,1] .= I.image
	else
		_ = mat2img(I.image, stretch=true, img8=view(img,:,:,1), scale_only=1)
	end
	I = isa(bndG, GMT.GMTimage) ? bndG : gmtread(bndG)
	@assert size(I,1) == size(img,1) && size(I,2) == size(img,2)
	if (eltype(I) == UInt8)
		img[:,:,2] .= I.image
	else
		_ = mat2img(I.image, stretch=true, img8=view(img,:,:,2), scale_only=1)
	end
	I = isa(bndB, GMT.GMTimage) ? bndB : gmtread(bndB)
	@assert size(I,1) == size(img,1) && size(I,2) == size(img,2)
	if (eltype(I) == UInt8)
		img[:,:,3] .= I.image
	else
		_ = mat2img(I.image, stretch=true, img8=view(img,:,:,3), scale_only=1)
	end
	Io = mat2img(img, I);
	Io.layout = (isa(bndR, GMT.GMTimage)) ? "T" * bndR.layout[2] * "Ba" : "TRBa"	# This is shitty fragile
	(isa(bndR, GMT.GMTimage) && startswith(bndR.layout, "BC")) && (Io.layout = "BCBa")	# Horrible patch that needs to know why.
	Io
end

truecolor(cube::GMT.GMTimage{UInt16, 3}, layers::Vector{Int}) = truecolor(cube, layers=layers)
function truecolor(cube::GMT.GMTimage{UInt16, 3}; layers::Vector{Int}=Int[])
	(length(layers) != 3) && error("For an RGB composition 'bands' must be a 3 elements array and not $(length(layers))")
	(cube.layout[3] != 'B') && error("For an RGB composition the image object must be Band interleaved and not $(cube.layout)")
	img = Array{UInt8, 3}(undef, size(cube,1), size(cube,2), 3)
	layers = find_layers(cube, layers, 3)
	_ = mat2img(@view(cube[:,:,layers[1]]), stretch=true, img8=view(img,:,:,1), scale_only=1)
	_ = mat2img(@view(cube[:,:,layers[2]]), stretch=true, img8=view(img,:,:,2), scale_only=1)
	_ = mat2img(@view(cube[:,:,layers[3]]), stretch=true, img8=view(img,:,:,3), scale_only=1)
	Io = mat2img(img, cube);	Io.layout = "TRBa"
	Io
end

truecolor(cube::GMT.GMTgrid{Float32, 3}, layers::Vector{Int}) = truecolor(cube, layers=layers)
function truecolor(cube::GMT.GMTgrid{Float32, 3}; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[],
	               bandnames::Vector{String}=String[], type::DataType=UInt8)
	(isempty(bands) && isempty(bandnames) && isempty(layers)) && (bandnames = ["red", "green", "blue"])
	isempty(layers) && (layers = find_layers(cube, bandnames, bands))
	(length(layers) != 3) && error("For an RGB composition 'bands' must be a 3 elements array and not $(length(layers))")
	img = Array{type, 3}(undef, size(cube,1), size(cube,2), 3)
	img[:,:,1] = rescale(@view(cube[:,:,layers[1]]), stretch=true, type=type)
	img[:,:,2] = rescale(@view(cube[:,:,layers[2]]), stretch=true, type=type)
	img[:,:,3] = rescale(@view(cube[:,:,layers[3]]), stretch=true, type=type)
	Io = mat2img(img, cube);	Io.layout = "TRBa"
	Io
end

function truecolor(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], raw::Bool=false)
	# The `raw` option returns a GMTimage{UInt16, 3} and does not convert to UInt8 with auto-stretch (default)
	(isempty(bands) && isempty(bandnames) && isempty(layers)) && (bandnames = ["red", "green", "blue"])
	rgb = subcube(cube, bands=bands, layers=layers, bandnames=bandnames)
	(raw) && return rgb
	return (eltype(rgb) <: AbstractFloat) ? truecolor(rgb, layers=[1,2,3]) : mat2img(rgb, stretch=:auto)
end

# ----------------------------------------------------------------------------------------------------------
"""
    subcube(cube::String; bands=Int[], bandnames=String[], layers=Int[])

Extracts a subcube from `cube` with the layers in the `bands` vector, case in which we will search for bands
named "Band band[k]", or those whose names correspond (even partially and case insensitive) to the descriptions
in `bandnames` string vector. This means that the options `bands` and `bandnames` can only be used in 'cubes'
with bands description. The `layers` option blindly extract the `cube` planes listed in the `layer` vector.

Returns a GMTimage

    subcube(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}}; bands=Int[], bandnames=String[], layers=Int[])

Does the same but from an already in memory cube. Returns a type equal to the input type. No views, a data copy.

### Example
Extracts the Red, Green and Blue layers from a Landsat 8 cube created with `cutcube`

```
Irgb = subcube("LC08__cube.tiff", bandnames = ["red", "green", "blue"])
```
"""
function subcube(cube::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[])
	# This is also used as a helper function in 'truecolor' and the radiometric indices functions.
	#(isempty(bands) && isempty(bandnames) && isempty(layers)) && (bands = collect(1:length(reportbands(cube))))	# Read them all
	alllayers = (isempty(bands) && isempty(bandnames) && isempty(layers)) ? true : false	# Read them all
	(isempty(layers)) && (layers = find_layers(cube, bands=bands, bandnames=bandnames, alllayers=alllayers)[1])
	gmtread(cube, band=layers, layout="TRBa")		# This one is still UInt16
end

function subcube(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}};
                 bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[])
	(isempty(bands) && isempty(bandnames) && isempty(layers)) && return cube	# Stupid but silently ignore it.
	(isempty(layers)) && (layers = find_layers(cube, bandnames, bands))
	slicecube(cube, layers)
end

# ----------------------------------------------------------------------------------------------------------
function find_layers(cube::Union{GMT.GMTimage{UInt16, 3}, AbstractArray{<:AbstractFloat, 3}}, list::Vector{Int}, n_layers::Int)
	# This function is not finished
	if (maximum(list) < 200)		# The list of bands to pass to the caling fun
		(maximum(list) > size(cube,3)) && error("Not enough bands to satisfy the bands list request.")
		bands = list
	#=
	elseif (!isempty(cube.v))		# Must match the wavelength in `list` with those of cube.v
		bands, n = zeros(Int, length(list)), 0
		for w in list
			d, ind = findmin(abs.(w .- cube.v))
			if (d < tol)
				bands[n += 1] = ind
				continue
			end
		end
		(n != length(list)) && error("Some of the wavelegth in $(list) are not present in the `cube.v`")
		bands
	=#
	else
		error("The `cube` object does not have a frequencies (`v` coordinates) vector as required here.")
	end
	(length(bands) != n_layers) && error("Need $(n_layers) bands but got $(length(bands))")
	bands
end

#function find_layers(cube::GMT.GMTimage{UInt16, 3}, bandnames::Vector{String}=String[], bands::Vector{Int}=Int[])
function find_layers(cube::AbstractArray, bandnames::Vector{String}=String[], bands::Vector{Int}=Int[])
	# Find the layers corresponding to the (parts of) contents of "bandnames"
	(!isa(cube, GMT.GMTimage{UInt16, 3}) && !isa(cube, GMT.GMTgrid{Float32, 3})) && error("'cube' must be a 3D GMTimage or GMTgrid")
	(isempty(cube.names)) && error("The `cube` object does not have a `names` (band names) assigned field as required here.")
	(!isempty(bands)) && (bandnames = ["Band $(bands[k])" for k = 1:length(bands)])		# Create a bandnames vector
	helper_find_layers(lowercase.(cube.names), bandnames)
end

function find_layers(fname::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], alllayers::Bool=false)
	# 'layers' just return itself after checking that the cube file actually contain that many layers
	# 'bands' search for "Band bands[k]"
	# 'bandnames' search the bands description for the first layer that contains bandnames[k]
	# Returns the numeric layers (1-based) corresponding to the search criteria and the bands description
	(!alllayers && isempty(bands) && isempty(layers) && isempty(bandnames)) &&
		error("Must use either the 'bands' OR the 'bandnames' option.")
	GMT.ressurectGDAL()		# Yep, shit, sometimes its needed.
	ds = GMT.Gdal.unsafe_read(fname)
	nbands = GMT.Gdal.nraster(ds)
	msg = ""
	(nbands < 2) && (msg = "This file ($fname) does not contain cube data (more than one layer).")
	(!isempty(layers) && maximum(layers) > nbands) && (msg = "Asked for more 'layers' than this cube contains.")
	(msg != "") && (GMT.Gdal.GDALClose(ds.ptr); error(msg))
	desc = fill("", nbands)
	[desc[k] = GMT.Gdal.GDALGetDescription(GMT.Gdal.GDALGetRasterBand(ds.ptr, k)) for k = 1:nbands]
	GMT.Gdal.GDALClose(ds.ptr)

	(!isempty(layers)) && return layers, desc
	(all(desc .== "")) && error("This cube file has no band descriptions so cannot use the 'band' or 'bandnames' options.")
	(alllayers) && return collect(1:nbands), desc	# OK, just return them ALL

	(!isempty(bands)) && (bandnames = ["Band$(bands[k])" for k = 1:length(bands)])		# Create a bandnames vector
	_layers = helper_find_layers(lowercase.(desc), bandnames)
	return _layers, desc
end

helper_find_layers(desc::Vector{String}, band_names::String) = helper_find_layers(desc, [band_names])
function helper_find_layers(desc::Vector{String}, band_names::Vector{String})::Vector{Int}
	# Note, 'desc' is supposed to be in lowercase already.
	bnd_names = lowercase.(band_names)
	layers, n = zeros(Int, length(bnd_names)), 0
	for bnd_name in bnd_names
		((ind = findfirst(contains.(desc, bnd_name))) !== nothing) && (layers[n += 1] = ind)
	end
	(n != length(bnd_names)) && error("All or some of the names in $(band_names) are not present in the `cube` names/description.\n\tUse the `reportbands()` function to check the bands description.")
	layers
end

# ----------------------------------------------------------------------------------------------------------
"""
    reportbands(in; [layers=Int[]])
or

    reportbands(in, layer;)

Report the Bands description of the `in` input argument. This can be a GMTimage, a GMTgrid or a file name (a String) of 
a 'cube' file. Normally one made with the `cutcube` function. When the use conditions of this function are not met,
either a warning or an error message (if too deep to be caught as a warning) will be issued.

- `layers`: When this optional parameter is used, report the description of the bands in the vector `layers`
- `layer`: A scalar with a unique band number. Alternative form to `reportbands(in, layers=[layer])`

Returns a string vector.
"""
reportbands(in, layer::Int) = reportbands(in, layers=[layer])
function reportbands(in; layers::Vector{Int}=Int[])
	(!isa(in, GMT.GMTimage) && !isa(in, GMT.GMTgrid) && !isa(in, String)) &&
		error("Bad input. Must be a GMTimage, a GMTgrid or a file name. Not $(typeof(in))")
	if (isa(in, GMT.GMTimage) || isa(in, GMT.GMTgrid))
		isempty(in.names) && (println("This image object does not have a `names` assigned field"); return nothing)
		return (isempty(layers)) ? in.names : in.names[layers]
	else
		layers, desc = find_layers(in, layers=layers, alllayers=isempty(layers))
		return desc[layers]
	end
end

# ----------------------------------------------------------------------------------------------------------
"""
    cutcube(names=String[], bands=Int[], template="", region=nothing, extension=".TIF", description=String[], mtl="", sentinel2=0, save="")

Cut a 3D cube out of a Landsat/Sentinel scene within a subregion `region` and a selection of bands.

- `names`: (optional) A vector with the individual bands full file name
- `bands`: When `names` is not provided give a vector of integers corresponding to the choosen bands.
           This works well for Landsat and most of Sentinel bands. However, in later case, there are also
           bands that contain characters, for example band 8A. In this case `bands` should be a vector of
           strings including the extension. _e.g._ ["02.jp2", "8A.jp2"]
- `template`: Goes together with the `bands` option. They are both composed a template * band[n] to recreate
           the full file name of each band.
- `region` Is the region to extract and must contain the extracting region limits as [W, E, S, N] or a
           GMT style -R string (without the leading "-R").
- `extension`: In case the `bands` is numeric but file extensions are not "*.TIF" (case insensitive),
           use the extension passed by this option.
- `description`: A vector of strings (as many as bands) with a description for each band. If not provided and
           the file is recognized as a Landasat 8, band description is added automatically, otherwise
           we build one with the bands file names. This info will saved if data is written to a file.
- `mtl`:   If reading from Landsat and the MTL file is not automatically found (you get an error) use this
           option to pass the full name of the MTL file.
- `sentinel2`: ESA is just unconsistent and names change with time and band numbers can have character (e.g. 8A)
           hence we need help to recognize Sentinel files so the known description can be assigned.
           Use `sentinel=10`, or `=20` or `=60` to indicate Sentinel files at those resolutions.
- `save`:  The file name where to save the output. If not provided, a GMTimage is returned.

Return: `nothing` if the result is written in file or a GMTimage otherwise.

## Examples  
  
```julia
# Cut a Landsat 8 scene for a small region (in UTM) and return a GMTimage with 3 bands in UInt16.
temp = "C:\\SIG_AnaliseDadosSatelite\\SIG_ADS\\DadosEx2\\LC82040332015145LGN00\\LC82040332015145LGN00_B";
cube = cutcube(bands=[2,3,4], template=temp, region=[479670,492720,4282230,4294500])

# The same example as above but save the data in a GeoTIFF disk file and use a string for `region`
cutcube(bands=[2,3,4], template=temp, region="479670/492720/4282230/4294500", save="landsat_cube.tif")
```
"""
function cutcube(; names::Vector{String}=String[], bands::AbstractVector=Int[], template::String="",
                   region=nothing, extension::String=".TIF", description::Vector{String}=String[], save::String="", sentinel2::Int=0, mtl::String="")

	if (isempty(names))
		(isempty(bands) || template == "") && error("When band file `names` are not provided, MUST indicate `bands` AND `template`")
		!isa(bands, Vector{<:Integer}) && !isa(bands, Vector{<:String}) && error("`bands` must a vector of Int or Strings.")
		names = Vector{String}(undef, length(bands))
		if isa(bands, Vector{Int})
			[names[k] = @sprintf("%s%d%s", template, bands[k], extension) for k = 1:length(bands)]
			!isfile(names[1]) &&	# Landsat uses "B2.TIF" and Sentinel "B02.jp2", try again.
				[names[k] = @sprintf("%s%.02d%s", template, bands[k], extension) for k = 1:length(bands)]
		else
			[names[k] = @sprintf("%s%s", template, bands[k]) for k = 1:length(bands)]
		end
	end
	# Now test if any of the file names, either passed in or generated here, do not exist
	k = 1
	for name in names
		if (!isfile(name))
			if (startswith(name, "LC0") && template[end-3:end-2] == "SR")	# For Landsat 9 termal bands are now _ST_B10
				name = replace(name, "_SR_B" => "_ST_B")
				if (isfile(name))
					names[k] = name
				else
					error("File name $name does not exist. Must stop here.")
				end
			else
				error("File name $name does not exist. Must stop here.")
			end
		end
		k += 1
	end
	MTL = read_mtl(names[1], mtl, get_full=true)
	(MTL !== nothing) && (MTL = ["MTL=" * join(MTL, "\n")])

	# Little parsing of the -R string but does not test if W < E & S < N
	(region === nothing) && (region = grdinfo(names[1], C=true)[1:4])		# Swallow the entire region
	_region::String = isa(region, String) ? region :
	                  @sprintf("%.12g/%.12g/%.12g/%.12g", region[1], region[2], region[3], region[4])
	startswith(_region, "-R") && (_region = _region[3:end])		# Tolerate a region that starts with "-R"

	desc = assign_description(names, description, sentinel2)
	cube = grdcut(names[1], R=_region)
	mat = isa(cube, GMTimage) ? cube.image : cube.z
	for k = 2:length(bands)
		B = grdcut(names[k], R=_region)
		mat = cat(mat, isa(cube, GMTimage) ? B.image : B.z, dims=3)
	end
	cube = isa(cube, GMTimage) ? mat2img(mat, cube, names=desc) : mat2grid(mat, reg=cube.registration, x=cube.x, y=cube.y, proj4=cube.proj4, wkt=cube.wkt, names=desc)
	if (save != "")
		_, ext = splitext(save)
		if (lowercase(ext) == "nc")		# Let save as a nc cube as well
			gdalwrite(cube, save, bands, dim_name="bands")
		else
			gdaltranslate(cube, dest=save, meta=MTL)
		end
	end
	return (save != "") ? nothing : cube
end

function assign_description(names::Vector{String}, description::Vector{String}, sentinel2::Int=0)
	# Create a description for each band. If 'description', the cutcube() kwarg, is provided we use it as is.
	# Next we try to find if 'names' indicate a Landsat8 origin and if yes we use the known names & frequencies
	# Otherwise we use the file names as descriptors.
	(!isempty(description) && length(names) != length(description)) &&
		error("'description' and 'names' vectors must have the same length")
	desc = (!isempty(description)) ? description : Vector{String}(undef, length(names))

	t = splitext(splitdir(names[1])[2])[1]
	if (isempty(description) && (startswith(t, "LC")))		# Have Landsat data
		for k = 1:length(names)
			t = splitext(splitdir(names[k])[2])[1]
			ind = findfirst("_B", t)
			bnd = parse(Int, t[ind[end]+1:end])
			desc[k] = Lsat8_desc[bnd]
		end
	elseif (sentinel2 == 10 || sentinel2 == 20 || sentinel2 == 60)
		for k = 1:length(names)
			t = splitdir(names[k])[2]
			ind = findfirst("_B", t)
			bnd = tryparse(Int, t[ind[end]+1:ind[end]+2])
			(bnd === nothing && t[ind[end]+1:ind[end]+2] == "8A") && (bnd = 8)	# ESA is crazzy. Fck character
			desc[k] = (sentinel2 == 10) ? Sentinel2_10m_desc[bnd] : Sentinel2_20m_desc[bnd]	# 20m applyies also for 60m
		end
	else
		if (isempty(description))
			[desc[k] = splitext(splitdir(names[k])[2])[1] for k = 1:length(names)]
		else
			desc = description
		end
	end
	return desc
end

# ----------------------------------------------------------------------------------------------------------
"""
read_mtl(band_name::String, mtl::String=""; get_full=false)

Use the `band_name` of a Landsat8 band to find the MTL file with the scene parameters at which that band
belongs and read the params needed to compute Brightness temperature, radiance at top of atmosphere, etc.
If the MTL file does not lieve next to the band file, send its name via the `mtl` argument.

The `get_full` option makes this function return a tring with contents of the MTL file or `nothing` if
the MTL file is not found.

### Returns a tuple with:

(band=band, rad_mul=rad_mul, rad_add=rad_add, rad_max=rad_max, reflect_mul=reflect_mul, reflect_add=reflect_add, reflect_max=reflect_max, sun_azim=sun_azim, sun_elev=sun_elev, sun_dis=sun_azim, K1=K1, K2=K2)

or a string with MTL contents (or nothing if MTL file is not found)
"""
function read_mtl(fname::String, mtl::String=""; get_full::Bool=false)
	pato,_fname = splitdir(fname)
	((ind = findfirst("_B", _fname)) === nothing && !get_full) &&
		error("This $(fname) is not a valid Landsat8 band file name or of a Landsat 8 cube file.")
	(ind === nothing) && return nothing		# Only happens when get_full = true and name is no Landsat
	if (mtl == "")
		mtl = joinpath(pato, _fname[1:ind[1]] * "MTL.txt")
		if (!isfile(mtl))
			lst = filter(x -> endswith(x, "MTL.txt"), readdir((pato == "") ? "." : pato))
			if (length(lst) == 1 && startswith(lst[1], _fname[1:16]))
				mtl = joinpath(pato, lst[1])
				(!isfile(mtl) && get_full) && return nothing		# Not fatal error in this case
				!isfile(mtl) && (@warn("MTL file was not transmitted in input and I couldn't find it."); return nothing)
			end
		end
	end
	!isfile(mtl) && return nothing			# Not fatal error (if we can live without a MTL).

	(!get_full) && (__fname = splitext(_fname)[1];	band = parse(Int, __fname[ind[1]+2:end]))
	f = open(mtl);	lines = readlines(f);	close(f)
	return (get_full) ? lines : parse_mtl(lines, band)
end

function parse_mtl(mtl::Vector{<:AbstractString}, band::Int)
	# Parse the 'mtl' string vector and extract the info relevant for 'band'
	# Make this a separate function so it can be called from read_mtl() or from the 
	# MTL metadata stored in a cube file created with cutcube(). 

	function get_par(str)
		ind = findfirst(findfirst.(str, mtl) .!== nothing)
		parse(Float64, split(mtl[ind], "=")[2])
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
	MTL_short(band, rad_mul, rad_add, rad_max, reflect_mul, reflect_add, reflect_max, sun_azim, sun_elev, sun_dist, K1, K2)
end

# ----------------------------------------------------------------------------------------------------------
function helper1_sats(fname::String, band_layer::Int)
	I::GMT.GMTimage{UInt16, 2} = (band_layer == 0) ? gmtread(fname) : gmtread(fname, layer=band_layer, layout="TRB")
	indNaN = isnodata(I)
	o = Matrix{Float32}(undef, size(I))
	return I, indNaN, o
end

function parse_lsat8_file(fname::String; band::Int=0, layer::Int=0, mtl::String="")
	# See if 'fname' is of a plain Landsat8 .tif file or of a cube created with cutcube().
	# Depending on the case find the MTL info from file or from the cube's Metadata. Former case
	# still accepts that the MTL file name be transmitted via the 'mtl' option.
	# If no errors so far, parse the MTL and extract the parameters concerning the wished 'band'.
	# This band number will be fetch from the band file name (full Landsat8 product name), or must be
	# transmitted via 'band' option when reading a cube file.
	# In alternative to 'band' we can use 'layer' and fetch the Band number that is stored in that layer.
	GMT.ressurectGDAL()				# Another black-hole plug attempt.
	ds = GMT.Gdal.unsafe_read(fname)
	nbands = GMT.Gdal.nraster(ds)
	meta = GMT.Gdal.GDALGetMetadata(ds.ptr, C_NULL)
	if (nbands > 1)
		desc = Vector{String}(undef,0)
		[append!(desc, [GMT.Gdal.GDALGetDescription(GMT.Gdal.GDALGetRasterBand(ds.ptr, bd))]) for bd = 1:nbands]
	else
		desc = [""]
	end
	GMT.Gdal.GDALClose(ds.ptr)

	if (nbands > 1)
		if (0 < layer < 12)			# Try to find which Band is at layer 'layer'
			((_band = tryparse(Int, split(desc[layer])[1][5:end])) === nothing) && @warn("Failed to find Band with description in layer $layer")
			(_band !== nothing) && (band = _band)
		end
		(band < 1) && error("The `band` option must contain the wished Landsat 8 band number. Use `band=N` to set it.")
		(band > 11) && throw(ArgumentError("Bad Landsat 8 band number $band"))
		MTL = ((ind = findfirst(startswith.(meta, "MTL=GROUP ="))) !== nothing) ? meta[ind][5:end] :
			error("Data in a cube (> one layer) must contain the MTL info in Metadata and this one does not.")

		((band_layer = findfirst(startswith.(desc, "Band$band"))) === nothing) && error("Band $band not found in this cube")
		pars = parse_mtl(split(MTL, "\n"), band)
	else
		pars = read_mtl(fname, mtl)		# No 'band' here because it's suposed to be findable in 'fname'
		band_layer = 0
	end
	band_layer, pars
end

# ----------------------------------------------------------------------------------------------------------
function helper_dns(fname::String, band::Int, bandname::String, bandnames::String, mtl::String)
	# Helper function for the dn2temperature, dn2radiance, dn2... functions
	(band == 0 && bandname == "" && bandnames == "") && error("No information provided on what Band to process.")
	(bandname == "" && bandnames != "") && (bandname = bandnames)	# Accept singular and plural name
	if (bandname != "")
		layers, = find_layers(fname, bandnames=[bandname])
		band_layer, pars = parse_lsat8_file(fname, layer=layers[1], mtl=mtl)
		(band_layer != layers[1]) && error("Internal error. $band_layer and $(layers[1]) should be equal.")
	else
		band_layer, pars = parse_lsat8_file(fname, band=band, mtl=mtl)
	end
	band_layer, pars
end

const dns_doc = "
- `fname`: The name of either a ``LANDSAT_PRODUCT_ID`` geotiff band, or the name of a cube file created with
  the `cutcube` function. In the first case, if the companion ``...MTL.txt`` file is not in the same directory
  as `fname` one can still pass it via the `mtl=path-to-MTL-file` option. In the second case it is mandatory
  to use one of the following two options.
- `band`: _cubes_ created with [`cutcube`](@ref) assign descriptions starting with \"Band 1 ...\" an so on
  the other bands. So when `band` is used we search for the band named \"Band N\", where N = `band`.
- `bandname`: When we know the common designation of a band, for example \"Green\", or any part of a band
  description, for example \"NIR\", we can use that info to create a `bandname` string that will be
  matched against the cube's bands descriptions. We can use the `reportbands` function to see the bands description.
- `save`:  The file name where to save the output. If not provided, a GMTgrid is returned.

Returns a Float32 GMTgrid
"

# ----------------------------------------------------------------------------------------------------------
"""
    R = dn2temperature(fname::String; band::Int=0, mtl::String="", save::String)

Computes the brigthness temperature of Landasat8 termal band (10 or 11)

$dns_doc

# Example:
Compute the brightness temperature of Band 10 stored in a `cube`

```
T = dn2temperature(cube, band=10)
```
"""
function dn2temperature(fname::String; band::Int=0, mtl::String="", save::String="")
	band_layer, pars = parse_lsat8_file(fname, band=band, mtl=mtl)
	(pars.band < 10) && throw(ArgumentError("Brightness temperature is only for bands 10 or 11. Not this one: $(pars.band)"))
	I, indNaN, o = helper1_sats(fname, band_layer)
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
		o[k] = indNaN[k] ? NaN32 : pars.K2 / (log(pars.K1 / (I.image[k] * pars.rad_mul + pars.rad_add) + 1.0)) - 273.15
	end
	G = mat2grid(o, I)
	(save != "") && gdaltranslate(G, dest=save)
	return (save != "") ? nothing : G
end

# ----------------------------------------------------------------------------------------------------------
function helper_dns_op(I::GMTimage, mul::Float64, add::Float64, indNaN::Matrix{Bool}, o::Matrix{Float32})
	# Helper function for the dn2radiance, dn2reflectance functions. Avoid IF branch if no NaNs
	if any(indNaN)
		@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
			o[k] = indNaN[k] ? NaN32 : I.image[k] * mul + add
		end
	else
		@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)  o[k] = I.image[k] * mul + add  end
	end
	mat2grid(o, I)
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = dn2radiance(fname::String, [band::Int, bandname::String, mtl::String, save::String])

Computes the radiance at TopOfAtmosphere of a Landsat 8 file

$dns_doc

# Example:
Compute the radiance TOA of Band 2 file.
```
R = dn2temperature("LC08_L1TP_204033_20210525_20210529_02_T1_B2.TIF")
```
"""
function dn2radiance(fname::String; band::Int=0, bandname::String="", bandnames::String="", mtl::String="", save::String="")
	dn2aux(fname, "radiance"; band=band, bandname=bandname, bandnames=bandnames, mtl=mtl, save=save)
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = dn2reflectance(fname::String, [band::Int, bandname::String, mtl::String, save::String])

Computes the TopOfAtmosphere planetary reflectance of a Landsat8 file

$dns_doc

# Example:
Compute the reflectance TOA of Red Band stored in a `cube`
```
R = dn2reflectance(cube, bandname="red")
```
"""
function dn2reflectance(fname::String; band::Int=0, bandname::String="", bandnames::String="", mtl::String="", save::String="")
	dn2aux(fname, "reflect"; band=band, bandname=bandname, bandnames=bandnames, mtl=mtl, save=save)
end

function dn2aux(fname::String, fun::String; band::Int=0, bandname::String="", bandnames::String="", mtl::String="", save::String="")
	# Most of dn2reflectance & dn2radiance codes are similar, so gather it here under a common function.
	function inn_helper_rad(fname, band, bandname, bandnames, mtl)
		band_layer, pars = helper_dns(fname, band, bandname, bandnames, mtl)
		I, indNaN, o = helper1_sats(fname, band_layer)
		helper_dns_op(I, pars.rad_mul, pars.rad_add, indNaN, o)
	end
	function inn_helper_ref(fname, band, bandname, bandnames, mtl)
		band_layer, pars = helper_dns(fname, band, bandname, bandnames, mtl)
		(pars.band >= 10) && error("Computing Reflectance for Thermal bands is not defined.")
		I, indNaN, o = helper1_sats(fname, band_layer)
		s_elev = sin(pars.sun_elev * pi/180)
		fact_x = pars.reflect_mul / s_elev
		fact_a = pars.reflect_add / s_elev
		helper_dns_op(I, fact_x, fact_a, indNaN, o)
	end

	helper_fun = (fun == "radiance") ? inn_helper_rad : inn_helper_ref		# Which helper function to use.
	if (band == 0 && bandname == "" && bandnames == "")
		bdnames = reportbands(fname)
		if (fun == "reflect")
			((ind = findfirst(contains.(bdnames, "Band 10"))) !== nothing) && (deleteat!(bdnames, ind))
			((ind = findfirst(contains.(bdnames, "Band 11"))) !== nothing) && (deleteat!(bdnames, ind))
		end
		G = helper_fun(fname, 0, bdnames[1], "", mtl)
		mat = G.z
		for k = 2:length(bdnames)
			t = helper_fun(fname, 0, bdnames[k], "", mtl)
			mat = cat(mat, t.z, dims=3)
		end
		G = mat2grid(mat, G)
		G.names = bdnames
	else
		G = helper_fun(fname, band, bandname, bandnames, mtl)
	end
	(save != "") && gdaltranslate(G, dest=save)
	return (save != "") ? nothing : G
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = reflectance_surf(fname::String, [band::Int, bandname::String, mtl::String, save::String])

Computes the radiance-at-surface of Landsat8 band using the COST model.

$dns_doc
"""
function reflectance_surf(fname::String; band::Int=0, mtl::String="", save::String="")
	band_layer, pars = parse_lsat8_file(fname, band=band, mtl=mtl)
	(pars.band >= 10) && error("Computing Surface Reflectance for Thermal bands is not defined.")
	I, indNaN, o = helper1_sats(fname, band_layer)

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
		indNaN[k] && (o[k] = NaN32)
	end
	G = mat2grid(o, I)
	(save != "") && gdaltranslate(G, dest=save)
	return (save != "") ? nothing : G
end
