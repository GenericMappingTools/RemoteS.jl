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

const Lsat8_bd_desc = Dict(
	1 => "Band 1 - Coastal aerosol [0.43-0.45]",
	2 => "Band 2 - Blue [0.45-0.51]",
	3 => "Band 3 - Green [0.53-0.59]",
	4 => "Band 4 - Red [0.64-0.67]",
	5 => "Band 5 - Near Infrared [0.85-0.88]",
	6 => "Band 6 - SWIR 1 [1.57-1.65]",
	7 => "Band 7 - SWIR 2 [2.11-2.29]",
	9 => "Band 9 - Cirrus [1.36-1.38]",
	10 => "Band 10 - Thermal IR 1 [10.6-11.19]",
	11 => "Band 11 - Thermal IR 2 [11.50-12.51]")

# ----------------------------------------------------------------------------------------------------------
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
function truecolor(cube::GMT.GMTimage{UInt16, 3}, wavelength)
	img = Array{UInt8}(undef, size(cube,1), size(cube,2), 3)
	layers = find_layers(cube, wavelength, 3)
	_ = mat2img(@view(cube[:,:,layers[1]]), stretch=true, img8=view(img,:,:,1), scale_only=1)
	_ = mat2img(@view(cube[:,:,layers[2]]), stretch=true, img8=view(img,:,:,2), scale_only=1)
	_ = mat2img(@view(cube[:,:,layers[3]]), stretch=true, img8=view(img,:,:,3), scale_only=1)
	Io = mat2img(img, cube);	Io.layout = "TRBa"
	Io
end

# ----------------------------------------------------------------------------------------------------------
function find_layers(cube::GMT.GMTimage{UInt16, 3}, list::Vector{Int}, n_layers::Int)
	# This function is not finished
	if (maximum(list) < 200)		# The list of bands to pass to the caling fun
		(maximum(list) > size(cube,3)) && error("Not enough bands to satisfy the bands list request.")
		bands = list
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
	else
		error("The `cube` object does not have a frequencies (`v` coordinates) vector as required here.")
	end
	(length(bands) != n_layers) && error("Need $(n_layers) bands but got $(length(bands))")
	bands
end

function find_layers(cube::GMT.GMTimage{UInt16, 3}, fun_bnd_names::Vector{String})
	(isempty(cube.names)) && error("The `cube` object does not have a `names` (band names) assigned field as required here.")
	helper_find_layers(cube.names, fun_bnd_names)
end

function find_layers(fname::String; bands::Vector{Int}=Int[], layers::Vector{Int}=Int[], bandnames::Vector{String}=String[], alllayers::Bool=false)
	# 'layers' just return itself after checking that the cube file actually contain that many layers
	# 'bands' search for "Band bands[k]"
	# 'bandnames' search the bands description for the first layer that contains bandnames[k]
	# Returns the numeric layers (1-based) corresponding to the search criteria and the bands description
	(!alllayers && isempty(bands) && isempty(bandnames)) && error("Must use either the 'bands' OR the 'bandnames' option.")
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

	(!isempty(bands)) && (bandnames = ["Band $(bands[k])" for k = 1:length(bands)])		# Create a bandnames vector
	_layers = helper_find_layers(lowercase.(desc), bandnames)
	return _layers, desc
end

helper_find_layers(desc::Vector{String}, band_names::String) = helper_find_layers(desc, [band_names])
function helper_find_layers(desc::Vector{String}, band_names::Vector{String})::Vector{Int}
	# Note, 'desc' is supposed to be in lowercase already.
	bnd_names = lowercase.(band_names)
	bands, n = zeros(Int, length(bnd_names)), 0
	for bnd_name in bnd_names
		((ind = findfirst(contains.(desc, bnd_name))) !== nothing) && (bands[n += 1] = ind)
	end
	(n != length(bnd_names)) && error("Some of the names in $(band_names) are not present in the `cube` names/description.")
	bands
end

# ----------------------------------------------------------------------------------------------------------
"""
    reportbands(in; [bands::Vector{Int}])
or

    reportbands(in, band;)

Report the bands description of the `in` input argument. This can be a GMTimage or a file name (a String) of 
a 'cube' file. Normally one made with the `cutcube` function. When the use conditions of this function are not met,
either a warning or an error message (if too deep to be caught as a warning) will be issued.

- `bands`: When this optional parameter is used, report the description of the bands in the vector `bands`
- `band`: A scalar with a unique band number. Alternative form to `reportbands(in, bands=[band])`

Returns a string vector.
"""
reportbands(in, band::Int) = reportbands(in, bands=[band])
function reportbands(in; bands::Vector{Int}=Int[])
	if (isa(in, GMT.GMTimage))
		isempty(in.names) && (println("This image object does not have a `names` assigned field"); return nothing)
		isempty(bands) && return in.names		# Return them all
		return (isempty(bands)) ? in.names : in.names[bands]
	elseif (isa(in, String))
		layers, desc = find_layers(in, bands=bands, alllayers=isempty(bands))
		return desc[layers]
	else
		error("Bad input argument. Must be a GMTimage or a file name. Not $(typeof(in))")
	end
end

# ----------------------------------------------------------------------------------------------------------
"""
    cutcube(names=String[], bands=Int[], template="", region=nothing, extension=".TIF", description=String[], save="")

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
function cutcube(; names::Vector{String}=String[], bands::AbstractVector=Int[], template::String="", region=nothing, extension::String=".TIF", description::Vector{String}=String[], save::String="")
	(region === nothing) && error("The `region` option cannot be empty.")
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
	for name in names
		!isfile(name) && error("File name $name does not exist. Must stop here.")
	end
	MTL = read_mtl(names[1], get_full=true)
	(MTL !== nothing) && (MTL = ["MTL=" * join(MTL, "\n")])

	# Little parsing of the -R string but does not test if W < E & S < N
	_region::String = isa(region, String) ? region :
	                  @sprintf("%.12g/%.12g/%.12g/%.12g", region[1], region[2], region[3], region[4])
	startswith(_region, "-R") && (_region = _region[3:end])		# Tolerate a region that starts with "-R"
	(length(findall("/", _region)) != 3) && error("Badly formed region string: $_region")

	desc = assign_description(names, description)
	cube = grdcut(names[1], R=_region)
	mat = cube.image
	for k = 2:length(bands)
		B = grdcut(names[k], R=_region)
		mat = cat(mat, B.image, dims=3)
	end
	cube = mat2img(mat, cube, names=desc)
	(save != "") && gdaltranslate(cube, dest=save, meta=MTL)
	return (save != "") ? nothing : cube
end

function assign_description(names::Vector{String}, description::Vector{String})
	# Create a description for each band. If 'description', the cutcube() kwarg, is provided we use it as is.
	# Next we try to find if 'names' indicate a Landsat8 origin and if yes we use the known names & frequencies
	# Otherwise we use the file names as descriptors.
	(!isempty(description) && length(names) != length(description)) &&
		error("'Description' and 'names' vectors must have the same length")
	desc = (!isempty(description)) ? description : Vector{String}(undef, length(names))

	t = splitext(splitdir(names[1])[2])[1]
	if (isempty(description) && (startswith(t, "LC08_") || startswith(t, "LC8")))	# Have Landsat8 data
		for k = 1:length(names)
			t = splitext(splitdir(names[k])[2])[1]
			ind = findfirst("_B", t)
			bnd = parse(Int, t[ind[end]+1:end])
			desc[k] = Lsat8_bd_desc[bnd]
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
	_fname = splitext(fname)[1]
	((ind = findfirst("_B", fname)) === nothing && !get_full) &&
		error("This $(fname) is not a valid Landsat8 band file name or of a Landsat 8 cube file.")
	(ind === nothing) && return nothing		# Only happens when get_full = true and name is no Landsat
	if (mtl == "")
		mtl = fname[1:ind[1]] * "MTL.txt"
		if (!isfile(mtl))
			pato = splitdir(_fname)[1]
			lst = filter(x -> endswith(x, "MTL.txt"), readdir((pato == "") ? "." : pato))
			if (length(lst) == 1 && startswith(lst[1], _fname[1:16]))
				mtl = joinpath(pato, lst[1])
				(!isfile(mtl) && get_full) && return nothing		# Not fatal error in this case
				!isfile(mtl) && error("MTL file was not transmitted in input and I couldn't find it next the band file.")
			end
		end
	end

	(!get_full) && (band = parse(Int, _fname[ind[1]+2:end]))
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
	I::GMT.GMTimage{UInt16, 2} = (band_layer == 0) ? gmtread(fname) : gmtread(fname, layer=band_layer)
	indNaN = isnodata(I)
	o = Matrix{Float32}(undef, size(I))
	return I, indNaN, o
end

function isnodata(array::AbstractArray, val=0)
	nrows, ncols = size(array,1), size(array,2)
	indNaN = fill(false, nrows, ncols)
	@inbounds Threads.@threads for k = 1:nrows * ncols	# 5x faster than: indNaN = (I.image .== 0)
		(array[k] == val) && (indNaN[k] = true)
	end
	indNaN
end

function parse_lsat8_file(fname::String; band::Int=0, mtl::String="")
	# See if 'fname' is of a plain Landsat8 .tif file or of a cube created with cutcube().
	# Depending on the case find the MTL info from file or from the cube's Metadata. Former case
	# still accepts that the MTL file name be transmitted via the 'mtl' option.
	# If no errors so far, parse the MTL and extract the parameters concerning the wished 'band'.
	# This band number will be fetch from the band file name (full Landsat8 product name), or must be
	# transmitted via 'band' option when reading a cube file.
	(band < 0 || band > 11) && throw(ArgumentError("Bad Landsat 8 band number $band"))		# Must still accept 0 here
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
		(band == 0) && error("The `band` option must contain the wished Landsat 8 band number. Use `band=N` to set it.")
		MTL = ((ind = findfirst(startswith.(meta, "MTL=GROUP ="))) !== nothing) ? meta[ind][5:end] :
			error("Data in a cube (> one layer) must contain the MTL info in Metadata and this one does not.")

		((band_layer = findfirst(startswith.(desc, "Band $band"))) === nothing) && error("Band $band not found in this cube")
		pars = parse_mtl(split(MTL, "\n"), band)
	else
		pars = read_mtl(fname, mtl)		# No 'band' here because it's suposed to be findable in 'fname'
		band_layer = 0
	end
	band_layer, pars
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = dn2temperature(fname::String; band::Int=0, mtl::String="")

Returns a GMTgrid with the brigthness temperature of Landasat8 termal band (10 or 11)

Input can be either a file name of a LANDSAT_PRODUCT_ID geotiff band, or the name of a cube file created
with the `cutcube` function. In the first case, if the companion ...MTL.txt file is not in the same directory
as `fname` one can still pass it via the `mtl=path-to-MTL-file` option. In the second case it is mandatory
to use the `band=N` where N is the band number with the data to convert.
"""
function dn2temperature(fname::String; band::Int=0, mtl::String="")
	band_layer, pars = parse_lsat8_file(fname, band=band, mtl=mtl)
	(pars.band < 10) && throw(ArgumentError("Brightness temperature is only for bands 10 or 11. Not this one: $(pars.band)"))
	I, indNaN, o = helper1_sats(fname, band_layer)
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
		o[k] = pars.K2 / (log(pars.K1 / (I.image[k] * pars.rad_mul + pars.rad_add) + 1.0)) - 273.15
	end
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)  indNaN[k] && (o[k] = NaN)  end
	mat2grid(o, I)
end

#=
function dn2temperature(cube::GMT.GMTimage{UInt16, 3}, layer::Int=0, band::Int=0)
	if (layer == 0)
		((ind = findfirst(startswith.(cube.names, "Band $band"))) === nothing) && error("Band $band not found")
		layer = ind
	end
end
=#

# ----------------------------------------------------------------------------------------------------------
"""
    R = dn2radiance(fname::String; band::Int=0, mtl::String="")

Returns a GMTgrid with the radiance at TopOfAtmosphere for the Landsat8 band file `fname`

Input can be either a file name of a LANDSAT_PRODUCT_ID geotiff band, or the name of a cube file created
with the `cutcube` function. In the first case, if the companion ...MTL.txt file is not in the same directory
as `fname` one can still pass it via the `mtl=path-to-MTL-file` option. In the second case it is mandatory
to use the `band=N` where N is the band number with the data to convert.
"""
function dn2radiance(fname::String; band::Int=0, mtl::String="")
	band_layer, pars = parse_lsat8_file(fname, band=band, mtl=mtl)
	I, indNaN, o = helper1_sats(fname, band_layer)
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)
		o[k] = I.image[k] * pars.rad_mul + pars.rad_add
	end
	@inbounds Threads.@threads for k = 1:size(I,1)*size(I,2)  indNaN[k] && (o[k] = NaN)  end
	mat2grid(o, I)
end

# ----------------------------------------------------------------------------------------------------------
"""
    R = dn2reflectance(fname::String; band::Int=0, mtl::String="")

Returns a GMTgrid with the TopOfAtmosphere planetary reflectance for the Landsat8 band file `fname`

Input can be either a file name of a LANDSAT_PRODUCT_ID geotiff band, or the name of a cube file created
with the `cutcube` function. In the first case, if the companion ...MTL.txt file is not in the same directory
as `fname` one can still pass it via the `mtl=path-to-MTL-file` option. In the second case it is mandatory
to use the `band=N` where N is the band number with the data to convert.
"""
function dn2reflectance(fname::String; band::Int=0, mtl::String="")
	band_layer, pars = parse_lsat8_file(fname, band=band, mtl=mtl)
	(pars.band >= 10) && error("Computing Reflectance for Thermal bands is not defined.")
	I, indNaN, o = helper1_sats(fname, band_layer)
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
    R = reflectance_surf(fname::String; band::Int=0, mtl::String="")

Compute the radiance-at-surface of Landsat8 band using the COST model.

Returns a Float32 GMTgrid type

Input can be either a file name of a LANDSAT_PRODUCT_ID geotiff band, or the name of a cube file created
with the `cutcube` function. In the first case, if the companion ...MTL.txt file is not in the same directory
as `fname` one can still pass it via the `mtl=path-to-MTL-file` option. In the second case it is mandatory
to use the `band=N` where N is the band number with the data to convert.
"""
function reflectance_surf(fname::String; band::Int=0, mtl::String="")
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
		indNaN[k] && (o[k] = NaN)
	end
	mat2grid(o, I)
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

Enhanced vegetation index. Huete et al 1990

EVI = G * ((nir - red) / (nir + C1 * red - C2 * blue + Levi));
C1, C2, G, Levi = 6.0, 7.5, 2.5, 1.
"""
evi(blue, red, nir; kw...) = sp_indices(blue, red, nir; index="EVI", kw...)
evi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 3); index="EVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    EVI2 = evi(red, nir; kw...)

Two-band Enhanced vegetation index. Jiang et al 2008

EVI2 = G * ((nir - red) / (nir + 2.4 * red ))
"""
evi2(red, nir; kw...) = sp_indices(red, nir; index="EVI2", kw...)
evi2(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="EVI2", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    GNDVI = gndvi(green, nir; kw...)

green Normalized diff vegetation index: more sensitive to cholorphyll than ndvi. Gitelson, A., and M. Merzlyak

GNDVI = (nir - green) / (nir + green)
"""
gndvi(green, nir; kw...) = sp_indices(green, nir; index="GNDVI", kw...)
gndvi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="GNDVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    MNDWI = mndwi(green, swir2; kw...)

Modified Normalised Difference Water Index. Xu2006

MNDWI = (green-swir2) / (green+swir2)
"""
mndwi(green, swir2; kw...) = sp_indices(swir2, green; index="MNDWI", kw...)
mndwi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="MNDWI", kw...)

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

Modified soil adjusted vegetation index. Qi 1994

MSAVI = nir + 0.5 - (0.5 * sqrt(pow(2.0 * nir + 1.0, 2) - 8.0 * (nir - (2.0 * red))))
"""
msavi(red, nir; kw...) = sp_indices(red, nir; index="MSAVI", kw...)
msavi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="MSAVI", kw...)

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

Compute the NDVI vegetation index. Input can be either the bands file names, or GMTimage objects
with the band's data.

NDVI = (nir - red) / (nir + red)

Returns either a Float32 GMTgrid or a UInt8 GMTimage if the `mask` option is set to true.
"""
ndvi(red, nir; kw...) = sp_indices(red, nir; index="NDVI", kw...)
ndvi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="NDVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDWI = ndwi(green, nir; kw...)

Normalized difference water index. McFeeters 1996. NDWI => (green - nir)/(green + nir)

NDWI = (green - nir)/(green + nir)
"""
ndwi(green, nir; kw...) = sp_indices(nir, green; index="NDWI", kw...)
ndwi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="NDWI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    NDWI2 = ndwi2(nir, swir2; kw...)

Normalized difference water index. Gao 1996, Chen 2005 (also known as Normalized Difference Moisture Index
NDBI and LSWI)

NDWI2 = (nir - swir2)/(nir + swir2)
"""
ndwi2(nir, swir2; kw...) = sp_indices(swir2, nir; index="NDWI2", kw...)
ndwi2(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="NDWI2", kw...)

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

Soil adjusted vegetation index. Huete 1988

SAVI = (nir - red) * (1.0 + L) / (nir + red + L)
"""
savi(red, nir; kw...) = sp_indices(red, nir; index="SAVI", kw...)
savi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 2); index="SAVI", kw...)

# ----------------------------------------------------------------------------------------------------------
"""
    SLAVI = slavi(red, nir, swir2; kw...)

Specific Leaf Area Vegetation Index. Lymburger 2000

SLAVI = nir / (red + swir2)
"""
slavi(red, nir, swir2; kw...) = sp_indices(red, nir, swir2; index="SLAVI", kw...)
slavi(cube::GMT.GMTimage{UInt16, 3}, bnds; kw...) = sp_indices(cube, find_layers(cube, bnds, 3); index="SLAVI", kw...)

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
	# This version recieves the cube and a vector with the bands list and calls the worker with @view
	if (length(bands) == 2)
		sp_indices(@view(cube[:,:,bands[1]]), @view(cube[:,:,bands[2]]); index=index, kw...)
	else
		sp_indices(@view(cube[:,:,bands[1]]), @view(cube[:,:,bands[2]]), @view(cube[:,:,bands[3]]); index=index, kw...)
	end
end

function sp_indices(bnd1, bnd2, bnd3=nothing; index::String="", kwargs...)
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
