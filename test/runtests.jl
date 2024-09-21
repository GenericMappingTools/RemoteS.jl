using RemoteS, GMT, Dates
using Test, Statistics, Printf

include("../src/sat_tracks.jl")		# Even without SatelliteToolbox loaded we can still test some functions

truecolor(mat2img(rand(UInt16,128,128)), mat2img(rand(UInt16,128,128)), mat2img(rand(UInt16,128,128)));
@test RemoteS.guess_increment_from_coordvecs([1., 1, 1, 1], [1., 1, 1, 1]) == [1.0, 1.0]
@test RemoteS.helper_find_sds("AA", "xxxxxxxx:AA", findall("\n", @sprintf("aA\nbbbbbnnn\n"))) == "xxx:AA"

truecolor(mat2img(rand(UInt16, 16, 16, 3), noconv=1), [1,2,3]);
truecolor(mat2img(rand(UInt16, 16, 16), noconv=1), mat2img(rand(UInt16, 16, 16), noconv=1), mat2img(rand(UInt16, 16, 16), noconv=1));

@test RemoteS.assign_description(["LC08_B2.TIF", "LC08_B5.TIF"], String[])[1] == "Band2 - Blue [0.45-0.51]"
@test RemoteS.assign_description(["AC08_B2.TIF", "AC08_B5.TIF"], String[])  == ["AC08_B2", "AC08_B5"]
@test RemoteS.assign_description(["AC08_B2.TIF", "AC08_B5.TIF"], ["B2", "B5"]) == ["B2", "B5"]

@test reportbands(mat2img(rand(UInt16, 4,4,3), names=["Band1", "Band2", "Band3"]), 3)[1] == "Band3"
@test reportbands(mat2img(rand(UInt16, 4,4,3), names=["Band1", "Band2", "Band3"]), layers=[1,3]) == ["Band1", "Band3"]

#sat_tracks(position=true);
#sat_tracks(position=true, tle=["1 27424U 02022A   21245.83760660  .00000135  00000-0  39999-4 0  9997"; "2 27424  98.2123 186.0654 0002229  67.6025 313.3829 14.57107527 28342"]);

@test_throws ErrorException("Bad input type Symbol. Must be a DateTime, a String or a Tuple(Int)") getitDTime(:val)
@test get_sat_name(Dict(:sat => "TERRA")) == "TERRA"
@test get_MODIS_scene_name(datetime2julian(DateTime("2020-09-20")), "A") == "AQUA_MODIS.20200920T000001.L2.SST.NRT.nc"
@test get_MODIS_scene_name(datetime2julian(DateTime("2020-09-20")), "A", false) == "A2020264000000.L2_LAC_OC.nc"

	println("	reading cube")
cube = gmtread("LC08_cube.tiff", layout="TRB");
truecolor("LC08_cube.tiff");	println("	truecolor 1")
truecolor(cube, [3,2,1]);		println("	truecolor 2")
truecolor(cube, layers=[3,2,1]); println("	truecolor 3")
	println("	reportbands 1")
@test reportbands("LC08_cube.tiff", 2)[1] == "Band3 - Green [0.53-0.59]"
	println("	reportbands 2")
reportbands(cube);
	println("	subcube")
Irgb = subcube("LC08_cube.tiff", bandnames = ["red", "green", "blue"]);
RemoteS.parse_lsat8_file("LC08_cube.tiff", band=10);
dn2temperature("LC08_cube.tiff", band=10);	println("	dn2temperature")
dn2radiance("LC08_cube.tiff", band=5);		println("	dn2radiance")
dn2reflectance("LC08_cube.tiff", band=5);	println("	dn2reflectance")
reflectance_surf("LC08_cube.tiff", band=5); println("	reflectance_surf")

R = subcube("LC08_cube.tiff", bandnames=["red"]);	println("	subcube R")
G = subcube("LC08_cube.tiff", bandnames=["green"]);	println("	subcube G")
B = subcube("LC08_cube.tiff", bandnames=["blue"]);	println("	subcube B")
NIR = subcube("LC08_cube.tiff", bandnames=["nir"]);	println("	subcube NIR")
SW1 = subcube("LC08_cube.tiff", bandnames=["swir1"]);	println("	subcube SW1")
SW2 = subcube("LC08_cube.tiff", bandnames=["swir2"]);	println("	subcube SW2")

G1 = evi("LC08_cube.tiff"); println("	evi 1")
G2 = evi(B, R, NIR);		println("	evi 2")
G3 = evi(cube,layers=[1,3,4]);		println("	evi 3")
@test G1.range[5:6] == G2.range[5:6]
@test G1.range[5:6] == G3.range[5:6]

G1 = evi2("LC08_cube.tiff"); println("	evi2 1")
G2 = evi2(R, NIR);			println("	evi2 2")
G3 = evi2(cube);			println("	evi2 3")
@test G1.range[5:6] == G2.range[5:6]
@test G1.range[5:6] == G3.range[5:6]

G1 = gli(R,G,B); 			println("	gli 1")
G2 = gli("LC08_cube.tiff"); println("	gli 2")
G3 = gli(cube);				println("	gli 3")
@test G1.range[5:6] == G2.range[5:6]
@test G1.range[5:6] == G3.range[5:6]
G1 = gli(R,G,B, order="rbg"); println("	gli 4")

G1 = tgi(R,G,B); 			println("	tgi 1")
G2 = tgi("LC08_cube.tiff"); println("	tgi 2")
G3 = tgi(cube);				println("	tgi 3")
@test G1.range[5:6] == G2.range[5:6]
@test G1.range[5:6] == G3.range[5:6]

G1 = vari(R,G,B); 			println("	vari 1")
G2 = vari("LC08_cube.tiff"); println("	vari 2")
G3 = vari(cube);			println("	vari 3")
@test G1.range[5:6] == G2.range[5:6]
@test G1.range[5:6] == G3.range[5:6]

G1 = msavi("LC08_cube.tiff"); println("	msavi 1")
G2 = msavi(R, NIR);			println("	msavi 2")
G3 = msavi(cube);			println("	msavi 3")
@test G1.range[5:6] == G2.range[5:6]
@test G1.range[5:6] == G3.range[5:6]

G1 = ndvi("LC08_cube.tiff"); println("	ndvi 1")
G2 = ndvi(R, NIR);			println("	ndvi 2")
G3 = ndvi(cube);			println("	ndvi 3")
@test G1.range[5:6] == G2.range[5:6]
@test G1.range[5:6] == G3.range[5:6]

G1 = satvi("LC08_cube.tiff"); println("	satvi 1")
G2 = satvi(R, SW1, SW2);	println("	satvi 2")
@test G1.range[5:6] == G2.range[5:6]

G1 = savi("LC08_cube.tiff"); println("	savi 1")
G2 = savi(R, NIR);			println("	savi 2")
G3 = savi(cube);			println("	savi 3")
@test G1.range[5:6] == G2.range[5:6]
@test G1.range[5:6] == G3.range[5:6]

G1 = slavi("LC08_cube.tiff"); println("	slavi 1")
G2 = slavi(R, NIR, SW2);	println("	slavi 2")
G3 = slavi(cube);			println("	slavi 3")
@test G1.range[5:6] == G2.range[5:6]
@test G1.range[5:6] == G3.range[5:6]

G1 = mat2grid(rand(Float32, 16,16)); G2 = mat2grid(rand(Float32, 16,16)); G3 = mat2grid(rand(Float32, 16,16));
I = truecolor(G1,G2,G3);
