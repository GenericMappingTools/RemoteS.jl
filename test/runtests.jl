using RemoteS, GMT, Dates
using Test, Printf

include("../src/sat_tracks.jl")		# Even without SatelliteToolbox loaded we can still test some functions

truecolor(mat2img(rand(UInt16,128,128)), mat2img(rand(UInt16,128,128)), mat2img(rand(UInt16,128,128)));
@test RemoteS.guess_increment_from_coordvecs([1., 1, 1, 1], [1., 1, 1, 1]) == 1.0
@test RemoteS.helper_find_sds("AA", "xxxxxxxx:AA", findall("\n", @sprintf("aA\nbbbbbnnn\n"))) == "xxx:AA"

truecolor(mat2img(rand(UInt16, 16, 16, 3), noconv=1), [1,2,3]);
truecolor(mat2img(rand(UInt16, 16, 16), noconv=1), mat2img(rand(UInt16, 16, 16), noconv=1), mat2img(rand(UInt16, 16, 16), noconv=1));

@test RemoteS.assign_description(["LC08_B2.TIF", "LC08_B5.TIF"], String[])[1] == "Band 2 - Blue [0.45-0.51]"
@test RemoteS.assign_description(["AC08_B2.TIF", "AC08_B5.TIF"], String[])  == ["AC08_B2", "AC08_B5"]
@test RemoteS.assign_description(["AC08_B2.TIF", "AC08_B5.TIF"], ["B2", "B5"]) == ["B2", "B5"]

@test reportbands(mat2img(rand(UInt16, 4,4,3), names=["Band 1", "Band 2", "Band 3"]), 3)[1] == "Band 3"
@test reportbands(mat2img(rand(UInt16, 4,4,3), names=["Band 1", "Band 2", "Band 3"]), layers=[1,3]) == ["Band 1", "Band 3"]
@test_throws ErrorException("Bad input argument. Must be a GMTimage or a file name. Not Symbol") reportbands(:aa)

#sat_tracks(position=true);
#sat_tracks(position=true, tle=["1 27424U 02022A   21245.83760660  .00000135  00000-0  39999-4 0  9997"; "2 27424  98.2123 186.0654 0002229  67.6025 313.3829 14.57107527 28342"]);

@test_throws ErrorException("Bad input type Symbol. Must be a DateTime, a String or a Tuple(Int)") getitDTime(:val)
@test get_sat_name(Dict(:sat => "TERRA")) == "TERRA"
@test get_MODIS_scene_name(datetime2julian(DateTime("2020-09-20")), "A") == "AQUA_MODIS.20200920T000001.L2.SST.NRT.nc"
@test get_MODIS_scene_name(datetime2julian(DateTime("2020-09-20")), "A", false) == "A2020264000000.L2_LAC_OC.nc"