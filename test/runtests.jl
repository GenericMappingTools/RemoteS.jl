using RemoteS, GMT
using Test, Printf

truecolor(mat2img(rand(UInt16,128,128)), mat2img(rand(UInt16,128,128)), mat2img(rand(UInt16,128,128)));
@test RemoteS.guess_increment_from_coordvecs([1., 1, 1, 1], [1., 1, 1, 1]) == 1.0
@test RemoteS.helper_find_sds("AA", "xxxxxxxx:AA", findall("\n", @sprintf("aA\nbbbbbnnn\n"))) == "xxx:AA"

truecolor(mat2img(rand(UInt16, 16, 16, 3), noconv=1), [1,2,3]);
truecolor(mat2img(rand(UInt16, 16, 16), noconv=1), mat2img(rand(UInt16, 16, 16), noconv=1), mat2img(rand(UInt16, 16, 16), noconv=1));
