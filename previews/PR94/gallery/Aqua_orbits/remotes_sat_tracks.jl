### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ a0a65df0-1e4f-11ec-0c0c-87a3e3c07e73
using RemoteS, GMT, SatelliteToolbox

# The AQUA orbits TLE contents for the month of September 2021 
tle1 = "1 27424U 02022A   21245.83760660  .00000135  00000-0  39999-4 0  9997";
tle2 = "2 27424  98.2123 186.0654 0002229  67.6025 313.3829 14.57107527 28342";
orb = sat_tracks(tle=[tle1; tle2], duration="1D");

# The orbit track can be visualized with
imshow(orb,  proj=:Robinson, region=:global, coast=true)

# ╔═╡ 491c9bd0-1e63-11ec-0eca-33e88b5d9b74
using Dates

orb = sat_tracks(tle=[tle1; tle2], start=DateTime("2021-09-02T13:30:00"), stop=DateTime("2021-09-02T13:40:00"), step="1m");

Dsc = RemoteS.sat_scenes(orb, RemoteS.SCENE_HALFW["AQUA"], "AQUA");

# ╔═╡ 8295ae60-1e4f-11ec-1e88-35a4e3c1ee58
md"# Plot AQUA satellite tracks

Compute 1 day of AQUA satellite orbits starting at current local time. Note, this
will be accurate for the month of September 2021. For other dates it needs an updated TLE. Note also that repeating the commands below will produce different results since the current local time is used as the starting point.

Orbits calculation rely on the [SatelliteToolbox](https://github.com/JuliaSpace/SatelliteToolbox.jl) package, which is not a direct dependency (meaning, it's not loaded automatically) but need to be installed for this to work. Orbits are calculated with the help of the so called Two Line Element files that unfortunatelly have _accuracy validity_ quite short (around one month). They can be obtainded from the [Space Track](https://www.space-track.org) site.
"

# ╔═╡ ff8c1ec0-1e51-11ec-1605-f542f9898cfd
md"There is no point in plotting orbits for a lengthier duration because they would clutter the figure but we can compute them and display over a specific region on Earth.
"

# ╔═╡ 536ac320-1e52-11ec-343e-05e3a18e5b93
# Compute orbits during 2 days at 15 seconds intervals
orb = sat_tracks(tle=[tle1; tle2], duration="2D", step=15);

# and plot them on my "vicinity"
D = RemoteS.within_BB(orb, [-20, 10, 30, 50]);

imshow(D, proj=:guess, coast=true, dpi=200)

# ╔═╡ bbb879e0-1e57-11ec-32f6-ad3bc0a15aac
md"And imagine we would like to know the names of the AQUA scene files that cover a certain location? We can do it with the ``findscenes`` function to which we must provide the location of interest, the satellite (currently only AQUA or TERRA), the duration of the search and if we want Chlorophyl-a concentration or Sea Surface Temperature (see function help for details).
"

# ╔═╡ 63d520b0-1e58-11ec-23b8-c586b28ab99d
# Get the scene names of data with Chlorophyl-a concentration covering
# the point (-8,36) and for two days before "2021-09-07T17:00:00"
tle1 = "1 27424U 02022A   21245.83760660  .00000135  00000-0  39999-4 0  9997";
tle2 = "2 27424  98.2123 186.0654 0002229  67.6025 313.3829 14.57107527 28342";
findscenes(-8,36, start="2021-09-07T17:00:00", sat=:aqua, day=true, duration=-2, oc=1, tle=[tle1, tle2])

# ╔═╡ 6d7213e0-1e5c-11ec-18e4-99322dfd3273
md"
```
	2-element Vector{String}:
	A2021251125500.L2_LAC_OC.nc
	A2021252134000.L2_LAC_OC.nc
```
"

# ╔═╡ 12c0adf0-1e5a-11ec-024e-61169e4c5e09
md"Now you can download them from [```https://oceandata.sci.gsfc.nasa.gov/ob/getfile/A2021251125500.L2_LAC_OC.nc```](https://oceandata.sci.gsfc.nasa.gov/ob/getfile/A2021251125500.L2_LAC_OC.nc ) (but you will need to register first) and play with the ```grid_at_sensor``` function for interpolation and visualization.
"

# ╔═╡ 1dd36980-1e64-11ec-0753-617fe170d8ef
md"Another cool thing we can do is to plot the area extent of the AQUA or TERRA scenes. We do it by first computing and orbit with a `start` and `stop` times and with increments of 1 minute (attention, this is an important factor). With that orbit we use the `sat_tracks` function to compute the polygons delimiting those areas.
"

# ╔═╡ 2d9c2220-1e65-11ec-0b33-1f76ba219cb7
md"The scene names are stored in the `header` field of the `Dsc` GMTdatset

julia> Dsc[1].header

`AQUA_MODIS.20210902T133001.L2.SST.NRT.nc`

julia> Dsc[2].header

`AQUA_MODIS.20210902T133501.L2.SST.NRT.nc`
"

# ╔═╡ Cell order:
# ╠═8295ae60-1e4f-11ec-1e88-35a4e3c1ee58
# ╠═a0a65df0-1e4f-11ec-0c0c-87a3e3c07e73
# ╠═ff8c1ec0-1e51-11ec-1605-f542f9898cfd
# ╠═536ac320-1e52-11ec-343e-05e3a18e5b93
# ╠═bbb879e0-1e57-11ec-32f6-ad3bc0a15aac
# ╠═63d520b0-1e58-11ec-23b8-c586b28ab99d
# ╠═6d7213e0-1e5c-11ec-18e4-99322dfd3273
# ╠═12c0adf0-1e5a-11ec-024e-61169e4c5e09
# ╠═1dd36980-1e64-11ec-0753-617fe170d8ef
# ╠═491c9bd0-1e63-11ec-0eca-33e88b5d9b74
# ╠═2d9c2220-1e65-11ec-0b33-1f76ba219cb7
