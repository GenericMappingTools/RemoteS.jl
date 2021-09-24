# Exploring a _cube_ of Landsat 8 data
Here we will show examples of the type of things we can do with Landsat 8 data stored in a _cube_ for better organisation and easyness of access.

The data was downloaded from [EarthExplorer](https://earthexplorer.usgs.gov) and comprises the scene with Product ID ```LC08_L1TP_204033_20210525_20210529_02_T1```.

The _cube_ was made with this instructions (but they would only work if you had the scene files in your computer).

```path = C:/v/LC08_L1TP_204033_20210525_20210529_02_T1/LC08_L1TP_204033_20210525_20210529_02_T1_B;```

```cube = cutcube(bands=[2,3,4,5,6,7,10], template=pato, region=(485490,531060,4283280,4330290), save=\"LC08_L1TP_20210525_02_cube.tiff\")```

This creates a 3D GeoTIFF file with the companion MTL file saved in it as Metadata. We can see the band info by running the ```reportbands``` function. That information is quite handy because we can, for example, just refer to the _red_ band and it will figure out which layer of the cube contains the Red band.


```julia
using RemoteS, GMT
reportbands("c:/v/LC08_L1TP_20210525_02_cube.tiff")
7-element Vector{String}:
 "Band 2 - Blue [0.45-0.51]"
 "Band 3 - Green [0.53-0.59]"
 "Band 4 - Red [0.64-0.67]"
 "Band 5 - NIR [0.85-0.88]"
 "Band 6 - SWIR 1 [1.57-1.65]"
 "Band 7 - SWIR 2 [2.11-2.29]"
 "Band 10 - Thermal IR 1 [10.6-11.19]"
```

So to start our exploration the best is to generate a true color image. The ```truecolor``` function knows how to do that automatically.

```julia
Irgb = truecolor("c:/v/LC08_L1TP_20210525_02_cube.tiff");
imshow(Irgb)
```

```@raw html
<img src="LC08_L1TP_20210525_02_RGB.png" width="500" class="center"/>
```

And we can compute the brightness temperature in Celsius at the Top of Atmosphere (TOA) from Band 10.

```julia
T = dn2temperature("c:/v/LC08_L1TP_20210525_02_cube.tiff", band=10);
imshow(T, dpi=150, colorbar=true)
```

```@raw html
<img src="LC08_L1TP_20210525_02_Tbright.png" width="500" class="center"/>
```

Or the Radiance TOA for the Blue band.

```julia
Btoa = dn2radiance("c:/v/LC08_L1TP_20210525_02_cube.tiff", bandname="blue");
imshow(Btoa, dpi=150, color=:gray)
```

```@raw html
<img src="LC08_L1TP_20210525_02_Gtoa.png" width="500" class="center"/>
```

Hmmm, very dark. Let's look at its histogram.

```julia
histogram(Gtoa, show=1)
```

```@raw html
<img src="LC08_L1TP_20210525_02_Gtoa_histo.png" width="350" class="center"/>
```

Yes, the data is very concentraded in the low numbers. We will need to apply a contrast stretch.

To be continued.