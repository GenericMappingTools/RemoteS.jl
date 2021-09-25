# Spectral indices with Landsat 8 imagery
We will use here the same data _cube_ that was introduced in the
[_images_](https://www.generic-mapping-tools.org/RemoteS.jl/dev/gallery/L8cube_img/remotes_L8_cube_img/#Exploring-a-*cube*-of-Landsat-8-data-1) example.

A very popular spectral indice is the NDVI which is deffined as the a normalized difference between the
Near Infra Red (NIR) and the Red bands. Namely NDVI = (NIR - Red) / (NIR + Red). An heuristic says that
NDVI values ~greater than 0.4 indicate green vegetation as higher the indice (maximum = 1) larger is the green vegetation content.

Since the data _cube_ holds in it the information about each band, computing the NDVI is a trivial opration because we known where each band is in the data _cube_. Let's see it.


```Julia
using RemoteS, GMT
N = ndvi("c:/v/LC08_L1TP_20210525_02_cube.tiff");
imshow(N, colorbar=true)
```

```@raw html
<img src="../LC08_L1TP_20210525_02_NDVI_raw.png" width="500" class="center"/>
```

The spectral indices functions all have a ``threshold`` value that will NaNify all values ``< threshold``.
Below we wipe out all values < 0.4 to show only the _Green stuff_

```Julia
N = ndvi("c:/v/LC08_L1TP_20210525_02_cube.tiff", threshold=0.4);
imshow(N, dpi=150, colorbar=1)
```

```@raw html
<img src="../LC08_L1TP_20210525_02_NDVI_04.png" width="500" class="center"/>
```

_Cool_ but I would like to check that the result is correct.

The spectral indices functions have another option to obtain just a mask where the values are ``> threshold``,
so we can use it to mask out the
[true color image](https://www.generic-mapping-tools.org/RemoteS.jl/dev/gallery/L8cube_img/LC08_L1TP_20210525_02_RGB.png)
and see what we get. 

```Julia
# Compute a mask based on the condition that threshold >= 0.4
mask = ndvi("c:/v/LC08_L1TP_20210525_02_cube.tiff", threshold=0.4, mask=true);
```

To mask out the true color image the best way is to use the ```mask``` as the alpha band.
To make it easier we will recalculate the true color image here.

```Julia
# Recalculate the true color image
Irgb = truecolor("c:/v/LC08_L1TP_20210525_02_cube.tiff");

# Apply the mask
image_alpha!(Irgb, mask);

# And save it to disk so that we can visualize the result
gmtwrite("c:/v/rgb_masked.tiff", Irgb)

imshow("c:/v/rgb_masked.tiff")
```

```@raw html
<img src="../LC08_L1TP_20210525_02_RGB_masked.png" width="500" class="center"/>
```

While we can see that only the green color is present in the above image, it's not easy to see the details.
So let's zoom in but do also another thing. Let us see also what parts of the RGB image were **not** selected.
To see that we will calculate the ``inverse mask``, that is, the mask that retains all values ```< threshold```.
To achieve that we use the ``mask`` option with a negative number.

```Julia
# Compute a mask based on the condition that threshold >= 0.4
mask_inv = ndvi("c:/v/LC08_L1TP_20210525_02_cube.tiff", threshold=0.4, mask=-1);
image_alpha!(Irgb, mask_inv);
gmtwrite("c:/v/rgb_inv_masked.tiff", Irgb)
```

Plot the green vegetation that passed the NDVI threshold test and the other part, side by side.

```Julia
grdimage("c:/v/rgb_masked.tiff", region=(502380,514200,4311630,4321420), figsize=8, frame=:bare)
grdimage!("c:/v/rgb_inv_masked.tiff", figsize=8, projection=:linear, xshift=8, frame=:bare, show=true)
```

```@raw html
<img src="../LC08_L1TP_20210525_02_RGB_masked2.jpg" width="600" class="center"/>
```

As we can see only the greenest part, the one that is probably more irrigated considering that all the
river and cannals margins were retained, was extracted with the condition ```NDVI > 0.4```. Off course,
using different threshold values would lead to slightly different results.
