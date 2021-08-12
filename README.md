# RemoteS.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![codecov.io](http://codecov.io/github/GenericMappingTools/RemoteS.jl/coverage.svg?branch=main)](https://codecov.io/github/GenericMappingTools/RemoteS.jl?branch=main)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://www.generic-mapping-tools.org/RemoteS.jl/dev/)
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://joa-quim.github.io/RemoteS.jl/stable)
-->

Package to perform operations with satellite data. Easy to use in computing true color images with
automatic contrast stretch, many spectral indices and processing of MODIS L2 files. It depends heavily,
and so far only, on [GMT.jl](https://github.com/GenericMappingTools/GMT.jl). Since GMT relies on GDAL to
read many raster formats most of satellite and other remote sensing data can be read and processed by this package. 

Example, a true color image from Landsat8
=========================================

Load the Red, Green & Blue bands from a Landsat8 scene from AWS (takes a *little* time)

```
using RemoteS, GMT

Ir = gmtread("/vsicurl/http://landsat-pds.s3.amazonaws.com/c1/L8/037/034/LC08_L1TP_037034_20160712_20170221_01_T1/LC08_L1TP_037034_20160712_20170221_01_T1_B4.TIF");

Ig = gmtread("/vsicurl/http://landsat-pds.s3.amazonaws.com/c1/L8/037/034/LC08_L1TP_037034_20160712_20170221_01_T1/LC08_L1TP_037034_20160712_20170221_01_T1_B3.TIF");

Ib = gmtread("/vsicurl/http://landsat-pds.s3.amazonaws.com/c1/L8/037/034/LC08_L1TP_037034_20160712_20170221_01_T1/LC08_L1TP_037034_20160712_20170221_01_T1_B2.TIF");

Irgb = truecolor(Ir, Ig, Ib);

imshow(Irgb)
```
<img src="docs/src/figures/truecolor.png" width="600" class="center"/>

Install
=======

    ] add https://github.com/GenericMappingTools/RemoteS.jl
