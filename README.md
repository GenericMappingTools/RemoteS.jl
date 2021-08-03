# RemoteS.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![codecov.io](http://codecov.io/github/joa-quim/RemoteS.jl/coverage.svg?branch=master)](http://codecov.io/github/joa-quim/RemoteS.jl?branch=master)
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://joa-quim.github.io/RemoteS.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://joa-quim.github.io/RemoteS.jl/dev)
-->

Skeleton package to perform operations in satellite data. It depends heavily, and so far only, on
[GMT.jl](https://github.com/GenericMappingTools/GMT.jl). Since GMT relies on GDAL to read many raster formats
most of satellite and other remote sensing data can be read and processed by this package. 

Install
=======

    ] add https://github.com/GenericMappingTools/RemoteS.jl