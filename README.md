# exactextractr

This package provides routines to perform fast and accurate zonal statistics on `raster` and `sf` objects in R using the [`exactextract`](https://github.com/isciences/exactextract) library.

### Example

```r
require(raster)
require(sf)
require(exactextractr)

brazil <- st_as_sf(getData('GADM', country='BRA', level=2))
temp <- getData('worldclim', var='tmax', res=10)[[12]]

brazil$max_dec_temp <- exact_extract(temp, brazil, weighted.mean)

plot(brazil['max_dec_temp'])
```

### Installation

This package is not in CRAN and can be installed from source only. Using [`devtools`](https://CRAN.R-project.org/package=devtools), the package can be installed with:

```r
devtools::install_github('isciences/exactextractr')
```

On Windows, the R build tools ([Rtools](https://cran.r-project.org/bin/windows/Rtools/)) are also necessary.

### Dependencies

This package depends on the [GEOS](https://trac.osgeo.org/geos) library and is fastest with version 3.7 or later. On Windows, GEOS will be downloaded as part of the build process. On Linux, the library and headers should be installed beforehand using `apt-get install libgeos` or similar.


