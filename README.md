# exactextractr

[![Build Status](https://gitlab.com/isciences/exactextractr/badges/master/build.svg)](https://gitlab.com/isciences/exactextractr/pipelines)

This package provides routines to perform fast and accurate zonal statistics on `raster` and `sf` objects in R using the [`exactextract`](https://github.com/isciences/exactextract) library.
Relative to other implementations, such as the `extract` function in the `raster` package, it:

* computes the exact percentage of each raster cell that is covered by a polygon
* avoids the use of computationally expensive raster disaggregation or polygon clipping methods

Additional background and a description of the method is available [here](https://github.com/isciences/exactextract#background).

### Example

```r
require(raster)
require(sf)
require(exactextractr)

# Pull some administrative boundaries for Brazil
brazil <- st_as_sf(getData('GADM', country='BRA', level=2))

# Pull some climate data (max. temperature)
temp <- getData('worldclim', var='tmax', res=10)[[12]]

# Find the mean-max temperature for each administrative boundary
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

### Limitations

 * Raster and polygon inputs must be in the same coordinate reference system.
 * The portion of the raster that intersects any given polygon must fit into memory. (The entire raster does not need to fit into memory.)

