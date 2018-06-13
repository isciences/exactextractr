# exactextractr

This package provides routines to perform fast and accurate zonal statistics on `raster` and `sf` objects in R using the [`exactextract`](https://github.com/isciences/exactextract) library.

### Installation

This package is not in CRAN and can be installed from source only. Using [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html), the package can be installed with:

    devtools::install_github('isciences/exactextractr')

On Windows, the R build tools ([Rtools](https://cran.r-project.org/bin/windows/Rtools/)) are also necessary.

### Dependencies

This package depends on the [GEOS](https://trac.osgeo.org/geos) library and is fastest with version 3.7 or later. On Windows, GEOS will be downloaded as part of the build process. On Linux, the library and headers should be installed beforehand using `apt-get install libgeos` or similar.



