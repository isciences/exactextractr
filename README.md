# exactextractr

[![Build Status](https://gitlab.com/isciences/exactextractr/badges/master/build.svg)](https://gitlab.com/isciences/exactextractr/pipelines)
[![Build Status](https://ci.appveyor.com/api/projects/status/aixqdcq7e065eb2h/branch/master?svg=true)](https://ci.appveyor.com/project/dbaston1/exactextractr/branch/master)
[![coverage report](https://gitlab.com/isciences/exactextractr/badges/master/coverage.svg)](https://isciences.gitlab.io/exactextractr/coverage.html)

`exactextractr` is an R package that quickly and accurately summarizes raster
values over polygonal areas, commonly referred to as _zonal statistics_. Unlike
most zonal statistics implementations, it handles grid cells that are partially
covered by a polygon. Typical performance for real-world applications is orders
of magnitude faster than the
[`raster`](https://cran.r-project.org/web/packages/raster/index.html) package.

![Example Graphic](https://exactextractr.s3.us-east-2.amazonaws.com/brazil_precip.png).

Calculations are performed using the C++
[`exactextract`](https://github.com/isciences/exactextract) tool. Additional
background and a description of the method is available
[here](https://github.com/isciences/exactextract#background).
Full package reference documentation is available
[here](https://isciences.gitlab.io/exactextractr/reference).

### Basic Usage

The package provides an
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
method that operates analogously to the
[`extract`](https://www.rdocumentation.org/packages/raster/topics/extract)
method in the
[`raster`](https://cran.r-project.org/web/packages/raster/index.html) package.
The snippet below demonstrates the use of this function to compute a mean
December temperature for each municipality in Brazil.

```r
require(raster)
require(sf)
require(exactextractr)

# Pull municipal boundaries for Brazil
brazil <- st_as_sf(getData('GADM', country='BRA', level=2))

# Pull gridded temperature data
temp <- getData('worldclim', var='tmean', res=10)[[12]]

# Find the mean temperature for each administrative boundary
brazil$mean_temp <- exact_extract(temp, brazil, weighted.mean, na.rm=TRUE)
```

Because
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
takes into account the fraction of the cell that is covered by the polygon, the
summary function must take two arguments: the value of the raster in each cell
touched by the polygon, and the fraction of that cell area that is covered by
the polygon. This differs from
[`raster::extract`](https://www.rdocumentation.org/packages/raster/topics/extract),
where the summary function takes a single argument (the vector of raster values)
and effectively assumes that the coverage fraction is `1.0`.

An example of a built-in function with the appropriate signature is 
[`weighted.mean`](https://www.rdocumentation.org/packages/stats/topics/weighted.mean).
Some examples of custom summary functions are:

```r
# Number of cells covered by the polygon (raster values are ignored)
exact_extract(rast, poly, function(values, coverage_fraction)
                            sum(coverage_fraction))

# Sum of defined raster values within the polygon, accounting for coverage fraction
exact_extract(rast, poly, function(values, coverage_fraction)
                            sum(values * coverage_fraction, na.rm=TRUE))

# Number of distinct raster values within the polygon (coverage fractions are ignored)
exact_extract(rast, poly, function(values, coverage_fraction)
                            length(unique(values)))

# Number of distinct raster values in cells more than 10% covered by the polygon
exact_extract(rast, poly, function(values, coverage_fraction)
                            length(unique(values[coverage_fraction > 0.1])))
```

### Weighting / Geographic Coordinates

[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
also allows for calculation of summary statistics based on
multiple raster layers, such as a population-weighted temperature. If
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
is called with a 
[`RasterStack`](https://www.rdocumentation.org/packages/raster/topics/Raster-class)
instead of a 
[`RasterLayer`](https://www.rdocumentation.org/packages/raster/topics/Raster-class)
, the
summary function will be provided a matrix of raster values and a vector of
coverage fractions. Each column in the matrix represents values from each layer
in the stack, and the columns are named using the names of the layers in the stack.

One application of this feature is the calculation of zonal statistics on raster
data in geographic coordinates. The previous calculation of mean temperature
across Brazilian municipalities assumed that each raster cell covered the same
area, which is not correct for rasters in geographic coordinates
(latitude/longitude).

We can correct for varying cell areas by creating a two-layer `RasterStack`,
with the first layer containing the temperature in each cell, and the second
layer containing the area of each cell. (The
[`area`](https://www.rdocumentation.org/packages/raster/topics/area)
function from the `raster` package will calculate the cell areas for us.) We
continue to use 
[`weighted.mean`](https://www.rdocumentation.org/packages/stats/topics/weighted.mean)
to compute the mean temperature, but instead of using the coverage fractions as
weights, we use the product of the cell area and the coverage fraction.

```r
stk <- stack(list(temp=temp, area=area(temp))),

brazil$mean_temp_weighted <- 
  exact_extract(stk, brazil, function(values, coverage_frac)
                               weighted.mean(values[, 'temp'],
                                             values[, 'area']*coverage_frac,
                                             na.rm=TRUE))
```

With the relatively small polygons used in this example, the error introduced
by assuming fixed cell areas is neglible. However, for large polygons that 
span a wide range of latitudes, this may not be the case.

### Performance and Accuracy

An example benchmark using the example data is shown below. The mean execution
time for `exactextractr` was 3 seconds, vs 1250 for `raster`. Timing was
obtained from execution on an AWS `t2.medium` instance.

```r
microbenchmark(
  a <- exact_extract(temp, brazil, weighted.mean, na.rm=TRUE),
  b <- extract(temp, brazil, mean, na.rm=TRUE), times=5)
  
# Unit: seconds
#                    expr         min          lq        mean     median          uq        max neval
# a <- exact_extract(...)    2.999407    3.035462    3.055226    3.05491    3.089642    3.09671     5
#       b <- extract(...) 1195.111430 1199.662031 1250.310587 1209.81141 1221.649817 1425.31824     5
  
```
Results from `exactextractr` are more accurate than other methods because raster
pixels that are partially covered by polygons are considered. The significance
of partial coverage increases for polygons that are small or irregularly shaped.
For the 5500 Brazilian municipalities used in the example, the error introduced
by incorrectly handling partial coverage is less than 1% for 88% of
municipalties and reaches a maximum of 9%.

### Dependencies

Installation requires the [GEOS](https://geos.osgeo.org/) geometry processing
library. For best performance, it is recommended to use version 3.7, which
introduced some optimizations important to `exactextractr`. On Windows, GEOS
will be downloaded automatically as part of package install. On MacOS, it can be
installed using Homebrew (`brew install geos`). On Linux, it can be installed
from system package repositories (`apt-get install libgeos-dev` on
Debian/Ubuntu, or `yum install libgeos-devel` on CentOS/RedHat.)

### Limitations

 * Raster and polygon inputs must be in the same coordinate reference system.
 * The portion of the raster that intersects any given polygon must fit into 
   memory. (The entire raster does not need to fit into memory.)
