# exactextractr

[![Build Status](https://gitlab.com/isciences/exactextractr/badges/master/pipeline.svg)](https://gitlab.com/isciences/exactextractr/pipelines)
[![Build Status](https://ci.appveyor.com/api/projects/status/aixqdcq7e065eb2h/branch/master?svg=true)](https://ci.appveyor.com/project/dbaston1/exactextractr/branch/master)
[![coverage report](https://gitlab.com/isciences/exactextractr/badges/master/coverage.svg)](https://isciences.gitlab.io/exactextractr/coverage.html)
[![CRAN](http://www.r-pkg.org/badges/version/exactextractr)](https://cran.r-project.org/package=exactextractr)
[![cran checks](https://cranchecks.info/badges/worst/exactextractr)](https://cran.r-project.org/web/checks/check_results_exactextractr.html)

`exactextractr` is an R package that quickly and accurately summarizes raster
values over polygonal areas, commonly referred to as _zonal statistics_. Unlike
most zonal statistics implementations, it handles grid cells that are partially
covered by a polygon. Typical performance for real-world applications is orders
of magnitude faster than the
[`raster`](https://CRAN.R-project.org/package=raster) package.

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
[`raster`](https://CRAN.R-project.org/package=raster) package.
The snippet below demonstrates the use of this function to compute a mean
December precipitation for each municipality in Brazil.

```r
library(raster)
library(sf)
library(exactextractr)

# Pull municipal boundaries for Brazil
brazil <- st_as_sf(getData('GADM', country='BRA', level=2))

# Pull gridded precipitation data
prec <- getData('worldclim', var='prec', res=10)[[12]]

# Find the mean precipitation amount for each municipality
brazil$mean_prec <- exact_extract(prec, brazil, 'mean')

# Find min and max precipitation amount in a single pass
brazil[, c('min_prec', 'max_prec')] <- exact_extract(prec, brazil, c('min', 'max'))
```

#### Summary Operations

`exactextractr` can summarize raster values using several pre-defined operations as well
as arbitrary R functions. Pre-defined operations are specified by providing one or more
operation names to the `fun` parameter of
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html).

The following summary operations are supported:

| Name                   | Description    |                     
| ---------------------- |--------------- |
| `count`                | Sum of all cell coverage fractions. |
| `majority` (or `mode`) | The raster value with the largest sum of coverage fractions. |
| `max`                  | Maximum defined value of cells that intersect the polygon, ignoring coverage fractions. |
| `mean`                 | Mean defined value of cells that intersect the polygon, weighted by the percent of the cell that is covered. |
| `min`                  | Minimum defined value of cells that intersect the polygon, ignoring coverage fractions. |
| `minority`             | The raster value with the smallest sum of coverage fractions. |
| `sum`                  | Sum of defined values of raster cells that intersect the polygon, with each raster value weighted by its coverage fraction. |
| `variety`              | The number of distinct raster values in cells wholly or partially covered by the polygon. |

#### Summary Functions

In addition to the summary operations described above,
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
can accept an R function to summarize the cells covered by the polygon. Because
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
takes into account the fraction of the cell that is covered by the polygon, the
summary function must take two arguments: the value of the raster in each cell
touched by the polygon, and the fraction of that cell area that is covered by
the polygon. (This differs from
[`raster::extract`](https://www.rdocumentation.org/packages/raster/topics/extract),
where the summary function takes the vector of raster values as a single argument
and effectively assumes that the coverage fraction is `1.0`.)

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
summary function will be called with a data frame of raster values and a vector
of coverage fractions as arguments. Each column in the data frame represents
values from one layer in the stack, and the columns are named using the names
of the layers in the stack.

One application of this feature is the calculation of zonal statistics on
raster data in geographic coordinates. The previous calculation of mean
precipitation amount across Brazilian municipalities assumed that each raster
cell covered the same area, which is not correct for rasters in geographic
coordinates (latitude/longitude).

We can correct for varying cell areas by creating a two-layer
[`RasterStack`](https://www.rdocumentation.org/packages/raster/topics/Raster-class),
with the first layer containing the precipitation amount in each cell and the second
layer containing the area of each cell. (The
[`area`](https://www.rdocumentation.org/packages/raster/topics/area) function
from the `raster` package will calculate the cell areas for us.) We use
[`weighted.mean`](https://www.rdocumentation.org/packages/stats/topics/weighted.mean)
to compute the mean precipitation amount, using the product of the cell area and the
coverage fraction as a weight.

```r
stk <- stack(list(prec=prec, area=area(prec)))

brazil$mean_prec_weighted <-
  exact_extract(stk, brazil, function(values, coverage_frac)
                               weighted.mean(values$prec, values$area*coverage_frac, na.rm=TRUE))
```

With the relatively small polygons used in this example, the error introduced
by assuming fixed cell areas is negligible. However, for large polygons that 
span a wide range of latitudes, this may not be the case.

### Rasterization

`exactextractr` can also rasterize polygons though computation of the coverage
fraction in each cell. The
[`coverage_fraction`](https://isciences.gitlab.io/exactextractr/reference/coverage_fraction.html)
function returns a
[`RasterLayer`](https://www.rdocumentation.org/packages/raster/topics/Raster-class)
with values from 0 to 1 indicating the fraction of each cell that is covered by
the polygon. Because this function generates a
[`RasterLayer`](https://www.rdocumentation.org/packages/raster/topics/Raster-class)
for each feature in the input dataset, it can quickly consume a large amount of
memory. Depending on the analysis being performed, it may be advisable to
manually loop over the features in the input dataset and combine the generated
rasters during each iteration.

### Performance and Accuracy

An example benchmark using the example data is shown below. The mean execution
time for `exactextractr` was 2.6 seconds, vs 136 for `raster`. Timing was
obtained from execution on an AWS `t2.medium` instance.

```r
microbenchmark(
  a <- exact_extract(prec, brazil, weighted.mean),
  b <- extract(prec, brazil, mean, na.rm=TRUE), times=5)
  
# Unit: seconds
#               expr         min          lq        mean     median          uq        max neval
# a <- exact_extract(...)    2.5674   2.586868   2.626761   2.587283   2.613296   2.778957     5
#       b <- extract(...)  136.1710 136.180563 136.741275 136.226435 136.773627 138.354764     5


Results from `exactextractr` are more accurate than other methods because raster
pixels that are partially covered by polygons are considered. The significance
of partial coverage increases for polygons that are small or irregularly shaped.
For the 5500 Brazilian municipalities used in the example, the error introduced
by incorrectly handling partial coverage is less than 1% for 88% of
municipalities and reaches a maximum of 9%.
```

### Dependencies

Installation requires version 3.5 or greater of the
[GEOS](https://geos.osgeo.org/) geometry processing library.  It is recommended
to use the most recent released version (3.8) for best performance. On Windows,
GEOS will be downloaded automatically as part of package install. On MacOS, it
can be installed using Homebrew (`brew install geos`). On Linux, it can be
installed from system package repositories (`apt-get install libgeos-dev` on
Debian/Ubuntu, or `yum install libgeos-devel` on CentOS/RedHat.)
