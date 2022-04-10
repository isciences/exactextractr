# exactextractr

[![Build Status](https://gitlab.com/isciences/exactextractr/badges/master/pipeline.svg)](https://gitlab.com/isciences/exactextractr/-/pipelines)
[![Build Status](https://ci.appveyor.com/api/projects/status/aixqdcq7e065eb2h/branch/master?svg=true)](https://ci.appveyor.com/project/dbaston1/exactextractr/branch/master)
[![coverage report](https://gitlab.com/isciences/exactextractr/badges/master/coverage.svg)](https://isciences.gitlab.io/exactextractr/coverage.html)
[![CRAN](http://www.r-pkg.org/badges/version/exactextractr)](https://cran.r-project.org/package=exactextractr)
[![cran checks](https://cranchecks.info/badges/worst/exactextractr)](https://cran.r-project.org/web/checks/check_results_exactextractr.html)

`exactextractr` is an R package that quickly and accurately summarizes raster
values over polygonal areas, commonly referred to as _zonal statistics_. Unlike
most zonal statistics implementations, it handles grid cells that are partially
covered by a polygon. Despite this, it performs faster other packages for many
real-world applications.

![Example Graphic](https://gitlab.com/isciences/exactextractr/-/raw/assets/readme/brazil_precip.png).

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
The snippet below demonstrates the use of this function to compute monthly mean precipitation for each municipality in Brazil.

```r
library(raster)
library(sf)
library(exactextractr)

# Pull municipal boundaries for Brazil
brazil <- st_as_sf(getData('GADM', country='BRA', level=2))

# Pull gridded precipitation data
prec <- getData('worldclim', var='prec', res=10)

# Calculate vector of mean December precipitation amount for each municipality
brazil$mean_dec_prec <- exact_extract(prec[[12]], brazil, 'mean')

# Calculate data frame of min and max precipitation for all months
brazil <- cbind(brazil, exact_extract(prec, brazil, c('min', 'max')))
```

#### Summary Operations

`exactextractr` can summarize raster values using several named operations as well
as arbitrary R functions. Where applicable, a named operation will provide
better performance and reduced memory usage relative to an equivalent R function.
Named operations are specified by providing a character vector with one or more 
operation names to the `fun` parameter of [`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html).

The following summary operations are supported:

| Name                   | Description    |                     
| ---------------------- |--------------- |
| `count`                | Sum of all cell coverage fractions. |
| `majority` (or `mode`) | The raster value with the largest sum of coverage fractions. |
| `max`                  | Maximum value of cells that intersect the polygon, ignoring coverage fractions. |
| `mean`                 | Mean value of cells that intersect the polygon, weighted by the fraction of the cell that is covered. |
| `median`               | Median value of cells that intersect the polygon, weighted by the fraction of the cell that is covered. |
| `quantile`             | Arbitrary quantile value of cells that intersect the polygon, weighted by the fraction of the cell that is covered. |
| `min`                  | Minimum value of cells that intersect the polygon, ignoring coverage fractions. |
| `minority`             | The raster value with the smallest sum of coverage fractions. |
| `sum`                  | Sum of values of raster cells that intersect the polygon, with each raster value weighted by its coverage fraction. |
| `variety`              | The number of distinct raster values in cells wholly or partially covered by the polygon. |
| `variance`             | The population variance of cell values, weighted by the fraction of each cell that is covered by the polygon. |
| `stdev`                | The population standard deviation of cell values, weighted by the fraction of each cell that is covered by the polygon. |
| `coefficient_of_variation` | The population coefficient of variation of cell values, weighted by the fraction of each cell that is covered by the polygon. |
| `frac`                 | Fraction of covered cells that are occupied by each distinct raster value. |

Three additional summary operations require the use of a second weighting raster,
provided in the `weights` argument to 
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html):

| Name                   | Description    |                     
| ---------------------- |--------------- |
| `weighted_mean`        | Mean defined value of cells that intersect the polygon, weighted by the product of the coverage fraction and the value of a second weighting raster. |
| `weighted_sum`         | Sum of defined values of raster cells that intersect the polygon, multiplied by the coverage fraction and the value of a second weighting raster. |
| `weighted_frac`        | Fraction of covered cells that are occupied by each distinct raster value, with coverage fractions multiplied by the value of a second weighting raster. |

Weighted usage is discussed in more detail [below](#weighted-usage).

Undefined (`NA`) values are ignored in all of the named summary operations when 
they occur in the value raster. When they occur in the weighting raster, they 
cause the result of the summary operation to be `NA`.

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

### Weighted Usage

[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
allows for calculation of summary statistics based on
multiple raster layers, such as a population-weighted temperature.
The weighting raster must use the same coordinate system as the primary raster,
and it must use a grid that is compatible with the primary raster. (The resolutions and
extents of the rasters need not be the same, but the higher resolution must must be an 
integer multiple of the lower resolution, and the cell boundaries of both rasters must
coincide with cell boundaries in the higher-resolution grid.)

One application of this feature is the calculation of zonal statistics on
raster data in geographic coordinates. The previous calculation of mean
precipitation amount across Brazilian municipalities assumed that each raster
cell covered the same area, which is not correct for rasters in geographic
coordinates (latitude/longitude).

We can correct for varying cell areas by creating a weighting raster with the area of
each cell in the primary raster using the 
[`area`](https://www.rdocumentation.org/packages/raster/topics/area) function
from the `raster` package.

#### Weighted Summary Operations

Performing a weighted summary with the `weighted_mean` and `weighted_sum` operations
is as simple as providing a weighting
[`RasterLayer`](https://www.rdocumentation.org/packages/raster/topics/Raster-class)
or
[`RasterStack`](https://www.rdocumentation.org/packages/raster/topics/Raster-class)
to the `weights` argument of
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html).

The area-weighted mean precipitation calculation can be expressed as:

```r
brazil$mean_dec_prec_weighted <- exact_extract(prec[[12]], brazil, 'weighted_mean', weights = area(prec))
```

With the relatively small polygons used in this example, the error introduced
by assuming constant cell area is negligible. However, for large polygons that 
span a wide range of latitudes, this may not be the case.

#### Weighted Summary Functions

A weighting raster can also be provided when an R summary function is used.
When a weighting raster is provided, the summary function must accept a third
argument containing the values of the weighting raster.

An equivalent to the `weighted_mean` usage above could be written as:

```r
brazil$mean_dec_prec_weighted <- 
  exact_extract(prec[[12]], brazil, function(values, coverage_frac, weights) {
    weighted.mean(values, coverage_frac * weights)
  }, weights = area(prec))
```

Or, to calculate the area-weighted mean precipitation for all months:

```r
brazil <- cbind(brazil,
  exact_extract(prec, brazil, function(values, coverage_frac, weights) {
    weighted.mean(values, coverage_frac * weights)
  },
  weights = area(prec),
  stack_apply = TRUE))
```

In this example, the `stack_apply` argument is set to `TRUE` so that the summary function
will be applied to each layer of `prec` independently. (If `stack_apply = FALSE`,
the summary function will be called with all values of `prec` in a 12-column
data frame.)

### Additional Usages

#### Multi-Raster Summary Functions

A multi-raster summary function can also be written to implement complex
behavior that requires that multiple layers in a
[`RasterStack`](https://www.rdocumentation.org/packages/raster/topics/Raster-class)
be considered simultaneously.

Here, we compute an area-weighted average temperature by calling
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
with a
[`RasterStack`](https://www.rdocumentation.org/packages/raster/topics/Raster-class)
of minimum and maximum temperatures, and a
[`RasterLayer`](https://www.rdocumentation.org/packages/raster/topics/Raster-class),
of cell areas.

```r
tmin <- getData('worldclim', var = 'tmin', res = 10)
tmax <- getData('worldclim', var = 'tmax', res = 10)

temp <- stack(tmin[[12]], tmax[[12]])

brazil$tavg_dec <- exact_extract(temp, brazil,
  function(values, coverage_fraction, weights) {
    tavg <- 0.5*(values$tmin12 + values$tmax12)
    weighted.mean(tavg, coverage_fraction * weights)
  }, weights = area(prec))
```

When 
[`exact_extract`](https://isciences.gitlab.io/exactextractr/reference/exact_extract.html)
is called with a 
[`RasterStack`](https://www.rdocumentation.org/packages/raster/topics/Raster-class)
of values or weights and `stack_apply = FALSE` (the default), the values or
weights from each layer of the 
[`RasterStack`](https://www.rdocumentation.org/packages/raster/topics/Raster-class)
will be provided to the summary function as a data frame.

In the example above, the summary function is provided with a data frame of
values (containing the values for each layer in the `temp` stack), a vector of
coverage fractions, and a vector of weights.

#### Multi-Valued Summary Functions

In some cases, it is desirable for a summary function to return multiple values
for each input feature. A common application is to summarize the fraction of
each polygon that is covered by a given class of a categorical raster.
This can be accomplished by writing a summary function that returns a one-row
data frame for each input feature. The data frames for each feature will be 
combined into a single data frame using using `rbind` or, if it is available, `dplyr::bind_rows`.

In this example, the mean temperature for each municipality is returned for
each altitude category.

```r
altitude <- getData('alt', country = 'BRA')

prec_for_altitude <- exact_extract(prec[[12]], brazil, function(prec, frac, alt) {
  # ignore cells with unknown altitude
  prec <- prec[!is.na(alt)]
  frac <- frac[!is.na(alt)]
  alt <- alt[!is.na(alt)]
  
  low <- !is.na(alt) & alt < 500
  high <- !is.na(alt) & alt >= 500

  data.frame(
    prec_low_alt = weighted.mean(prec[low], frac[low]),
    prec_high_alt = weighted.mean(prec[high], frac[high])
  )
}, weights = altitude)
```

### Rasterization

`exactextractr` can rasterize polygons though computation of the coverage
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

### Performance

For typical applications, `exactextractr` is much faster than the `raster`
package and somewhat faster than the `terra` package. An example benchmark
is below:

```r
brazil <- st_as_sf(getData('GADM', country='BRA', level=1))
brazil_spat <- as(brazil, 'SpatVector')

prec_rast <- getData('worldclim', var='prec', res=10)
prec_terra <- rast(prec_rast) 
prec12_rast <- prec_rast[[12]]
prec12_terra <- rast(prec_rast[[12]])

microbenchmark(
  extract(prec_rast, brazil, mean, na.rm = TRUE),
  extract(prec_terra, brazil_spat, mean, na.rm = TRUE),
  exact_extract(prec_rast, brazil, 'mean', progress = FALSE),
  exact_extract(prec_terra, brazil, 'mean', progress = FALSE),
  extract(prec12_rast, brazil, mean, na.rm = TRUE),
  extract(prec12_terra, brazil_spat, mean, na.rm = TRUE),
  exact_extract(prec12_rast, brazil, 'mean', progress = FALSE),
  exact_extract(prec12_terra, brazil, 'mean', progress = FALSE),
  times = 5)

```
| Package       | Raster Type | Layers | Expression                                                     | Time (ms)    |                     
| ------------  | ----------- | ------ |---------------------------------------- |--------------------- |
| raster        | RasterLayer | 1      |          `extract(prec_rast, brazil, mean, na.rm = TRUE)`| 48708           | 
| terra         | SpatRaster  | 1      | `extract(prec_terra, brazil_spat, mean, na.rm = TRUE)`|   436           | 
| exactextractr | RasterLayer | 1      | `exact_extract(prec_rast, brazil, "mean", progress = FALSE)`|  1541           | 
| exactextractr | SpatRaster  | 1      | `exact_extract(prec_terra, brazil, "mean", progress = FALSE)`|   129           | 
| raster        | RasterStack | 12     |            `extract(prec12_rast, brazil, mean, na.rm = TRUE)`| 10148           | 
| terra         | SpatRaster  | 12     |      `extract(prec12_terra, brazil_spat, mean, na.rm = TRUE)`|   266           | 
| exactextractr | RasterLayer | 12     |`exact_extract(prec12_rast, brazil, "mean", progress = FALSE)`|   222           | 
| exactextractr | SpatRaster  | 12     |`exact_extract(prec12_terra, brazil, "mean", progress = FALSE)`|   112           | 

Actual performance is a complex topic that can vary dramatically depending on
factors such as:

- the number of layers in the input raster(s)
- the data type of input rasters (for best performance, use a `terra::SpatRaster`)
- the raster file format (GeoTIFF, netCDF, etc)
- the chunking strategy used by the raster file (striped, tiled, etc.)
- the relative size of the area to be read and the GDAL block cache

If `exact_extract` is called with `progress = TRUE`, messages will be emitted
if the package detects a situation that could lead to poor performance, such
as a raster chunk size that is too large to allow caching of blocks between
vector features.

If performance is poor, it may be possible to improve performance by:

- increasing the `max_cells_in_memory` parameter
- increasing the size of the GDAL block cache
- rewriting the input rasters to use a different chunking scheme
- processing inputs as batches of nearby polygons


# Accuracy

Results from `exactextractr` are more accurate than other common
implementations because raster pixels that are partially covered by polygons
are considered.  The significance of partial coverage increases for polygons
that are small or irregularly shaped. For the 5500 Brazilian municipalities
used in the example, the error introduced by incorrectly handling partial
coverage is less than 1% for 88% of municipalities and reaches a maximum of 9%.

### Dependencies

Installation requires version 3.5 or greater of the
[GEOS](https://libgeos.org/) geometry processing library.  It is recommended
to use the most recent released version for best performance. On Windows,
GEOS will be downloaded automatically as part of package install. On MacOS, it
can be installed using Homebrew (`brew install geos`). On Linux, it can be
installed from system package repositories (`apt-get install libgeos-dev` on
Debian/Ubuntu, or `yum install libgeos-devel` on CentOS/RedHat.)
