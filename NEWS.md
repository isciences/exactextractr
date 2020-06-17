# version 0.4.0

- Add named summary operations for variance, standard deviation, coefficient of variation
- Improve performance, most dramatically for named summary operations on RasterStack or RasterBrick inputs.
- Avoid floating point robustness error when computing stats on certain polygons using
  different resolution value and weighting rasters.

# version 0.3.0

- Add include_cell argument to exact_extract (Michael Sumner)
- Fix error thrown by exact_extract when a polygon is partially outside the extent of the value raster but fully within the extent of the weighting raster.
- Avoid error when processing polygons with a Z dimension on GEOS 3.6
- Support sfc_GEOMETRY inputs if all features are polygonal

# version 0.2.1

- Fix incorrect results for some polygons that exactly follow boundaries of raster grid cells
- Update tests for compatibility with sf >= 0.9.0

# version 0.2.0

- Add weighted mean and weighted sum operations to summarize one raster using the values of another (e.g., population-weighted mean temperature)
- Support use of named summary options on RasterStack and RasterBrick inputs (reduces verbosity, improves performance)
- Fix memory leak in CPP_coverage_fraction

# version 0.1.2

- Generate zero-row data frame for polygons not intersecting the raster

# version 0.1.1

- Attempt to revert to static linking against GEOS if dynamic linking fails
- Avoid undefined behavior when calling `exact_extract` on polygon that does not
  intersect raster.

# version 0.1.0

- Initial release to CRAN
