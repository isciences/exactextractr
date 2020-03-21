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
