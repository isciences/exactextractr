# Copyright (c) 2018-2022 ISciences, LLC.
# All rights reserved.
#
# This software is licensed under the Apache License, Version 2.0 (the "License").
# You may not use this file except in compliance with the License. You may
# obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

context('exact_extract eager loading')

test_that("message emitted when working area doesn't fit in memory", {
  rast_fname <- system.file(file.path('sao_miguel', 'clc2018_v2020_20u1.tif'),
                            package = 'exactextractr')

  poly_fname <- system.file(file.path('sao_miguel', 'concelhos.gpkg'),
                            package = 'exactextractr')

  r <- terra::rast(rast_fname)
  polys <- st_read(poly_fname, quiet = TRUE)

  # no output when everything fits in memory
  capture.output({
    msg <- capture_messages({
      exact_extract(r, polys, 'mode', progress = TRUE, max_cells_in_memory = 1e7)
    })
  })
  expect_equal(msg, character())

  # message emitted when it doesn't fit
  capture.output({
    expect_message(
      exact_extract(r, polys, 'mode', progress = TRUE, max_cells_in_memory = 1e6),
      'Cannot preload'
    )
  })

  # if progress is disabled, so are hints
  expect_silent(
    exact_extract(r, polys, 'mode', progress = FALSE, max_cells_in_memory = 1e6)
  )

  # get additional warning by blowing out the GDAL block cache
  prevCacheSize <- terra::gdalCache()
  terra::gdalCache(1)

  capture.output({
    expect_message(
      exact_extract(r, polys, 'mode', max_cells_in_memory = 1e6),
      'GDAL block size cache is only 1 MB'
    )
  })

  # get additional warning if we are using a RasterStack
  capture.output({
    expect_message(
      exact_extract(stack(r), polys, 'mode', max_cells_in_memory = 1e6),
      'It is recommended to use a SpatRaster'
    )
  })

  terra::gdalCache(prevCacheSize)
})

test_that('cropping does not introduce grid incompatibility', {
  rast_fname <- system.file(file.path('sao_miguel', 'clc2018_v2020_20u1.tif'),
                            package = 'exactextractr')

  poly_fname <- system.file(file.path('sao_miguel', 'concelhos.gpkg'),
                            package = 'exactextractr')

  weight_fname <- system.file(file.path('sao_miguel', 'gpw_v411_2020_density_2020.tif'),
                              package = 'exactextractr')

  r <- terra::rast(rast_fname)
  p <- st_read(poly_fname, quiet = TRUE)
  w <- terra::rast(weight_fname)

  expect_silent({
    exact_extract(r, p, weights = w, grid_compat_tol = 1e-3, progress = FALSE)
  })
})

test_that("eager loading does not change values", {
  # this will fail if terra::crop is not called with snap = 'out'

  rast_fname <- system.file(file.path('sao_miguel', 'clc2018_v2020_20u1.tif'),
                            package = 'exactextractr')

  poly_fname <- system.file(file.path('sao_miguel', 'concelhos.gpkg'),
                            package = 'exactextractr')

  weight_fname <- system.file(file.path('sao_miguel', 'gpw_v411_2020_density_2020.tif'),
                              package = 'exactextractr')

  r <- terra::rast(rast_fname)
  p <- st_read(poly_fname, quiet = TRUE)
  w <- terra::rast(weight_fname)

  no_eager_load <- exact_extract(r, p, weights = w,
                                 include_xy = TRUE,
                                 include_cell = TRUE,
                                 max_cells_in_memory = 2000,
                                 progress = FALSE)

  eager_load <- exact_extract(r, p, weights = w,
                              include_xy = TRUE,
                              include_cell = TRUE,
                              progress = FALSE)

  expect_equal(eager_load, no_eager_load, tol = 2e-7)
})

test_that('eager loading does not error when geometry is outside extent of raster', {
  ras <- terra::rast(matrix(1:100, nrow=10))

  touches_corner <- make_rect(xmin = 10, xmax = 20, ymin = 10, ymax = 20,
                              crs = sf::st_crs(ras))

  loaded <- .eagerLoad(ras, touches_corner, Inf, '')
  expect_equal(
    nrow(exact_extract(loaded, touches_corner)[[1]]),
    0)
})
