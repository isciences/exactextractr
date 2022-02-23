# Copyright (c) 2020-2022 ISciences, LLC.
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

context('exact_resample (terra)')

test_that("exact_resample supports SpatRaster arguments", {
  set.seed(123)

  # generate a random raster with a strange extent and resolution
  src <- raster::raster(matrix(runif(10000), nrow=100),
                               xmn=runif(1), xmx=runif(1) + 9,
                               ymn=runif(1), ymx=runif(1) + 9)

  # resample it to a raster with a larger grid with different resolution
  dst <- raster::raster(xmn=0, xmx=10, ymn=0, ymx=10, res=c(1, 2),
                        crs=raster::crs(src))

  dst <- exact_resample(src, dst, 'sum')

  dst_terra <- exact_resample(terra::rast(src), terra::rast(dst), 'sum')

  expect_equal(terra::rast(dst), dst_terra)
})

test_that("resampling can be weighted with coverage areas instead of coverage fractions", {
  r <- terra::rast(nrows = 10, ncols = 10,
                   xmin = 0, xmax = 10, ymin = 60, ymax = 70, crs = 'EPSG:4326')
  terra::values(r) <- seq(1, terra::ncell(r))

  r2 <- terra::rast(nrows = 1, ncols = 1,
                    xmin = 0, xmax = 10, ymin = 60, ymax = 70, crs = 'EPSG:4326')

  unweighted <- exact_resample(r, r2, 'mean')

  area_weighted <- exact_resample(r, r2, 'mean', coverage_area = TRUE)

  expect_true(area_weighted[1] > unweighted[1])
})

test_that("an R function can be used for resampling", {
  r1 <- make_square_rast(1:100)

  r2 <- terra::rast(nrows = 4, ncols = 4,
                    xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                    crs = terra::crs(r1))

  r2_rfun <- exact_resample(r1, r2, function(value, cov_frac) {
    sum(value * cov_frac)
  })

  r2_stat <- exact_resample(r1, r2, 'sum')

  expect_equal(terra::values(r2_rfun),
               terra::values(r2_stat))
})

test_that("a multi-layer SpatRaster can be provided to an R summary function", {
  r1 <- make_square_rast(1:100, crs = 'EPSG:4326')

  r2 <- terra::rast(nrows = 4, ncols = 4,
                    xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                    crs = terra::crs(r1))

  # calculate an area-weighed mean by putting areas in a second layer
  r1_area <- terra::cellSize(r1)
  r1_stk <- terra::rast(list(r1, r1_area))
  result_a <- exact_resample(r1_stk, r2, function(values, coverage_fraction) {
    weighted.mean(values[,1], values[,2] * coverage_fraction)
  })

  # compare this to the more straightforward method of setting coverage_area = TRUE
  result_b <- exact_resample(r1, r2, 'mean', coverage_area = TRUE)

  expect_equal(
    terra::values(result_a),
    terra::values(result_b),
    tolerance = 1e-3
  )

  expect_error(
    exact_resample(r1_stk, r2, 'mean'),
    'must have a single layer'
  )
})

test_that("error thrown if R function returns non-scalar value", {
  r1 <- make_square_rast(1:100)

  r2 <- terra::rast(nrows = 4, ncols = 4,
                    xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                    crs = terra::crs(r1))

  expect_error(
    exact_resample(r1, r2, function(value, cov_frac) {
      return(1:2)
    }),
    'must return a single value'
  )

  expect_error(
    exact_resample(r1, r2, function(value, cov_frac) {
      return(numeric())
    }),
    'must return a single value'
  )

  expect_error(
    exact_resample(r1, r2, function(value, cov_frac) {
      return(NULL)
    }),
    'Not compatible'
  )

  expect_error(
    exact_resample(r1, r2, function(value, cov_frac) {
      'abc'
    }),
    'Not compatible'
  )
})

test_that("error thrown if R function has wrong signature", {
  r1 <- make_square_rast(1:100)

  r2 <- terra::rast(nrows = 4, ncols = 4,
                    xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                    crs = terra::crs(r1))

  expect_error(
    exact_resample(r1, r2, sum),
    'does not appear to be of the form'
  )
})
