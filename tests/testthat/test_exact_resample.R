# Copyright (c) 2020 ISciences, LLC.
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

context('exact_resample')

test_that("exact_resample preserves values", {
  set.seed(123)

  # generate a random raster with a strange extent and resolution
  src <- raster::raster(matrix(runif(10000), nrow=100),
                               xmn=runif(1), xmx=runif(1) + 9,
                               ymn=runif(1), ymx=runif(1) + 9)

  # resample it to a raster with a larger grid with different resolution
  dst <- raster::raster(xmn=0, xmx=10, ymn=0, ymx=10, res=c(1, 2),
                        crs=raster::crs(src))

  dst <- exact_resample(src, dst, 'sum')

  # total values should be preserved
  expect_equal(cellStats(src, 'sum'),
               cellStats(dst, 'sum'))

  # resample it to a raster with a larger grid and a smaller resolution
  dst <- raster::raster(xmn=0, xmx=10, ymn=0, ymx=10, res=c(0.01, 0.02),
                        crs=raster::crs(src))

  dst <- exact_resample(src, dst, 'sum')

  # total values should be preserved
  expect_equal(cellStats(src, 'sum'),
               cellStats(dst, 'sum'))
})

test_that("error thrown if multiple stats provided", {
  src <- make_square_raster(1:100)
  dst <- make_square_raster(1:4)

  expect_error(
    exact_resample(src, dst, c('sum', 'mean')),
    'Only a single')
})

test_that("error thrown if rasters have different CRS", {
  src <- make_square_raster(1:100, crs='+init=epsg:4326')
  dst <- make_square_raster(1:100, crs='+init=epsg:4269')

  expect_error(
    exact_resample(src, dst, 'sum'),
    'same CRS')
})

test_that("warning raised if one CRS undefined", {
  a <- make_square_raster(1:100, crs='+init=epsg:4326')
  b <- make_square_raster(1:100, crs=NA)

  expect_warning(
    exact_resample(a, b, 'sum'),
    'No CRS specified for destination'
  )

  expect_warning(
    exact_resample(b, a, 'sum'),
    'No CRS specified for source'
  )
})
