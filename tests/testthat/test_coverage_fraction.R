# Copyright (c) 2018-2020 ISciences, LLC.
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

context('coverage_fraction')

test_that("Coverage fraction function works", {
  # This test just verifies a successful journey from R
  # to C++ and back. The correctness of the algorithm
  # is tested at the C++ level.
  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
        ncol=2,
        byrow=TRUE))))

  rast <- raster::raster(xmn=0, xmx=3, ymn=0, ymx=3, nrows=3, ncols=3, crs=NA)

  weights <- coverage_fraction(rast, square)[[1]]

  expect_equal(as.matrix(weights),
               rbind(
                 c(0.25, 0.5, 0.25),
                 c(0.50, 1.0, 0.50),
                 c(0.25, 0.5, 0.25)
               ), check.attributes=FALSE)
})

test_that("Output can be cropped to the extent of the input feature", {
  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
        ncol=2,
        byrow=TRUE))))

  rast <- raster::raster(xmn=0, xmx=10, ymn=0, ymx=10, nrows=10, ncols=10, crs=NA)

  weights <- coverage_fraction(rast, square, crop=TRUE)[[1]]

  expect_equal(raster::res(weights), raster::res(rast))
  expect_equal(raster::crs(weights), raster::crs(rast))
  expect_equal(raster::extent(weights), raster::extent(0, 3, 0, 3))
})

test_that("When output is not cropped, cells outside of the processed area are 0, not NA", {
  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
        ncol=2,
        byrow=TRUE))))

  rast <- raster::raster(xmn=0, xmx=10, ymn=0, ymx=10, nrows=10, ncols=10, crs=NA)

  weights <- coverage_fraction(rast, square, crop=TRUE)[[1]]

  expect_false(any(is.na(as.matrix(weights))))
})

test_that('Raster returned by coverage_fraction has same properties as the input', {
  r <- raster::raster(xmn=391030, xmx=419780, ymn=5520000, ymx=5547400, crs=NA)
  raster::res(r) = c(100, 100)
  raster::values(r) <- 1:ncell(r)

  p <- sf::st_as_sfc('POLYGON((397199.680921053 5541748.05921053,402813.496710526 5543125.03289474,407103.299342105 5537246.41447368,398470.733552632 5533962.86184211,397199.680921053 5541748.05921053))')

  w <- coverage_fraction(r, p)

  expect_length(w, 1)
  expect_is(w[[1]], 'RasterLayer')

  expect_equal(raster::res(r),    raster::res(w[[1]]))
  expect_equal(raster::extent(r), raster::extent(w[[1]]))
  expect_equal(raster::crs(r),    raster::crs(w[[1]]))
})

test_that('Coverage fractions are exact', {
  r <- raster::raster(xmn=391030, xmx=419780, ymn=5520000, ymx=5547400, crs=NA)
  raster::res(r) = c(100, 100)
  raster::values(r) <- 1:ncell(r)

  p <- sf::st_as_sfc('POLYGON((397199.680921053 5541748.05921053,402813.496710526 5543125.03289474,407103.299342105 5537246.41447368,398470.733552632 5533962.86184211,397199.680921053 5541748.05921053))')

  w <- coverage_fraction(r, p)

  cell_area <- prod(raster::res(w[[1]]))
  ncells <- raster::cellStats(w[[1]], 'sum')

  expect_equal(sf::st_area(sf::st_geometry(p)),
               ncells*cell_area)
})

test_that('Warning is raised on CRS mismatch', {
  rast <- raster::raster(matrix(1:100, nrow=10),
                         xmn=-75, xmx=-70, ymn=41, ymx=46,
                         crs='+proj=longlat +datum=WGS84')

  poly <- sf::st_buffer(
    sf::st_as_sfc('POINT(442944.5 217528.7)', crs=32145),
    150000)

  expect_warning(coverage_fraction(rast, poly),
                 'transformed to raster')
})

test_that('Warning is raised on undefined CRS', {
  rast <- raster::raster(matrix(1:100, nrow=10),
                         xmn=0, xmx=10, ymn=0, ymx=10)

  poly <- sf::st_buffer(sf::st_as_sfc('POINT(8 4)'), 0.4)

  # neither has a defined CRS
  expect_silent(coverage_fraction(rast, poly))

  # only raster has defined CRS
  raster::crs(rast) <- '+proj=longlat +datum=WGS84'
  expect_warning(coverage_fraction(rast, poly),
                 'assuming .* same CRS .* raster')

  # both have defined crs
  sf::st_crs(poly) <- sf::st_crs(rast)
  expect_silent(coverage_fraction(rast, poly))

  # only polygons have defined crs
  raster::crs(rast) <- NULL
  expect_warning(coverage_fraction(rast, poly),
                 'assuming .* same CRS .* polygon')

})

test_that('Z dimension is ignored, if present', {
  # see https://github.com/isciences/exactextractr/issues/26
  polyz <- st_as_sfc('POLYGON Z ((1 1 0, 4 1 0, 4 4 0, 1 1 0))')
  poly <- st_as_sfc('POLYGON ((1 1, 4 1, 4 4, 1 1))')
  values <- raster(matrix(1:25, nrow=5, ncol=5, byrow=TRUE),
                   xmn=0, xmx=5, ymn=0, ymx=5)

  expect_equal(coverage_fraction(values, poly),
               coverage_fraction(values, polyz))
})

