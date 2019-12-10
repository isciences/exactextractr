# Copyright (c) 2018-2019 ISciences, LLC.
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

library(testthat)
library(exactextractr)
context('exact_extract')

test_that("Basic stat functions work", {
  # This test just verifies a successful journey from R
  # to C++ and back. The correctness of the algorithm
  # is tested at the C++ level.
  data <- matrix(1:9, nrow=3, byrow=TRUE)
  rast <- raster::raster(data,
                         xmn=0, xmx=3, ymn=0, ymx=3,
                         crs='+proj=longlat +datum=WGS84')

  square <- sf::st_sfc(
    sf::st_polygon(
      list(
        matrix(
          c(0.5, 0.5,
            2.5, 0.5,
            2.5, 2.5,
            0.5, 2.5,
            0.5, 0.5),
          ncol=2,
          byrow=TRUE))),
    crs=sf::st_crs(rast))

  dat <- exact_extract(rast, square)

  # Calling without a function returns a data frame with values and coverage fractions
  expect_equal(dat[[1]],
    data.frame(value=1:9,
               coverage_fraction=c(0.25, 0.5, 0.25, 0.5, 1, 0.5, 0.25, 0.5, 0.25))
  )

  # Calling with a function(w, v) returns the result of the function
  expect_equal(exact_extract(rast, square, fun=weighted.mean),
               5)

  # Calling with a string computes a named operation from the C++ library
  expect_equal(exact_extract(rast, square, fun='count'), 4)
  expect_equal(exact_extract(rast, square, fun='mean'), 5)
  expect_equal(exact_extract(rast, square, fun='min'), 1)
  expect_equal(exact_extract(rast, square, fun='max'), 9)
  expect_equal(exact_extract(rast, square, fun='mode'), 5)
  expect_equal(exact_extract(rast, square, fun='majority'), 5)
  expect_equal(exact_extract(rast, square, fun='minority'), 1)
  expect_equal(exact_extract(rast, square, fun='variety'), 9)

  # Can also do multiple stats at once
  expect_equal(exact_extract(rast, square, fun=c('min', 'max', 'mode')),
               data.frame(min=1, max=9, mode=5))

  expect_equal(exact_extract(rast, c(square, square), fun=c('min', 'max', 'mode')),
               data.frame(min=c(1, 1), max=c(9, 9), mode=c(5, 5)))
})

test_that('Raster NA values are correctly handled', {
  data <- matrix(1:100, nrow=10, byrow=TRUE)
  data[7:10, 1:4] <- NA # cut out lower-left corner
  rast <- raster::raster(data,
                         xmn=0, xmx=10, ymn=0, ymx=10,
                         crs='+proj=longlat +datum=WGS84')

  # check polygon entirely within NA region
  circ <- sf::st_sfc(sf::st_buffer(sf::st_point(c(2,2)), 0.9),
                     crs=sf::st_crs(rast))

  expect_equal(0, exact_extract(rast, circ, 'count'))
  expect_equal(NA_real_, exact_extract(rast, circ, 'mean'))
  expect_equal(NA_real_, exact_extract(rast, circ, weighted.mean))

  # check polygon partially within NA region
  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(3.5, 3.5,
          4.5, 3.5,
          4.5, 4.5,
          3.5, 4.5,
          3.5, 3.5),
        ncol=2,
        byrow=TRUE))),
    crs=sf::st_crs(rast))

  expect_equal(43.5, exact_extract(rast, square, 'sum'))
  expect_equal(NA_real_, exact_extract(rast, square, weighted.mean))
  expect_equal(58, exact_extract(rast, square, weighted.mean, na.rm=TRUE))
})

test_that('MultiPolygons also work', {
  data <- matrix(1:100, nrow=10, byrow=TRUE)
  rast <- raster::raster(data,
                         xmn=0, xmx=10, ymn=0, ymx=10,
                         crs='+proj=longlat +datum=WGS84')

  multipoly <- sf::st_sfc(
    sf::st_multipolygon(list(
      sf::st_polygon(
        list(
          matrix(
            c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
            ncol=2,
            byrow=TRUE))),
      sf::st_polygon(
        list(
          matrix(
            4 + c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
            ncol=2,
            byrow=TRUE))))),
    crs=sf::st_crs(rast))

  expect_equal(exact_extract(rast, multipoly, fun='variety'), 18)
})

test_that('We ignore portions of the polygon that extend outside the raster', {
  rast <- raster::raster(matrix(1:(360*720), nrow=360),
                         xmn=-180, xmx=180, ymn=-90, ymx=90,
                         crs='+proj=longlat +datum=WGS84')

  square <- sf::st_sfc(
    sf::st_polygon(
      list(
        matrix(
          c(179.5, 0,
            180.5, 0,
            180.5, 1,
            179.5, 1,
            179.5, 0),
          ncol=2,
          byrow=TRUE))),
    crs=sf::st_crs(rast))

  cells_included <- exact_extract(rast, square, include_xy=TRUE)[[1]][, c('x', 'y')]

  expect_equal(cells_included,
               data.frame(x=179.75, y=c(0.75, 0.25)),
               check.attributes=FALSE)
})

test_that('Additional arguments can be passed to fun', {
  data <- matrix(1:9, nrow=3, byrow=TRUE)
  rast <- raster::raster(data,
                         xmn=0, xmx=3, ymn=0, ymx=3,
                         crs='+proj=longlat +datum=WGS84')

  square <- sf::st_sfc(
    sf::st_polygon(
      list(
        matrix(
          c(0.5, 0.5,
            2.5, 0.5,
            2.5, 2.5,
            0.5, 2.5,
            0.5, 0.5),
          ncol=2,
          byrow=TRUE))),
    crs=sf::st_crs(rast))

  exact_extract(rast, square, function(x, w, custom) {
    expect_equal(custom, 6)
  }, progress=FALSE, 6)
})

test_that('Incorrect argument types are handled gracefully', {
  data <- matrix(1:9, nrow=3, byrow=TRUE)
  rast <- raster::raster(data,
                         xmn=0, xmx=3, ymn=0, ymx=3,
                         crs='+proj=longlat +datum=WGS84')

  point <- sf::st_sfc(sf::st_point(1:2),
                      crs=sf::st_crs(rast))

  linestring <- sf::st_sfc(sf::st_linestring(matrix(1:4, nrow=2)),
                           crs=sf::st_crs(rast))

  multipoint <- sf::st_sfc(sf::st_multipoint(matrix(1:4, nrow=2)),
                           crs=sf::st_crs(rast))

  multilinestring <- sf::st_sfc(
    sf::st_multilinestring(list(
      matrix(1:4, nrow=2),
      matrix(5:8, nrow=2)
    )),
    crs=sf::st_crs(rast))

  geometrycollection <- sf::st_sfc(
    sf::st_geometrycollection(list(
      sf::st_geometry(point)[[1]],
      sf::st_geometry(linestring)[[1]])),
    crs=sf::st_crs(rast))

  expect_error(exact_extract(rast, point), 'unable to find.* method')
  expect_error(exact_extract(rast, linestring, 'unable to find.* method'))
  expect_error(exact_extract(rast, multipoint, 'unable to find.* method'))
  expect_error(exact_extract(rast, multilinesetring, 'unable to find.* method'))
  expect_error(exact_extract(rast, geometrycollection, 'unable to find.* method'))
})

test_that('We can extract values from a RasterStack', {
  rast <- raster::raster(matrix(1:16, nrow=4, byrow=TRUE),
                         xmn=0, xmx=4, ymn=0, ymx=4,
                         crs='+proj=longlat +datum=WGS84')

  stk <- raster::stack(rast, sqrt(rast))

  square <- sf::st_sfc(
    sf::st_polygon(
      list(
        matrix(
          c(0.5, 0.5,
            2.5, 0.5,
            2.5, 2.5,
            0.5, 2.5,
            0.5, 0.5),
          ncol=2,
          byrow=TRUE))),
    crs=sf::st_crs(rast))

  extracted <- exact_extract(stk, square)[[1]]

  expect_equal(names(extracted),
               c('layer.1', 'layer.2', 'coverage_fraction'))
  expect_equal(extracted[, 'layer.2'],
               sqrt(extracted[, 'layer.1']))
  expect_equal(extracted[extracted$coverage_fraction==0.25, 'layer.1'],
               c(5, 7, 13, 15))
  expect_equal(extracted[extracted$coverage_fraction==0.50, 'layer.1'],
               c(6, 9, 11, 14))
  expect_equal(extracted[extracted$coverage_fraction==1.00, 'layer.1'],
               10)
})

test_that('We can pass extracted RasterStack values to an R function', {
  population <- raster::raster(matrix(1:16, nrow=4, byrow=TRUE),
                               xmn=0, xmx=4, ymn=0, ymx=4,
                               crs='+proj=longlat +datum=WGS84')
  income <- sqrt(population)

  square <- sf::st_sfc(
    sf::st_polygon(
      list(
        matrix(
          c(0.5, 0.5,
            2.5, 0.5,
            2.5, 2.5,
            0.5, 2.5,
            0.5, 0.5),
          ncol=2,
          byrow=TRUE))),
    crs=sf::st_crs(population))

  mean_income <- exact_extract(raster::stack(list(population=population, income=income)), square, function(vals, weights) {
    weighted.mean(vals[, 'population']*vals[, 'income'], weights)
  })

  expect_equal(mean_income, 32.64279, tolerance=1e-5)
})

test_that('We get an error when trying to pass extracted RasterStack values to a C++ function', {
  rast <- raster::raster(matrix(runif(16), nrow=4),
                         xmn=0, xmx=4, ymn=0, ymx=4,
                         crs='+proj=longlat +datum=WGS84')

  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(0.5, 0.5,
          2.5, 0.5,
          2.5, 2.5,
          0.5, 2.5,
          0.5, 0.5),
        ncol=2,
        byrow=TRUE))),
    crs=sf::st_crs(rast))

  expect_error(
    exact_extract(raster::stack(rast, sqrt(rast)), square, 'variety'),
    'only available for single-layer raster'
  )
})

test_that('We get acceptable default values when processing a polygon that does not intersect the raster', {
  rast <- raster::raster(matrix(runif(100), nrow=5),
                         xmn=-180, xmx=180, ymn=-65, ymx=85,
                         crs='+proj=longlat +datum=WGS84') # extent of GPW

  poly <- sf::st_sfc(st_polygon(
    list(
      matrix(
        c(-180, -90,
          180, -90,
          180, -65.5,
          -180, -65.5,
          -180, -90),
        ncol=2,
        byrow=TRUE
      )
    )
  ), crs=sf::st_crs(rast)) # extent of Antarctica in Natural Earth

  # RasterLayer
  expect_equal(list(data.frame(value=numeric(), coverage_fraction=numeric())),
               exact_extract(rast, poly))

  expect_equal(list(data.frame(value=numeric(), x=numeric(), y=numeric(), coverage_fraction=numeric())),
               exact_extract(rast, poly, include_xy=TRUE))

  expect_equal(0, exact_extract(rast, poly, function(x, c) sum(x)))
  expect_equal(0, exact_extract(rast, poly, 'count'))
  expect_equal(0, exact_extract(rast, poly, 'sum'))
  expect_equal(0, exact_extract(rast, poly, 'variety'))
  expect_equal(NA_real_, exact_extract(rast, poly, 'majority'))
  expect_equal(NA_real_, exact_extract(rast, poly, 'minority'))
  expect_equal(NA_real_, exact_extract(rast, poly, 'minority'))
  expect_equal(NA_real_, exact_extract(rast, poly, 'mean'))
  expect_equal(NA_real_, exact_extract(rast, poly, 'min'))
  expect_equal(NA_real_, exact_extract(rast, poly, 'max'))

  # RasterStack
  rast2 <- as.integer(rast)
  raster::dataType(rast2) <- 'INT4S'

  stk <- raster::stack(list(q=rast, xi=rast2, area=raster::area(rast)))

  expect_equal(list(data.frame(q=numeric(), xi=integer(), area=numeric(), coverage_fraction=numeric())),
               exact_extract(stk, poly))

  expect_equal(list(data.frame(q=numeric(), xi=integer(), area=numeric(), x=numeric(), y=numeric(), coverage_fraction=numeric())),
               exact_extract(stk, poly, include_xy=TRUE))

  exact_extract(stk, poly, function(values, cov) {
    expect_equal(values, data.frame(q=numeric(), xi=integer(), area=numeric()))
    expect_equal(cov, numeric())
  })
})

test_that('We can optionally get cell center coordinates included in our output', {
  rast <- raster::raster(matrix(1:100, nrow=10),
                         xmn=0, xmx=10, ymn=0, ymx=10,
                         crs='+proj=longlat +datum=WGS84')

  poly <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(3.5, 4.4, 7.5, 4.5, 7.5, 6.5, 3.5, 6.5, 3.5, 4.4),
        ncol=2,
        byrow=TRUE
      )
    )
  ), crs=sf::st_crs(rast))

  results <- exact_extract(rast, poly, include_xy=TRUE)[[1]]

  # check that correct ranges of X,Y values are output
  expect_equal( c(3.5, 4.5, 5.5, 6.5, 7.5), sort(unique(results[, 'x'])))
  expect_equal( c(4.5, 5.5, 6.5),           sort(unique(results[, 'y'])))

  # check the XY values of an individal cell with a known coverage fraction
  expect_equal( results[results[, 'x']==3.5 & results[,'y']==4.5, 'coverage_fraction'],
                0.2968749999999998,
                tolerance=1e-8,
                check.attributes=FALSE)

  # we can also send the weights to a callback
  exact_extract(rast, st_sf(data.frame(id=1), geom=poly), include_xy=TRUE, fun=function(values, weights) {
    expect_equal(3, ncol(values))
  })
})

test_that('Warning is raised on CRS mismatch', {
  rast <- raster::raster(matrix(1:100, nrow=10),
                         xmn=-180, xmx=180, ymn=-90, ymx=90,
                         crs='+proj=longlat +datum=WGS84')

  poly <- sf::st_buffer(
    sf::st_as_sfc('POINT(442944.5 217528.7)', crs=32145),
    150000)

  expect_warning(exact_extract(rast, poly, weighted.mean, na.rm=TRUE),
                 'transformed from .*32145.* to .*4326')
})

test_that('Warning is raised on undefined CRS', {
  rast <- raster::raster(matrix(1:100, nrow=10),
                         xmn=0, xmx=10, ymn=0, ymx=10)

  poly <- sf::st_buffer(sf::st_as_sfc('POINT(8 4)'), 0.4)

  # neither has a defined CRS
  expect_silent(exact_extract(rast, poly, 'sum'))

  # only raster has defined CRS
  raster::crs(rast) <- '+proj=longlat +datum=WGS84'
  expect_warning(exact_extract(rast, poly, 'sum'),
                 'assuming .* same CRS .* raster')

  # both have defined crs
  sf::st_crs(poly) <- sf::st_crs(rast)
  expect_silent(exact_extract(rast, poly, 'sum'))

  # only polygons have defined crs
  raster::crs(rast) <- NULL
  expect_warning(exact_extract(rast, poly, 'sum'),
                 'assuming .* same CRS .* polygon')
})


test_that('Error is raised if function has unexpected signature', {
  rast <- raster::raster(matrix(1:100, nrow=10),
                         xmn=0, xmx=10, ymn=0, ymx=10,
                         crs='+proj=longlat +datum=WGS84')

  poly <- suppressWarnings(sf::st_buffer(
    sf::st_sfc(
      sf::st_point(c(5,5)),
      crs=sf::st_crs(rast)),
    3))

  for (fun in c(length, sum, median, mean, sd)) {
    expect_error(exact_extract(rast, poly, fun),
                 'function .* not .* of the form')
  }

  expect_silent(exact_extract(rast, poly, weighted.mean))
})

test_that('Error is raised for unknown summary operation', {
  rast <- raster::raster(matrix(1:100, nrow=10),
                         xmn=0, xmx=10, ymn=0, ymx=10,
                         crs='+proj=longlat +datum=WGS84')

  poly <- suppressWarnings(sf::st_buffer(
    sf::st_sfc(
      sf::st_point(c(5,5)),
      crs=sf::st_crs(rast)),
    3))

  expect_error(exact_extract(rast, poly, 'whatimean', 'Unknown stat'))
})

test_that('Error is raised if arguments passed to summary operation', {
  rast <- raster::raster(matrix(1:100, nrow=10),
                         xmn=0, xmx=10, ymn=0, ymx=10,
                         crs='+proj=longlat +datum=WGS84')

  poly <- suppressWarnings(sf::st_buffer(
    sf::st_sfc(
      sf::st_point(c(5,5)),
      crs=sf::st_crs(rast)),
    3))

  expect_error(exact_extract(rast, poly, 'sum', na.rm=TRUE),
               'does not accept additional arguments')
})
