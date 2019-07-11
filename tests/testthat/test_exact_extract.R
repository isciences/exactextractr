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
  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
        ncol=2,
        byrow=TRUE))))

  data <- matrix(1:9, nrow=3, byrow=TRUE)

  rast <- raster::raster(data, xmn=0, xmx=3, ymn=0, ymx=3)

  dat <- exact_extract(rast, square)

  # Calling without a function returns a matrix of weights and values
  expect_equal(dat[[1]],
    cbind(values=1:9, weights=c(0.25, 0.5, 0.25, 0.5, 1, 0.5, 0.25, 0.5, 0.25))
  )

  # Calling with a function(w, v) returns the result of the function
  expect_equal(exact_extract(rast, square, fun=weighted.mean),
               5)

  # Calling with a string computes a named stat from the C++ library
  expect_equal(exact_extract(rast, square, fun='count'), 4)
  expect_equal(exact_extract(rast, square, fun='mean'), 5)
  expect_equal(exact_extract(rast, square, fun='min'), 1)
  expect_equal(exact_extract(rast, square, fun='max'), 9)
  expect_equal(exact_extract(rast, square, fun='mode'), 5)
  expect_equal(exact_extract(rast, square, fun='minority'), 1)
  expect_equal(exact_extract(rast, square, fun='variety'), 9)

  # Can also do multiple stats at once
  expect_equal(exact_extract(rast, square, fun=c('min', 'max', 'mode')),
               matrix(c(1, 9, 5), nrow=1))
})

test_that('Raster NA values are correctly handled', {
  data <- matrix(1:100, nrow=10, byrow=TRUE)
  data[7:10, 1:4] <- NA # cut out lower-left corner
  rast <- raster::raster(data, xmn=0, xmx=10, ymn=0, ymx=10)

  # check polygon entirely within NA region
  circ <- sf::st_sfc(sf::st_buffer(sf::st_point(c(2,2)), 0.9))

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
        byrow=TRUE))))

  expect_equal(43.5, exact_extract(rast, square, 'sum'))
  expect_equal(NA_real_, exact_extract(rast, square, weighted.mean))
  expect_equal(58, exact_extract(rast, square, weighted.mean, na.rm=TRUE))
})

test_that('MultiPolygons also work', {
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
            byrow=TRUE))))))

  data <- matrix(1:100, nrow=10, byrow=TRUE)

  rast <- raster::raster(data, xmn=0, xmx=10, ymn=0, ymx=10)

  dat <- exact_extract(rast, multipoly)

  expect_equal(exact_extract(rast, multipoly, fun='variety'), 18)
})

test_that('We ignore portions of the polygon that extend outside the raster', {
  rast <- raster::raster(matrix(1:(360*720), nrow=360),
                         xmn=-180,
                         xmx=180,
                         ymn=-90,
                         ymx=90)

  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(179.5, 0,
          180.5, 0,
          180.5, 1,
          179.5, 1,
          179.5, 0),
        ncol=2,
        byrow=TRUE))))

  cells_included <- exact_extract(rast, square, include_xy=TRUE)[[1]][, c('x', 'y')]

  expect_equal(cells_included,
               rbind(c(179.75, 0.75),
                     c(179.75, 0.25)),
               check.attributes=FALSE)
})

test_that('Additional arguments can be passed to fun', {
  data <- matrix(1:9, nrow=3, byrow=TRUE)
  rast <- raster::raster(data, xmn=0, xmx=3, ymn=0, ymx=3)

  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
        ncol=2,
        byrow=TRUE))))

  exact_extract(rast, square, function(x, w, custom) {
    expect_equal(custom, 6)
  }, 6)
})

test_that('Incorrect argument types are handled gracefully', {
  data <- matrix(1:9, nrow=3, byrow=TRUE)
  rast <- raster::raster(data, xmn=0, xmx=3, ymn=0, ymx=3)

  point <- sf::st_sfc(sf::st_point(1:2))
  linestring <- sf::st_sfc(sf::st_linestring(matrix(1:4, nrow=2)))
  multipoint <- sf::st_sfc(sf::st_multipoint(matrix(1:4, nrow=2)))
  multilinestring <- sf::st_sfc(sf::st_multilinestring(list(
    matrix(1:4, nrow=2),
    matrix(5:8, nrow=2)
  )))
  geometrycollection <- sf::st_sfc(sf::st_geometrycollection(list(
    sf::st_geometry(point)[[1]],
    sf::st_geometry(linestring)[[1]])))

  expect_error(exact_extract(rast, point))
  expect_error(exact_extract(rast, linestring))
  expect_error(exact_extract(rast, multipoint))
  expect_error(exact_extract(rast, multilinesetring))
  expect_error(exact_extract(rast, geometrycollection))
})

test_that('We can extract values from a RasterStack', {
  rast <- raster(matrix(1:16, nrow=4, byrow=TRUE), xmn=0, xmx=4, ymn=0, ymx=4)
  stk <- raster::stack(rast, sqrt(rast))

  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
        ncol=2,
        byrow=TRUE))))

  extracted <- exact_extract(stk, square)[[1]]

  expect_equal(dim(extracted), c(9, 3))
  expect_equal(extracted[, 'layer.2'], sqrt(extracted[, 'layer.1']))
})

test_that('We can pass extracted RasterStack values to an R function', {
  population <- raster(matrix(1:16, nrow=4, byrow=TRUE), xmn=0, xmx=4, ymn=0, ymx=4)
  income <- sqrt(population)

  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
        ncol=2,
        byrow=TRUE))))

  mean_income <- exact_extract(raster::stack(list(population=population, income=income)), square, function(vals, weights) {
    weighted.mean(vals[, 'population']*vals[, 'income'], weights)
  })

  expect_equal(mean_income, 32.64279, tolerance=1e-5)
})

test_that('We get an error when trying to pass extracted RasterStack values to a C++ function', {
  rast <- raster(matrix(runif(16), nrow=4), xmn=0, xmx=4, ymn=0, ymx=4)

  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
        ncol=2,
        byrow=TRUE))))

  expect_error(
    exact_extract(raster::stack(rast, sqrt(rast)), square, 'variety'),
    'only available for single-layer raster'
  )
})

test_that('We can optionally get cell center coordinates included in our output', {
  rast <- raster(matrix(1:100, nrow=10), xmn=0, xmx=10, ymn=0, ymx=10)

  poly <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(3.5, 4.4, 7.5, 4.5, 7.5, 6.5, 3.5, 6.5, 3.5, 4.4),
        ncol=2,
        byrow=TRUE
      )
    )
  ))

  results <- exact_extract(rast, poly, include_xy=TRUE)[[1]]

  # check that correct ranges of X,Y values are output
  expect_equal( c(3.5, 4.5, 5.5, 6.5, 7.5), sort(unique(results[, 'x'])))
  expect_equal( c(4.5, 5.5, 6.5),           sort(unique(results[, 'y'])))

  # check the XY values of an individal cell with a known weight
  expect_equal( results[results[, 'x']==3.5 & results[,'y']==4.5, 'weights'],
                0.2968749999999998,
                tolerance=1e-8,
                check.attributes=FALSE)

  # we can also send the weights to a callback
  exact_extract(rast, st_sf(data.frame(id=1), geom=poly), include_xy=TRUE, fun=function(values, weights) {
    expect_equal(3, ncol(values))
  })
})

test_that('Error is thrown on CRS mismatch', {
  rast <- raster::raster(
    matrix(1:100, nrow=10),
    xmn=-180,
    xmx=180,
    ymn=-90,
    ymx=90,
    crs='+proj=longlat +datum=WGS84'
  )

  poly <- sf::st_buffer(
    sf::st_as_sfc('POINT(442944.5 217528.7)', crs=32145),
    100)

  expect_error(exact_extract(rast, poly, weighted.mean, na.rm=TRUE),
               'must be .* same .* system')
})

test_that('Error is raised if function has unexpected signature', {
  rast <- raster::raster(
    matrix(1:100, nrow=10),
    xmn=0,
    xmx=10,
    ymn=0,
    ymx=10
  )

  poly <- sf::st_buffer(
    sf::st_sfc(
      sf::st_point(c(5,5))),
    3)

  for (fun in c(length, sum, median, mean, sd)) {
    expect_error(exact_extract(rast, poly, fun),
                 'function .* not .* of the form')
  }

  expect_silent(exact_extract(rast, poly, weighted.mean))
})
