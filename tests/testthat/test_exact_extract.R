# Copyright (c) 2018-2020 ISciences, LLC.
# All rights reserved.
#
# This software is licensed under the Apache License, Version 2.0 (the "License").
# You may not use this file except in compliance with the License. You may
# obtain a copy of the License ta http://www.apache.org/licenses/LICENSE-2.0.
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

  square <- make_rect(0.5, 0.5, 2.5, 2.5, sf::st_crs(rast))

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
  expect_equal(exact_extract(rast, square, fun='median'), 5)
  expect_equal(exact_extract(rast, square, fun='quantile', quantiles=0.25), 3.5)
  expect_equal(exact_extract(rast, square, fun='quantile', quantiles=0.75), 6.5)
  expect_equal(exact_extract(rast, square, fun='min'), 1)
  expect_equal(exact_extract(rast, square, fun='max'), 9)
  expect_equal(exact_extract(rast, square, fun='mode'), 5)
  expect_equal(exact_extract(rast, square, fun='majority'), 5)
  expect_equal(exact_extract(rast, square, fun='minority'), 1)
  expect_equal(exact_extract(rast, square, fun='variety'), 9)
  expect_equal(exact_extract(rast, square, fun='variance'), 5)
  expect_equal(exact_extract(rast, square, fun='stdev'), sqrt(5))
  expect_equal(exact_extract(rast, square, fun='coefficient_of_variation'), sqrt(5)/5)

  # Can also do multiple stats at once
  expect_equal(exact_extract(rast, square, fun=c('min', 'max', 'mode')),
               data.frame(min=1, max=9, mode=5))

  expect_equal(exact_extract(rast, c(square, square), fun=c('min', 'max', 'mode'), progress = FALSE),
               data.frame(min=c(1, 1), max=c(9, 9), mode=c(5, 5)))
})

test_that('Weighted stat functions work', {
  data <- matrix(1:9, nrow=3, byrow=TRUE)
  rast <- raster::raster(data,
                         xmn=0, xmx=3, ymn=0, ymx=3,
                         crs='+proj=longlat +datum=WGS84')

  equal_weights <- raster::raster(matrix(1, nrow=3, ncol=3),
                                  xmn=0, xmx=3, ymn=0, ymx=3,
                                  crs='+proj=longlat +datum=WGS84')

  bottom_row_only <- raster::raster(rbind(c(0, 0, 0), c(0, 0, 0), c(1, 1, 1)),
                                    xmn=0, xmx=3, ymn=0, ymx=3,
                                    crs='+proj=longlat +datum=WGS84')

  square <- make_rect(0.5, 0.5, 2.5, 2.5, sf::st_crs(rast))

  expect_equal(exact_extract(rast, square, 'weighted_mean', weights=equal_weights),
               exact_extract(rast, square, 'mean'))

  expect_equal(exact_extract(rast, square, 'weighted_sum', weights=equal_weights),
               exact_extract(rast, square, 'sum'))

  expect_equal(exact_extract(rast, square, 'weighted_mean', weights=bottom_row_only),
               (0.25*7 + 0.5*8 + 0.25*9)/(0.25 + 0.5 + 0.25))

  expect_equal(exact_extract(rast, square, 'weighted_sum', weights=bottom_row_only),
               (0.25*7 + 0.5*8 + 0.25*9))
})

test_that('Error thrown if weighted stat requested but weights not provided', {
  rast <- make_square_raster(1:9)
  square <- make_circle(2, 2, 0.5, sf::st_crs(rast))

  for (stat in c('weighted_mean', 'weighted_sum')) {
    expect_error(exact_extract(rast, square, stat),
                 'no weights provided')
  }
})

test_that('Warning raised if weights provided but weighted stat not requested', {
  rast <- make_square_raster(1:9)
  square <- make_circle(2, 2, 0.5, sf::st_crs(rast))

  for (stat in c('count', 'sum', 'mean', 'min', 'max', 'minority', 'majority', 'mode', 'variety')) {
    expect_warning(exact_extract(rast, square, stat, weights=rast),
                   'Weights provided but no.*operations use them')
  }
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
  square <- make_rect(3.5, 3.5, 4.5, 4.5, sf::st_crs(rast))

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

test_that('Generic sfc_GEOMETRY works if the features are polygonal', {
  rast <- make_square_raster(1:100)
  polys <- st_as_sfc(c('POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))',
                       'MULTIPOLYGON (((2 2, 4 2, 4 4, 2 4, 2 2)), ((4 4, 8 4, 8 8, 4 8, 4 4)))'),
                     crs=sf::st_crs(rast))

  expect_equal(exact_extract(rast, polys, 'count', progress = FALSE),
               c(4, 4+16))
})

test_that('Generic sfc_GEOMETRY fails if a feature is not polygonal', {
  rast <- make_square_raster(1:100)
  features <- st_as_sfc(c('POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))',
                          'POINT (2 7)'), crs=sf::st_crs(rast))

  expect_error(exact_extract(rast, features, 'sum', progress = FALSE),
               'must be polygonal')
})

test_that('We ignore portions of the polygon that extend outside the raster', {
  rast <- raster::raster(matrix(1:(360*720), nrow=360),
                         xmn=-180, xmx=180, ymn=-90, ymx=90,
                         crs='+proj=longlat +datum=WGS84')

  rect <- make_rect(179.5, 0, 180.5, 1, sf::st_crs(rast))

  cells_included <- exact_extract(rast, rect, include_xy=TRUE)[[1]][, c('x', 'y')]

  expect_equal(cells_included,
               data.frame(x=179.75, y=c(0.75, 0.25)),
               check.attributes=FALSE)


  index_included <- exact_extract(rast, rect, include_xy=TRUE, include_cell = TRUE)[[1]][, c('x', 'y', 'cell')]
  expect_equivalent(as.matrix(cells_included[c("x", "y")]),
               raster::xyFromCell(rast, index_included$cell))
  expect_equal(index_included$cell,
               raster::cellFromXY(rast, cbind(cells_included$x, cells_included$y)))
})

test_that('Additional arguments can be passed to fun', {
  data <- matrix(1:9, nrow=3, byrow=TRUE)
  rast <- raster::raster(data,
                         xmn=0, xmx=3, ymn=0, ymx=3,
                         crs='+proj=longlat +datum=WGS84')

  square <- make_rect(0.5, 0.5, 2.5, 2.5, sf::st_crs(rast))

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

  square <- make_rect(0.5, 0.5, 2.5, 2.5, sf::st_crs(rast))

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

  square <- make_rect(0.5, 0.5, 2.5, 2.5, sf::st_crs(population))

  mean_income <- exact_extract(raster::stack(list(population=population, income=income)), square, function(vals, weights) {
    weighted.mean(vals[, 'population']*vals[, 'income'], weights)
  })

  expect_equal(mean_income, 32.64279, tolerance=1e-5)
})

test_that('We can pass extracted RasterStack values to a C++ function', {
  rast <- raster::raster(matrix(runif(16), nrow=4),
                         xmn=0, xmx=4, ymn=0, ymx=4,
                         crs='+proj=longlat +datum=WGS84')

  square <- make_rect(0.5, 0.5, 2.5, 2.5, sf::st_crs(rast))

  stk <- raster::stack(list(a=rast, b=sqrt(rast)))
  brk <- raster::brick(stk)

  for (input in c(stk, brk)) {
    expect_equal(
      exact_extract(input, square, 'variety'),
      data.frame(variety.a=9, variety.b=9)
    )

    twostats <- exact_extract(input, square, c('variety', 'mean'))

    expect_equal(nrow(twostats), 1)
    expect_named(twostats, c('variety.a', 'variety.b', 'mean.a', 'mean.b'))
  }
})

test_that('We can apply the same function to each layer of a RasterStack', {
  set.seed(123)

  stk <- raster::stack(list(a = make_square_raster(runif(100)),
                            b = make_square_raster(runif(100))))

  circles <- c(
    make_circle(5, 4, 2, sf::st_crs(stk)),
    make_circle(3, 1, 1, sf::st_crs(stk)))

  # by default layers are processed together
  expect_error(
    exact_extract(stk, circles, weighted.mean, progress=FALSE),
    'must have the same length'
  )

  # but we can process them independently with stack_apply
  means <- exact_extract(stk, circles, weighted.mean, progress=FALSE, stack_apply=TRUE)

  expect_named(means, c('fun.a', 'fun.b'))

  # results are same as we would get by processing layers independently
  for (i in 1:raster::nlayers(stk)) {
    expect_equal(means[, i], exact_extract(stk[[i]], circles, weighted.mean, progress=FALSE))
  }
})

test_that('We can use the stack_apply argument with include_xy and include_cols', {
  set.seed(123)

  stk <- raster::stack(list(a = make_square_raster(runif(100)),
                            b = make_square_raster(runif(100))))

  circles <- c(
    make_circle(5, 4, 2, sf::st_crs(stk)),
    make_circle(3, 1, 1, sf::st_crs(stk)))

  result <- exact_extract(stk, circles, include_xy = TRUE, stack_apply = TRUE, progress = FALSE,
                          function(df, frac) {
                            weighted.mean(df$value[df$y > 1],
                                          frac[df$y > 1])
                          })

  expect_named(result, c('fun.a', 'fun.b'))
})

test_that('We can summarize a RasterStack / RasterBrick using weights from a RasterLayer', {
  set.seed(123)

  stk <- raster::stack(list(a = make_square_raster(1:100),
                            b = make_square_raster(101:200)))

  weights <- make_square_raster(runif(100))


  circle <- make_circle(5, 4, 2, sf::st_crs(stk))

  # same weights get used for both
  expect_equal(exact_extract(stk, circle, 'weighted_mean', weights=weights),
               data.frame(weighted_mean.a = 63.0014,
                          weighted_mean.b = 163.0014),
               tolerance=1e-6)

  # error when trying to use a stack as weights
  expect_error(exact_extract(stk, circle, 'weighted_mean', weights=stk),
               "Weighting raster must have only a single layer")

  # error when trying to use a non-raster as weights
  expect_error(exact_extract(stk, circle, 'weighted_mean', weights='stk'),
               "Weights must be a Raster")
})

test_that('We get an error trying to use weights without a named summary operation', {
  rast <- make_square_raster(1:100)
  weights <- make_square_raster(runif(100))
  circle <- make_circle(5, 4, 2, sf::st_crs(rast))

  expect_error(exact_extract(rast, circle, function(value, coverage_fraction) { value }, weights=weights),
               'Weighting raster can only be used with named summary operations')
})

test_that('We get acceptable default values when processing a polygon that does not intersect the raster', {
  rast <- raster::raster(matrix(runif(100), nrow=5),
                         xmn=-180, xmx=180, ymn=-65, ymx=85,
                         crs='+proj=longlat +datum=WGS84') # extent of GPW

  poly <- make_rect(-180, -90, 180, -65.5, sf::st_crs(rast)) # extent of Antarctica in Natural Earth

  # RasterLayer
  expect_equal(list(data.frame(value=numeric(), coverage_fraction=numeric())),
               exact_extract(rast, poly))

  expect_equal(list(data.frame(value=numeric(),
                               x=numeric(),
                               y=numeric(),
                               cell=numeric(),
                               coverage_fraction=numeric())),
               exact_extract(rast, poly, include_xy=TRUE, include_cell=TRUE))

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

  expect_equal(list(data.frame(q=numeric(),
                               xi=integer(),
                               area=numeric(),
                               x=numeric(),
                               y=numeric(),
                               cell=numeric(),
                               coverage_fraction=numeric())),
               exact_extract(stk, poly, include_xy=TRUE, include_cell=TRUE))

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

  results <- exact_extract(rast, poly, include_xy=TRUE, include_cell=TRUE)[[1]]

  # check that correct ranges of X,Y values are output
  expect_equal( c(3.5, 4.5, 5.5, 6.5, 7.5), sort(unique(results[, 'x'])))
  expect_equal( c(4.5, 5.5, 6.5),           sort(unique(results[, 'y'])))

  expect_equal(results[, 'cell'], raster::cellFromXY(rast, results[, c('x', 'y')]))

  # check the XY values of an individual cell with a known coverage fraction
  expect_equal( results[results[, 'x']==3.5 & results[,'y']==4.5, 'coverage_fraction'],
                0.2968749999999998,
                tolerance=1e-8,
                check.attributes=FALSE)

  # we can also send the weights to a callback
  exact_extract(rast, sf::st_sf(data.frame(id=1), geom=poly), include_xy=TRUE, fun=function(values, weights) {
    expect_equal(3, ncol(values))
  }, progress=FALSE)
})

test_that('Warning is raised on CRS mismatch', {
  rast <- raster::raster(matrix(1:100, nrow=10),
                         xmn=-180, xmx=180, ymn=-90, ymx=90,
                         crs='+proj=longlat +datum=WGS84')

  poly <- sf::st_buffer(
    sf::st_as_sfc('POINT(442944.5 217528.7)', crs=32145),
    150000)

  expect_warning(exact_extract(rast, poly, weighted.mean, na.rm=TRUE),
                 'transformed to raster')
})

test_that('Warning is raised on undefined CRS', {
  rast <- raster::raster(matrix(1:100, nrow=10),
                         xmn=0, xmx=10, ymn=0, ymx=10)

  weights <- raster::raster(matrix(runif(100), nrow=10),
                            xmn=0, xmx=10, ymn=0, ymx=10)

  poly <- make_circle(8, 4, 0.4, crs=NA_integer_)

  # neither has a defined CRS
  expect_silent(exact_extract(rast, poly, 'sum'))

  # only raster has defined CRS
  raster::crs(rast) <- '+proj=longlat +datum=WGS84'
  expect_warning(exact_extract(rast, poly, 'sum'),
                 'assuming .* same CRS .* raster')

  # weights have no defined CRS
  expect_warning(exact_extract(rast, poly, 'weighted_mean', weights=weights),
                 'No CRS .* weighting raster.* assuming .* same CRS')

  # both have defined crs
  sf::st_crs(poly) <- sf::st_crs(rast)
  expect_silent(exact_extract(rast, poly, 'sum'))

  # only polygons have defined crs
  raster::crs(rast) <- NULL
  expect_warning(exact_extract(rast, poly, 'sum'),
                 'assuming .* same CRS .* polygon')
})

test_that('Error thrown if value raster and weighting raster have different crs', {
   values <- make_square_raster(runif(100), crs=NA)
   weights <- make_square_raster(runif(100), crs=NA)

   poly <- make_circle(8, 4, 1.5, crs=NA_real_)

   # no CRS for values or weights
   exact_extract(values, poly, 'weighted_mean', weights=weights)

   # values have defined CRS, weights do not
   raster::crs(values) <- '+proj=longlat +datum=WGS84'
   raster::crs(weights) <- '+proj=longlat +datum=NAD83'
   expect_error(
     exact_extract(values, poly, 'weighted_mean', weights=weights),
     'Weighting raster does not have .* same CRS as value raster')
})


test_that('Error thrown if value raster and weighting raster have incompatible grids', {
  poly <- make_circle(5, 4, 2, NA_integer_)

  values <- raster::raster(matrix(runif(10*10), nrow=10),
                           xmn=0, xmx=10, ymn=0, ymx=10)


  # weights have same extent as values, higher resolution
  weights <- raster::raster(matrix(runif(100*100), nrow=100),
                            xmn=0, xmx=10, ymn=0, ymx=10)

  exact_extract(values, poly, 'weighted_mean', weights=weights)

  # weights have same extent as values, lower resolution
  weights <- raster::raster(matrix(1:4, nrow=2),
                            xmn=0, xmx=10, ymn=0, ymx=10)

  exact_extract(values, poly, 'weighted_mean', weights=weights)

  # weights have offset extent from values, same resolution, compatible origin
  weights <- raster::raster(matrix(runif(10*10), nrow=2),
                           xmn=1, xmx=11, ymn=2, ymx=12)

  exact_extract(values, poly, 'weighted_mean', weights=weights)

  # weights have offset extent from values, same resolution, incompatible origin
  weights <- raster::raster(matrix(runif(10*10), nrow=2),
                           xmn=0.5, xmx=10.5, ymn=2, ymx=12)

  expect_error(exact_extract(values, poly, 'weighted_mean', weights=weights),
               'Incompatible extents')
})

test_that('Error is raised if function has unexpected signature', {
  rast <- make_square_raster(1:100)

  poly <- make_circle(5, 5, 3, sf::st_crs(rast))

  for (fun in c(length, sum, median, mean, sd)) {
    expect_error(exact_extract(rast, poly, fun),
                 'function .* not .* of the form')
  }

  expect_silent(exact_extract(rast, poly, weighted.mean))
})

test_that('Error is raised for unknown summary operation', {
  rast <- make_square_raster(1:100)

  poly <- make_circle(5, 5, 3, sf::st_crs(rast))

  expect_error(exact_extract(rast, poly, 'whatimean'),
               'Unknown stat')
})

test_that('Error is raised if arguments passed to summary operation', {
  rast <- make_square_raster(1:100)

  poly <- make_circle(5, 5, 3, sf::st_crs(rast))

  expect_error(exact_extract(rast, poly, 'sum', na.rm=TRUE),
               'does not accept additional arguments')
})

test_that('Error is raised for invalid max_cells_in_memory', {
  rast <- make_square_raster(1:100)
  poly <- make_circle(5, 5, 3, sf::st_crs(rast))

  expect_error(exact_extract(rast, poly, 'mean', max_cells_in_memory=-123),
               'Invalid.*max_cells')
})

test_that('Correct results obtained when max_cells_in_memory is limited', {
  rast <- make_square_raster(1:100)
  poly <- make_circle(5, 5, 3, sf::st_crs(rast))

  expect_equal(exact_extract(rast, poly, 'mean'),
               exact_extract(rast, poly, 'mean', max_cells_in_memory=1))
})

test_that('Weighted stats work when polygon is contained in weight raster but only partially contained in value raster', {
  values <- raster(matrix(1:15, nrow=3, ncol=5, byrow=TRUE),
                   xmn=0, xmx=5, ymn=2, ymx=5)
  weights <- raster(sqrt(matrix(1:25, nrow=5, ncol=5, byrow=TRUE)),
                    xmn=0, xmx=5, ymn=0, ymx=5)
  poly <- make_circle(2.1, 2.1, 1, NA_real_)

  value_tbl <- exact_extract(values, poly, include_xy=TRUE)[[1]]
  weight_tbl <- exact_extract(weights, poly, include_xy=TRUE)[[1]]

  tbl <- merge(value_tbl, weight_tbl, by=c('x', 'y'))

  expect_equal(
    exact_extract(values, poly, 'weighted_mean', weights=weights),
    weighted.mean(tbl$value.x, tbl$coverage_fraction.x * tbl$value.y),
    tol=1e-6
  )
})

test_that('When part of a polygon is within the value raster but not the
           weighting raster, values for unweighted stats requested at the
           same time as weighted stats are correct', {

  values <- raster(matrix(1:25, nrow=5, ncol=5, byrow=TRUE),
                    xmn=0, xmx=5, ymn=0, ymx=5)
  weights <- raster(sqrt(matrix(1:15, nrow=3, ncol=5, byrow=TRUE)),
                   xmn=0, xmx=5, ymn=2, ymx=5)
  poly <- make_circle(2.1, 2.1, 1, NA_real_)

  expect_equal(
    exact_extract(values, poly, 'sum'),
    exact_extract(values, poly, c('sum', 'weighted_mean'), weights=weights)$sum
  )
})

test_that('When polygon is entirely outside the value raster and entirely
           within the weighting raster, we get NA instead of an exception', {
  values <- raster(matrix(1:25, nrow=5, ncol=5, byrow=TRUE),
                    xmn=5, xmx=10, ymn=5, ymx=10)
  weights <- raster(matrix(1:10, nrow=10, ncol=10, byrow=TRUE),
                   xmn=0, xmx=10, ymn=0, ymx=10)
  poly <- make_circle(2.1, 2.1, 1, NA_real_)

  expect_equal(NA_real_,
               exact_extract(values, poly, 'weighted_mean', weights=weights))
})

test_that('Z dimension is ignored, if present', {
  # see https://github.com/isciences/exactextractr/issues/26
  poly <- st_as_sfc('POLYGON Z ((1 1 0, 4 1 0, 4 4 0, 1 1 0))')
  values <- raster(matrix(1:25, nrow=5, ncol=5, byrow=TRUE),
                   xmn=0, xmx=5, ymn=0, ymx=5)

  expect_equal(exact_extract(values, poly, 'sum'), 70.5) # CPP code path
  expect_equal(exact_extract(values, poly, function(x,f) sum(x*f)), 70.5) # R code path
})

test_that('No error thrown when weighting with different resolution grid (regression)', {
   poly <- st_as_sfc(structure(list(
     '01060000000200000001030000000100000008000000065bb0055b7866401c222222223233c0454444242e776640338ee338842d33c0abaaaacac0776640338ee338962733c0676666469f776640a4aaaaaa362033c03a8ee3784f7866404f555555a41c33c0a64ffa840b7966406c1cc771522133c0454444645a796640f4a44ffa9c2b33c0065bb0055b7866401c222222223233c0010300000001000000080000004b9ff4499f7c6640a3aaaaaaaaaa32c0bdbbbb3b747a6640f8ffff7f549632c0ea933e09aa7b664004b6608b399132c0b1055bb0637e6640dc388e63278f32c0d9822d58827e6640dc388ee3109432c09a999979837c6640590bb660159c32c0676666867c7d664070777777039c32c04b9ff4499f7c6640a3aaaaaaaaaa32c0'),
     class='WKB'), EWKB=TRUE)

   v <- raster(matrix(1:360*720, nrow=360, ncol=720),
               xmn=-180, xmx=180, ymn=-90, ymx=90)

   w <- raster(matrix(1:360*720*36, nrow=360*6, ncol=720*6),
               xmn=-180, xmx=180, ymn=-90, ymx=90)

   exact_extract(v, poly, 'weighted_sum', weights=w)
   succeed()
})

test_that('We can get data frame output if we request it', {
  rast <- make_square_raster(1:100)
  names(rast) <- 'z'
  poly <- c(make_circle(5, 5, 3, sf::st_crs(rast)),
            make_circle(3, 1, 1, sf::st_crs(rast)))

  vals <- exact_extract(rast, poly, 'mean', progress=FALSE)

  # named summary operation
  vals_df <- exact_extract(rast, poly, 'mean', force_df=TRUE, progress=FALSE)
  expect_s3_class(vals_df, 'data.frame')
  expect_equal(vals, vals_df[['mean']])

  # R function
  vals2_df <- exact_extract(rast, poly, weighted.mean, force_df=TRUE, progress=FALSE)
  expect_s3_class(vals2_df, 'data.frame')
  expect_equal(vals, vals2_df[['result']], tol=1e-6)
})

test_that('We can have include the input raster name in column names even if
           the input raster has only one layer', {
  rast <- make_square_raster(1:100)
  names(rast) <- 'z'
  poly <- c(make_circle(5, 5, 3, sf::st_crs(rast)),
            make_circle(3, 1, 1, sf::st_crs(rast)))

  vals <- exact_extract(rast, poly, c('mean', 'sum'), progress=FALSE)
  expect_named(vals, c('mean', 'sum'))

  # named summary operations
  vals_named <- exact_extract(rast, poly, c('mean', 'sum'), full_colnames=TRUE, progress=FALSE)
  expect_named(vals_named, c('mean.z', 'sum.z'))
})

test_that('We can summarize a categorical raster by returning a data frame from a custom function', {
  set.seed(456) # smaller circle does not have class 5

  classes <- c(1, 2, 3, 5)

  rast <- raster::raster(xmn = 0, xmx = 10, ymn = 0, ymx = 10, res = 1)
  values(rast) <- sample(classes, length(rast), replace = TRUE)

  circles <- c(
    make_circle(5, 4, 2, sf::st_crs(rast)),
    make_circle(3, 1, 1, sf::st_crs(rast)))

  # approach 1: classes known in advance
  result <- exact_extract(rast, circles, function(x, c) {
    row <- lapply(classes, function(cls) sum(c[x == cls]))
    names(row) <- paste('sum', classes, sep='_')
    do.call(data.frame, row)
  }, progress = FALSE)

  expect_named(result, c('sum_1', 'sum_2', 'sum_3', 'sum_5'))

  # check a single value
  expect_equal(result[2, 'sum_3'],
               exact_extract(rast, circles[2, ], function(x, c) {
                sum(c[x == 3])
               }))

  if (requireNamespace('dplyr', quietly = TRUE)) {
    # approach 2: classes not known in advance (requires dplyr::bind_rows)
    result2 <- exact_extract(rast, circles, function(x, c) {
      found_classes <- unique(x)
      row <- lapply(found_classes, function(cls) sum(c[x == cls]))
      names(row) <- paste('sum', found_classes, sep='_')
      do.call(data.frame, row)
    }, progress = FALSE)

    for (colname in names(result)) {
      expect_equal(result[[colname]], dplyr::coalesce(result2[[colname]], 0))
    }
  }
})

test_that('Error is thrown when using include_* with named summary operation', {
  rast <- make_square_raster(1:100)

  circles <- st_sf(
    fid = c(2, 9),
    size = c('large', 'small'),
    geometry =  c(
    make_circle(5, 4, 2, sf::st_crs(rast)),
    make_circle(3, 1, 1, sf::st_crs(rast))))

  expect_error(exact_extract(rast, circles, 'sum', include_xy = TRUE),
               'include_xy must be FALSE')

  expect_error(exact_extract(rast, circles, 'sum', include_cell = TRUE),
               'include_cell must be FALSE')

  expect_error(exact_extract(rast, circles, 'sum', include_cols = 'fid'),
               'include_cols not supported')
})

test_that('We can append columns from the source data frame in the results', {
  rast <- make_square_raster(1:100)

  circles <- st_sf(
    fid = c(2, 9),
    size = c('large', 'small'),
    geometry =  c(
    make_circle(5, 4, 2, sf::st_crs(rast)),
    make_circle(3, 1, 1, sf::st_crs(rast))))

  result_1 <- exact_extract(rast, circles, 'mean', append_cols = c('size', 'fid'), progress = FALSE)
  expect_named(result_1, c('size', 'fid', 'mean'))

  result_2 <- exact_extract(rast, circles, weighted.mean, append_cols = c('size', 'fid'), progress = FALSE)
  # result_2 won't be identical to result_2 because the column names are different
  # instead, check that the naming is consistent with what we get from the force_df argument
  expect_identical(result_2,
                   cbind(sf::st_drop_geometry(circles[, c('size', 'fid')]),
                         exact_extract(rast, circles, weighted.mean, force_df = TRUE)))
})

test_that('We can include columns from the source data frame in returned data frames', {
  rast <- make_square_raster(1:100)

  circles <- st_sf(
    fid = c(2, 9),
    size = c('large', 'small'),
    geometry =  c(
    make_circle(5, 4, 2, sf::st_crs(rast)),
    make_circle(3, 1, 1, sf::st_crs(rast))))

  combined_result <- do.call(rbind, exact_extract(rast, circles, include_cols = 'fid', progress = FALSE))
  expect_named(combined_result, c('fid', 'value', 'coverage_fraction'))
})

test_that('We can get multiple quantiles with the "quantiles" argument', {
 rast <- make_square_raster(1:100)

  circles <- st_sf(
    fid = c(2, 9),
    size = c('large', 'small'),
    geometry =  c(
    make_circle(5, 4, 2, sf::st_crs(rast)),
    make_circle(3, 1, 1, sf::st_crs(rast))))

  result <- exact_extract(rast, circles, 'quantile', quantiles=c(0.25, 0.50, 0.75), progress=FALSE)
  expect_true(inherits(result, 'data.frame'))

  expect_named(result, c('q25', 'q50', 'q75'))
})

test_that('Error is thrown if quantiles not specified or not valid', {
 rast <- make_square_raster(1:100)
 square <- make_rect(2, 2, 4, 4, crs=sf::st_crs(rast))

 expect_error(exact_extract(rast, square, 'quantile'),
              'Quantiles not specified')

 expect_error(exact_extract(rast, square, 'quantile', quantiles=NA),
              'must be between 0 and 1')

 expect_error(exact_extract(rast, square, 'quantile', quantiles=c(0.5, 1.1)),
              'must be between 0 and 1')

 expect_error(exact_extract(rast, square, 'quantile', quantiles=numeric()),
              'Quantiles not specified')
})

test_that('Progress bar updates incrementally', {
  rast <- make_square_raster(1:100)

  npolys <- 13

  polys <- st_sf(fid = seq_len(npolys),
                 geometry = st_sfc(replicate(npolys, {
                   x <- runif(1, min=0, max=10)
                   y <- runif(1, min=0, max=10)
                   r <- runif(1, min=0, max=2)
                   make_circle(x, y, r, crs=sf::st_crs(rast))
                 }), crs=sf::st_crs(rast)))

  for (fun in list('sum', weighted.mean)) {
    for (input in list(polys, sf::st_geometry(polys))) {
      output <- capture.output(q <- exact_extract(rast, input, fun))
      lines <- strsplit(output, '\r', fixed=TRUE)[[1]]
      numlines <- lines[endsWith(lines, '%')]
      len <- nchar(numlines[1])
      pcts <- as.integer(substr(numlines, len - 3, len - 1))

      expect_length(pcts, 1 + npolys)
      expect_equal(pcts[1], 0)
      expect_equal(pcts[length(pcts)], 100)
      expect_false(is.unsorted(pcts))
    }
  }
})
