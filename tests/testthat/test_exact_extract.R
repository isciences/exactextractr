# Copyright (c) 2018-2021 ISciences, LLC.
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

test_that('sp inputs supported', {
  rast <- make_square_raster(1:100)

  circles <- c(
    make_circle(3, 2, 4, sf::st_crs(rast)),
    make_circle(7, 7, 2, sf::st_crs(rast))
  )
  circles_sf <- sf::st_sf(id = 1:2, geometry = circles)

  result <- exact_extract(rast, circles, 'mean', progress = FALSE)

  # SpatialPolygons
  circles_sp <- sf::as_Spatial(circles)
  result_sp <- exact_extract(rast, circles_sp, 'mean', progress = FALSE)
  expect_equal(result, result_sp)

  # SpatialPolygonsDataFrame
  circles_spdf <- sf::as_Spatial(circles_sf)
  result_spdf <- exact_extract(rast, circles_spdf, 'mean', progress = FALSE)
  expect_equal(result, result_spdf)
})

test_that('Generic sfc_GEOMETRY works if the features are polygonal', {
  rast <- make_square_raster(1:100)
  polys <- st_as_sfc(c('POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))',
                       'MULTIPOLYGON (((2 2, 4 2, 4 4, 2 4, 2 2)), ((4 4, 8 4, 8 8, 4 8, 4 4)))'),
                     crs=sf::st_crs(rast))

  expect_equal(exact_extract(rast, polys, 'count', progress = FALSE),
               c(4, 4+16))
})

test_that('GeometryCollections are supported if they are polygonal', {
  rast <- make_square_raster(1:100)
  gc <- st_as_sfc('GEOMETRYCOLLECTION( POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0)), POLYGON ((2 2, 4 2, 4 4, 2 4, 2 2)))',
                  crs = st_crs(rast))

  mp <- st_as_sfc('MULTIPOLYGON (((0 0, 2 0, 2 2, 0 2, 0 0)), ((2 2, 4 2, 4 4, 2 4, 2 2)))',
                  crs = st_crs(rast))

  expect_equal(
    exact_extract(rast, gc),
    exact_extract(rast, mp)
  )
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

  expect_named(means, c('weighted.mean.a', 'weighted.mean.b'))

  # results are same as we would get by processing layers independently
  for (i in 1:raster::nlayers(stk)) {
    expect_equal(means[, i], exact_extract(stk[[i]], circles, weighted.mean, progress=FALSE))
  }
})

test_that('Layers of a RasterBrick can be processed independently with stack_apply', {
  # https://github.com/isciences/exactextractr/issues/54

  data <- matrix(1:100, nrow=10, byrow=TRUE)
  data[7:10, 1:4] <- NA # cut out lower-left corner
  rast <- raster::raster(
    data,
    xmn=0, xmx=10, ymn=0, ymx=10,
    crs='+proj=longlat +datum=WGS84'
  )
  rast_brick <- brick(rast, rast)
  square <- make_rect(3.5, 3.5, 4.5, 4.5, sf::st_crs(rast))

  expect_equal(
    exact_extract(rast_brick, square, weighted.mean, stack_apply = T),
    data.frame(weighted.mean.layer.1 = NA_real_,
               weighted.mean.layer.2 = NA_real_))
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

  # error when trying to use a non-raster as weights
  expect_error(exact_extract(stk, circle, 'weighted_mean', weights='stk'),
               "Weights must be a Raster")
})

test_that('We get acceptable default values when processing a polygon that does
           not intersect the raster', {
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

test_that('Coverage area can be output instead of coverage fraction (projected)', {
  rast_utm <- disaggregate(make_square_raster(1:100), c(2, 3))

  circle <- make_circle(5, 5, 5, crs=st_crs(rast_utm))

  df_frac <- exact_extract(rast_utm, circle, include_area = TRUE)[[1]]

  df_area <- exact_extract(rast_utm, circle, coverage_area = TRUE)[[1]]
  expect_named(df_area, c('value', 'coverage_area'))

  expect_equal(df_frac$coverage_fraction * df_frac$area,
               df_area$coverage_area)
})

test_that('Coverage area can be output instead of coverage fraction (geographic)', {
  rast <- raster::raster(matrix(1:54000, ncol=360),
                         xmn=-180, xmx=180, ymn=-65, ymx=85,
                         crs='+proj=longlat +datum=WGS84')

  suppressMessages({
    circle <- make_circle(0, 45, 15, crs=st_crs(rast))
  })

  df_frac <- exact_extract(rast, circle, include_area = TRUE)[[1]]

  df_area <- exact_extract(rast, circle, coverage_area = TRUE)[[1]]

  expect_equal(df_frac$coverage_fraction * df_frac$area,
               df_area$coverage_area)
})

test_that('coverage_area argument can be used with named summary operations', {
  rast1 <- raster(matrix(1:54000, ncol=360),
                  xmn=-180, xmx=180, ymn=-65, ymx=85,
                  crs='+proj=longlat +datum=WGS84')
  rast2 <- sqrt(rast1)

  suppressMessages({
    circle <- make_circle(0, 45, 15, crs=st_crs(rast1))
  })

  # using only area as weighting
  expect_equal(exact_extract(rast1, circle, 'weighted_mean', weights = 'area'),
               exact_extract(rast1, circle, 'mean', coverage_area = TRUE))

  # using area x weight as weighting
  expect_equal(
    exact_extract(rast1, circle, 'weighted_mean', weights = rast2, coverage_area = TRUE),

    exact_extract(rast1, circle, fun = function(x, cov, w) {
      weighted.mean(x, cov * w$rast2 * w$area)
    }, weights = stack(list(rast2 = rast2, area = area(rast2)))),

    tol = 1e-2
  )
})

test_that('We can weight with cell areas (projected coordinates)', {
  rast_utm <- raster(matrix(1:100, ncol=10),
                     xmn=0, xmx=5,
                     ymn=0, ymx=5,
                     crs='+init=epsg:26918')
  circle1 <- make_circle(5, 5, 5, crs=st_crs(rast_utm))

  # for projected (Cartesian coordinates), means with cell area and
  # coverage fraction are the same
  expect_equal(exact_extract(rast_utm, circle1, 'mean'),
               exact_extract(rast_utm, circle1, 'weighted_mean', weights='area'))

  # same result with R summary function
  expect_equal(
    exact_extract(rast_utm, circle1, 'weighted_mean', weights='area'),
    exact_extract(rast_utm, circle1, function(x,c,w) {
      weighted.mean(x, c*w)
    }, weights='area'),
    1e-5
  )

  # name doesn't pop out in data frame columns
  expect_named(
    exact_extract(rast_utm, circle1, c('sum', 'weighted_mean'), weights='area', force_df = TRUE),
    c('sum', 'weighted_mean'))

  # sums differ by the cell area
  expect_equal(prod(res(rast_utm)) * exact_extract(rast_utm, circle1, 'sum'),
               exact_extract(rast_utm, circle1, 'weighted_sum', weights='area'))

  # when using area weighting, disaggregating does not affect the sum
  expect_equal(exact_extract(rast_utm, circle1, 'weighted_sum', weights='area'),
               exact_extract(disaggregate(rast_utm, 8), circle1, 'weighted_sum', weights='area'))
})

test_that('We can weight with cell areas (geographic coordinates)', {
  rast <- raster::raster(matrix(1:54000, ncol=360),
                         xmn=-180, xmx=180, ymn=-65, ymx=85,
                         crs='+proj=longlat +datum=WGS84')

  accuracy_pct_tol <- 0.01

  suppressMessages({
    circle <- make_circle(0, 45, 15, crs=st_crs(rast))
  })

  # result is reasonably close to what we get with raster::area, which uses
  # a geodesic calculation
  expected <- exact_extract(rast, circle, 'weighted_sum', weights = area(rast) * 1e6)
  actual <- exact_extract(rast, circle, 'weighted_sum', weights = 'area')

  expect_true(abs(actual - expected) / expected < accuracy_pct_tol)
})

test_that('Correct results obtained when max_cells_in_memory is limited', {
  rast <- make_square_raster(1:100)
  poly <- make_circle(5, 5, 3, sf::st_crs(rast))

  expect_equal(exact_extract(rast, poly, 'mean'),
               exact_extract(rast, poly, 'mean', max_cells_in_memory=1))
})

test_that('Weighted stats work when polygon is contained in weight raster but
          only partially contained in value raster', {
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

test_that('when force_df = TRUE, exact_extract always returns a data frame', {
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
               exact_extract(rast, circles[2], function(x, c) {
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
                         exact_extract(rast, circles, weighted.mean, force_df = TRUE, progress = FALSE)))
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

test_that('Both value and weighting rasters can be a stack', {
  vals <- stack(replicate(3, make_square_raster(runif(100))))
  names(vals) <- c('a', 'b', 'c')

  weights <- stack(replicate(2, make_square_raster(rbinom(100, 2, 0.5))))
  names(weights) <- c('w1', 'w2')

  circle <- make_circle(2, 7, 3, sf::st_crs(vals))

  extracted <- exact_extract(vals, circle, weights=weights)[[1]]

  expect_named(extracted, c('a', 'b', 'c', 'w1', 'w2', 'coverage_fraction'))

  # stack of values, stack of weights: both passed as data frames
  exact_extract(vals, circle, function(v, c, w) {
    expect_true(is.data.frame(v))
    expect_true(is.data.frame(w))
    expect_named(v, names(vals))
    expect_named(w, names(weights))
  }, weights = weights)

  # stack of values, single layer of weights: weights passed as vector
  exact_extract(vals, circle, function(v, c, w) {
    expect_true(is.data.frame(v))
    expect_true(is.vector(w))
  }, weights = weights[[1]])

  # single layer of values, stack of weights: values passed as vector
  exact_extract(vals[[1]], circle, function(v, c, w) {
    expect_true(is.vector(v))
    expect_true(is.data.frame(w))
  }, weights = weights)

  # single layer of values, single layer of weights: both passed as vector
  exact_extract(vals[[1]], circle, function(v, c, w) {
    expect_true(is.vector(v))
    expect_true(is.vector(w))
  }, weights = weights[[1]])
})

test_that('Named summary operations support both stacks of values and weights', {
  vals <- stack(replicate(3, make_square_raster(runif(100))))
  names(vals) <- c('v1', 'v2', 'v3')

  weights <- stack(replicate(3, make_square_raster(rbinom(100, 2, 0.5))))
  names(weights) <- c('w1', 'w2', 'w3')

  circle <- make_circle(2, 7, 3, sf::st_crs(vals))

  stats <- c('sum', 'weighted_mean')

  # stack of values, stack of weights: values and weights are applied pairwise
  result <- exact_extract(vals, circle, stats, weights=weights)
  expect_named(result, c(
    'sum.v1', 'sum.v2', 'sum.v3',
    'weighted_mean.v1.w1', 'weighted_mean.v2.w2', 'weighted_mean.v3.w3'))
  expect_equal(result$sum.v1, exact_extract(vals[[1]], circle, 'sum'))
  expect_equal(result$sum.v2, exact_extract(vals[[2]], circle, 'sum'))
  expect_equal(result$sum.v3, exact_extract(vals[[3]], circle, 'sum'))
  expect_equal(result$weighted_mean.v1, exact_extract(vals[[1]], circle, 'weighted_mean', weights=weights[[1]]))
  expect_equal(result$weighted_mean.v2, exact_extract(vals[[2]], circle, 'weighted_mean', weights=weights[[2]]))
  expect_equal(result$weighted_mean.v3, exact_extract(vals[[3]], circle, 'weighted_mean', weights=weights[[3]]))

  # stack of values, layer of weights: weights are recycled
  result <- exact_extract(vals, circle, stats, weights=weights[[1]])
  expect_named(result, c(
    'sum.v1', 'sum.v2', 'sum.v3',
    'weighted_mean.v1', 'weighted_mean.v2', 'weighted_mean.v3'))
  expect_equal(result$sum.v1, exact_extract(vals[[1]], circle, 'sum'))
  expect_equal(result$sum.v2, exact_extract(vals[[2]], circle, 'sum'))
  expect_equal(result$sum.v3, exact_extract(vals[[3]], circle, 'sum'))
  expect_equal(result$weighted_mean.v1, exact_extract(vals[[1]], circle, 'weighted_mean', weights=weights[[1]]))
  expect_equal(result$weighted_mean.v2, exact_extract(vals[[2]], circle, 'weighted_mean', weights=weights[[1]]))
  expect_equal(result$weighted_mean.v3, exact_extract(vals[[3]], circle, 'weighted_mean', weights=weights[[1]]))

  # layer of values, stack of weights: values are recycled
  result <- exact_extract(vals[[3]], circle, stats, weights=weights)
  expect_named(result, c('sum', 'weighted_mean.w1', 'weighted_mean.w2', 'weighted_mean.w3'))
  expect_equal(result$sum, exact_extract(vals[[3]], circle, 'sum'))
  expect_equal(result$weighted_mean.w1, exact_extract(vals[[3]], circle, 'weighted_mean', weights=weights[[1]]))
  expect_equal(result$weighted_mean.w2, exact_extract(vals[[3]], circle, 'weighted_mean', weights=weights[[2]]))
  expect_equal(result$weighted_mean.w3, exact_extract(vals[[3]], circle, 'weighted_mean', weights=weights[[3]]))
})

test_that('We can use stack_apply with both values and weights', {
  vals <- stack(replicate(3, make_square_raster(runif(100))))
  names(vals) <- c('v1', 'v2', 'v3')

  weights <- stack(replicate(3, make_square_raster(rbinom(100, 2, 0.5))))
  names(weights) <- c('w1', 'w2', 'w3')

  circle <- make_circle(2, 7, 3, sf::st_crs(vals))

  weighted_mean <- function(v, c, w) {
    expect_equal(length(v), length(c))
    expect_equal(length(v), length(w))

    weighted.mean(v, c*w)
  }

  # stack of values, stack of weights: values and weights are applied pairwise
  result <- exact_extract(vals, circle, weighted_mean, weights = weights, stack_apply = TRUE)
  expect_named(result, c('fun.v1.w1', 'fun.v2.w2', 'fun.v3.w3'))
  expect_equal(result$fun.v2.w2,
               exact_extract(vals[[2]], circle, 'weighted_mean', weights=weights[[2]]),
               tol = 1e-6)

  # stack of values, layer of weights: weights are recycled
  result <- exact_extract(vals, circle, weighted_mean, weights = weights[[2]], stack_apply = TRUE, full_colnames = TRUE)
  expect_named(result, c('fun.v1.w2', 'fun.v2.w2', 'fun.v3.w2'))
  expect_equal(result$fun.v1.w2,
               exact_extract(vals[[1]], circle, 'weighted_mean', weights=weights[[2]]),
               tol = 1e-6)

  # layer of values, stack of weights: values are recycled
  result <- exact_extract(vals[[3]], circle, weighted_mean, weights = weights, stack_apply = TRUE, full_colnames = TRUE)
  expect_named(result, c('fun.v3.w1', 'fun.v3.w2', 'fun.v3.w3'))
  expect_equal(result$fun.v3.w1,
               exact_extract(vals[[3]], circle, 'weighted_mean', weights=weights[[1]]),
               tol = 1e-6)
})

test_that('Layers are implicity renamed if value layers have same name as weight layers', {
  # this happens when a stack is created and no names are provided
  # raster package assigns layer.1, layer.2
  # here we assign our own identical names to avoid relying on raster package
  # implementation detail
  vals <- stack(replicate(2, make_square_raster(runif(100))))
  names(vals) <- c('a', 'b')

  weights <- stack(replicate(2, make_square_raster(runif(100))))
  names(weights) <- c('a', 'b')

  circle <- make_circle(2, 7, 3, sf::st_crs(vals))

  result <- exact_extract(vals, circle, weights=weights)[[1]]
  expect_named(result, c('a', 'b', 'a.1', 'b.1', 'coverage_fraction'))
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

test_that('generated column names follow expected pattern', {
  values <- c('v1', 'v2', 'v3')
  weights <- c('w1', 'w2', 'w3')

  stats <- c('mean', 'weighted_mean')

  test_mean <- function(x, c) { weighted.mean(x, c) }

  # layer of values, no weights
  # named summary operations
  expect_equal(.resultColNames(values[[2]], NULL, c('mean', 'sum'), TRUE),
               c('mean.v2', 'sum.v2'))
  expect_equal(.resultColNames(values[[2]], NULL, c('mean', 'sum'), FALSE),
               c('mean', 'sum'))

  # generic method (we can recover its name)
  expect_equal(.resultColNames(values[[2]], NULL, weighted.mean, TRUE),
               'weighted.mean.v2')
  expect_equal(.resultColNames(values[[2]], NULL, weighted.mean, FALSE),
               'weighted.mean')

  # regular function (we can't recover its name)
  expect_equal(.resultColNames(values[[2]], NULL, test_mean, TRUE),
               'fun.v2')
  expect_equal(.resultColNames(values[[2]], NULL, test_mean, FALSE),
               'fun')

  # stack of values, no weights
  for (full_colnames in c(TRUE, FALSE)) {
    expect_equal(.resultColNames(values, NULL, c('mean', 'sum'), full_colnames),
                 c('mean.v1', 'mean.v2', 'mean.v3',
                   'sum.v1', 'sum.v2', 'sum.v3'))
    expect_equal(.resultColNames(values, NULL, test_mean, full_colnames),
                 c('fun.v1', 'fun.v2', 'fun.v3'))
  }

  # values, weights processed in parallel
  for (full_colnames in c(TRUE, FALSE)) {
    expect_equal(.resultColNames(values, weights, stats, full_colnames),
                 c('mean.v1', 'mean.v2', 'mean.v3',
                   'weighted_mean.v1.w1', 'weighted_mean.v2.w2', 'weighted_mean.v3.w3'))
    expect_equal(.resultColNames(values, weights, test_mean, full_colnames),
                 c('fun.v1.w1', 'fun.v2.w2', 'fun.v3.w3'))
  }

  # values recycled (full names)
  expect_equal(.resultColNames(values[1], weights, stats, TRUE),
               c('mean.v1', 'mean.v1', 'mean.v1',
                 'weighted_mean.v1.w1', 'weighted_mean.v1.w2', 'weighted_mean.v1.w3'))
  expect_equal(.resultColNames(values[1], weights, test_mean, TRUE),
               c('fun.v1.w1', 'fun.v1.w2', 'fun.v1.w3'))

  # here the values are always the same so we don't bother adding them to the names
  expect_equal(.resultColNames(values[1], weights, stats, FALSE),
               c('mean', 'mean', 'mean',
                 'weighted_mean.w1', 'weighted_mean.w2', 'weighted_mean.w3'))
  expect_equal(.resultColNames(values[1], weights, test_mean, FALSE),
               c('fun.w1', 'fun.w2', 'fun.w3'))

  # weights recycled (full names)
  expect_equal(.resultColNames(values, weights[1], stats, TRUE),
               c('mean.v1', 'mean.v2', 'mean.v3',
                 'weighted_mean.v1.w1', 'weighted_mean.v2.w1', 'weighted_mean.v3.w1'))
  expect_equal(.resultColNames(values, weights[1], test_mean, TRUE),
               c('fun.v1.w1', 'fun.v2.w1', 'fun.v3.w1'))

  # here the weights are always the same so we don't bother adding them to the name
  expect_equal(.resultColNames(values, weights[1], stats, FALSE),
               c('mean.v1', 'mean.v2', 'mean.v3',
                 'weighted_mean.v1', 'weighted_mean.v2', 'weighted_mean.v3'))
  expect_equal(.resultColNames(values, weights[1], test_mean, FALSE),
               c('fun.v1', 'fun.v2', 'fun.v3'))
})

test_that('We can replace NA values in the value and weighting rasters with constants', {
  set.seed(05401)

  x <- runif(100)
  x[sample(length(x), 0.5*length(x))] <- NA

  y <- runif(100)
  y[sample(length(y), 0.5*length(y))] <- NA

  rx <- make_square_raster(x)
  ry <- make_square_raster(y)

  poly <- make_circle(4.5, 4.8, 4, crs=st_crs(rx))

  # manually fill the missing values with 0.5 and missing weights with 0.3
  rx_filled <- make_square_raster(ifelse(is.na(x), 0.5, x))
  ry_filled <- make_square_raster(ifelse(is.na(y), 0.3, y))

  expected <- exact_extract(rx_filled, poly, 'weighted_mean', weights = ry_filled)

  # fill values on the fly and verify that we get the same result
  expect_equal(
    exact_extract(rx, poly, 'weighted_mean', weights = ry, default_value = 0.5, default_weight = 0.3),
    expected)

  # check same calculation but using R summary function
  expect_equal(
    exact_extract(rx, poly, weights = ry, default_value = 0.5, default_weight = 0.3,
      fun = function(value, cov_frac, weight) {
        weighted.mean(value, cov_frac*weight)
      }),
    expected, 1e-6)

  # check substitution in raw returned values
  expect_equal(
    which(is.na(exact_extract(rx, poly)[[1]]$value)),
    which(44 == exact_extract(rx, poly, default_value = 44)[[1]]$value)
  )
})

test_that('All summary function arguments combined when summarize_df = TRUE', {
  rast <- make_square_raster(1:100)

  values <- stack(list(a = rast - 1,
                       b = rast,
                       c = rast + 1))

  weights <- sqrt(values)
  names(weights) <- c('d', 'e', 'f')

  circle <- st_sf(
    id = 77,
    make_circle(7.5, 5.5, 4, sf::st_crs(rast)))

  # in the tests below, we check names inside the R summary function
  # to verify that our checks were actually hit, we have the summary
  # function return NULL and check for it with `expect_null`.

  # values only
  expect_null(
    exact_extract(values, circle, summarize_df = TRUE, fun = function(df) {
      expect_named(df, c('a', 'b', 'c', 'coverage_fraction'))
      NULL
    })[[1]])

  expect_null(
    exact_extract(rast, circle, coverage_area = TRUE, summarize_df = TRUE, fun = function(df) {
      expect_named(df, c('value', 'coverage_area'))
      NULL
    })[[1]])

  expect_null(
    exact_extract(values[[1]], circle, coverage_area = TRUE, summarize_df = TRUE, fun = function(df) {
      expect_named(df, c('value', 'coverage_area'))
      NULL
    })[[1]])

  # values and weights
  expect_null(
    exact_extract(values, circle, summarize_df = TRUE, fun = function(df) {
      expect_named(df, c('a', 'b', 'c', 'd', 'e', 'f', 'coverage_fraction'))
      NULL
    }, weights = weights)[[1]])

  expect_null(
    exact_extract(values, circle,
                  include_cell = TRUE,
                  include_xy = TRUE,
                  include_area = TRUE,
                  include_cols = 'id',
                  summarize_df = TRUE,
                  fun = function(df, extra_arg) {
      expect_named(df,
        c('id', 'a', 'b', 'c', 'd', 'e', 'f', 'x', 'y', 'cell', 'area', 'coverage_fraction'))
      expect_equal(extra_arg, 600)
      NULL
    }, weights = weights, extra_arg = 600)[[1]])

  # values and weights, stack_apply = TRUE
  expect_equal(
    exact_extract(values, circle, weights = weights, summarize_df = TRUE, stack_apply = TRUE,
                  fun = function(df, extra_arg) {
                    expect_named(df, c('value', 'weight', 'coverage_fraction'))
                    extra_arg
                  }, extra_arg = 30809),
    data.frame(fun.a.d = 30809,
               fun.b.e = 30809,
               fun.c.f = 30809))
})

test_that('floating point errors do not cause an error that
          "logical subsetting requires vectors of identical size"', {
  rast <- raster(matrix(1:100, nrow=10), xm=0, xmx=1, ymn=0, ymx=1)
  poly <- make_rect(0.4, 0.7, 0.5, 0.8, crs = st_crs(rast))

  val <- exact_extract(rast, poly, weights = rast, fun = NULL, include_cell = TRUE)[[1]]

  expect_equal(val$value, rast[val$cell])
  expect_equal(val$weight, rast[val$cell])
})

test_that("append_cols works correctly when summary function returns multi-row data frame", {
  rast <- make_square_raster(1:100)

  circles <- st_sf(
    id = c('a', 'b'),
    geom = c(
      make_circle(3, 2, 4, sf::st_crs(rast)),
      make_circle(7, 7, 2, sf::st_crs(rast))
  ))

  expect_silent({
    result <- exact_extract(rast,
                            circles,
                            function(x, cov)
                              data.frame(x = 1:3,
                                         x2 = 4:6),
                            append_cols = 'id',
                            progress = FALSE)
  })

  expect_named(result, c('id', 'x', 'x2'))

  expect_equal(result$id, c('a', 'a', 'a', 'b', 'b', 'b'))
  expect_equal(result$x,  c(1:3, 1:3))
  expect_equal(result$x2,  c(4:6, 4:6))
})

test_that("append_cols works correctly when summary function returns vector with length > 1", {
  rast <- make_square_raster(1:100)

  circles <- st_sf(
    id = c('a', 'b'),
    geom = c(
      make_circle(3, 2, 4, sf::st_crs(rast)),
      make_circle(7, 7, 2, sf::st_crs(rast))
  ))

  expect_silent({
    result <- exact_extract(rast,
                            circles,
                            function(x, cov)
                              1:3,
                            append_cols = 'id',
                            progress = FALSE)
  })

  expect_named(result, c('id', 'result'))

  expect_equal(result$id, c('a', 'a', 'a', 'b', 'b', 'b'))
  expect_equal(result$result,  c(1:3, 1:3))
})

test_that("append_cols works correctly when summary function returns data frame with length 0", {
  rast <- make_square_raster(1:100)

  circles <- st_sf(
    id = c('a', 'b'),
    geom = c(
      make_circle(3, 2, 4, sf::st_crs(rast)),
      make_circle(7, 7, 2, sf::st_crs(rast))
  ))

  expect_silent({
    result <- exact_extract(rast,
                            circles,
                            function(x, cov)
                              data.frame(x = character(0),
                                         x2 = numeric(0)),
                            append_cols = 'id',
                            progress = FALSE)
  })

  expect_named(result, c('id', 'x', 'x2'))
  expect_equal(nrow(result), 0)
  expect_equal(class(result$id), class(circles$id))
  expect_equal(class(result$x), 'character')
  expect_equal(class(result$x2), 'numeric')
})

