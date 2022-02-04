# Copyright (c) 2018-2022 ISciences, LLC.
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
context('exact_extract include* arguments')

test_that('when include_xy = TRUE, center coordinates area included in the output', {
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

test_that('when include_area = TRUE, cell areas are included in output (geographic) and are accurate to 1%', {
  rast <- raster::raster(matrix(1:54000, ncol=360),
                         xmn=-180, xmx=180, ymn=-65, ymx=85,
                         crs='+proj=longlat +datum=WGS84')
  accuracy_pct_tol <- 0.01

  suppressMessages({
    circle <- make_circle(0, 45, 15, crs=st_crs(rast))
  })

  results <- exact_extract(rast, circle, include_cell = TRUE, include_area = TRUE)[[1]]

  expected_areas <- raster::area(rast)[results$cell]
  actual_areas <- results$area / 1e6

  expect_true(all(abs(actual_areas - expected_areas) / expected_areas < accuracy_pct_tol))
})

test_that('when include_area = TRUE, cell areas are included in output (projected)', {
  rast_utm <- make_square_raster(1:100)

  circle <- make_circle(5, 5, 5, crs=st_crs(rast_utm))

  areas <- exact_extract(rast_utm, circle, include_area = TRUE)[[1]]$area
  expect_true(all(areas == 1))
})

test_that('include_cols copies columns from the source data frame to the returned data frames', {
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

test_that('When disaggregating values, xy coordinates refer to disaggregated grid', {
  rast <- make_square_raster(1:100)
  rast2 <- raster::disaggregate(rast, 4)

  circle <- make_circle(7.5, 5.5, 0.4, sf::st_crs(rast))

  xy_disaggregated <- exact_extract(rast2, circle, include_xy = TRUE)[[1]][, c('x', 'y')]

  suppressWarnings({
    xy_weighted <- exact_extract(rast, circle, include_xy = TRUE, weights = rast2)[[1]][, c('x', 'y')]
    xy_weighted2 <- exact_extract(rast2, circle, include_xy = TRUE, weights = rast)[[1]][, c('x', 'y')]
  })

  expect_equal(xy_weighted, xy_disaggregated)
  expect_equal(xy_weighted2, xy_disaggregated)
})

test_that('When value and weighting rasters have different grids, cell numbers refer to value raster', {
  anom <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, res=10)
  values(anom) <- rnorm(length(anom))

  pop <- raster(xmn=-180, xmx=180, ymn=-65, ymx=85, res=5)
  values(pop) <- rlnorm(length(pop))

  circle <- make_circle(17, 21, 18, sf::st_crs(anom))

  suppressWarnings({
    extracted <- exact_extract(anom, circle, weights=pop, include_cell=TRUE)[[1]]
  })

  expect_equal(extracted$value, anom[extracted$cell])
})

test_that('include_ arguments supported with weighted summary function', {
  rast1 <- 5 + make_square_raster(1:100)
  rast2 <- make_square_raster(runif(100))

  circle <- st_sf(
    id = 77,
    make_circle(7.5, 5.5, 4, sf::st_crs(rast1)))

  x <- exact_extract(rast1, circle, function(v, c, w) {
    expect_is(v, 'data.frame')
    expect_named(v, c('value', 'id'))
    expect_true(all(v$id ==  77))

    expect_is(c, 'numeric')
    expect_is(w, 'numeric')
  }, weights=rast2, include_cols = 'id')

  x <- exact_extract(rast1, circle, function(v, c, w) {
    expect_is(v, 'data.frame')
    expect_named(v, c('value', 'id', 'x', 'y', 'cell'))
    expect_true(all(v$id ==  77))
    expect_equal(v$value, rast1[v$cell])
    expect_equal(w, rast2[v$cell])
    expect_equal(v$x, raster::xFromCell(rast1, v$cell))
    expect_equal(v$y, raster::yFromCell(rast1, v$cell))

    expect_is(c, 'numeric')
    expect_is(w, 'numeric')
  }, weights=rast2, include_cols = 'id', include_cell = TRUE, include_xy = TRUE)
})

test_that('we get a zero-row data frame for a polygon not intersecting a raster', {
  # https://github.com/isciences/exactextractr/issues/68

  rast <- raster(matrix(0, nrow = 100, ncol = 100))

  nonoverlap_poly <- st_sf(st_sfc(st_polygon(list(matrix(c(0, 0, 1, 0, 1,
                                                           -0.25, 0, -0.25, 0, 0),
                                                         ncol = 2, byrow = TRUE)))))

  df <- exact_extract(rast, nonoverlap_poly)[[1]]
  expect_named(df, c('value', 'coverage_fraction'))
  expect_equal(nrow(df), 0)

  df <- exact_extract(rast, nonoverlap_poly, include_xy = TRUE)[[1]]
  expect_named(df, c('value', 'x', 'y', 'coverage_fraction'))
  expect_equal(nrow(df), 0)

  df <- exact_extract(rast, nonoverlap_poly, include_cell = TRUE)[[1]]
  expect_named(df, c('value', 'cell', 'coverage_fraction'))
  expect_equal(nrow(df), 0)

  df <- exact_extract(rast, nonoverlap_poly, include_area = TRUE)[[1]]
  expect_named(df, c('value', 'area', 'coverage_fraction'))
  expect_equal(nrow(df), 0)
})
