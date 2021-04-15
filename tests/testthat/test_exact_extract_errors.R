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
context('exact_extract input validation')

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

test_that('Generic sfc_GEOMETRY fails if a feature is not polygonal', {
  rast <- make_square_raster(1:100)
  features <- st_as_sfc(c('POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))',
                          'POINT (2 7)'), crs=sf::st_crs(rast))

  expect_error(exact_extract(rast, features, 'sum', progress = FALSE),
               'must be polygonal')
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

  # unweighted, standard form
  for (fun in c(length, sum, median, mean, sd)) {
    expect_error(
      exact_extract(rast, poly, fun),
      'function .* not .* of the form')
  }

  expect_silent(exact_extract(rast, poly, weighted.mean))

  # unweighted, summarize_df
  expect_error(
    exact_extract(rast, poly, function() {}, summarize_df = TRUE),
    'function .* not .* of the form'
  )

  # weighted, standard form
  expect_error(
    exact_extract(rast, poly, weights = rast, fun = function(x, frac) {}),
    'function .* not .* of the form')
  expect_error(
    exact_extract(rast, poly, weights = rast, fun = function(x) {}),
    'function .* not .* of the form')
  expect_error(
    exact_extract(rast, poly, weights = rast, fun = function() {}),
    'function .* not .* of the form')

  # weighted, summarize_df
  expect_error(
    exact_extract(rast, poly, weights = rast, fun = function() {}, summarize_df = TRUE),
    'function .* not .* of the form'
  )
})

test_that('Error is raised for unknown summary operation', {
  rast <- make_square_raster(1:100)

  poly <- make_circle(5, 5, 3, sf::st_crs(rast))

  expect_error(exact_extract(rast, poly, 'whatimean'),
               'Unknown stat')
})

test_that('Error is raised if arguments passed without R summary function', {
  rast <- make_square_raster(1:100)

  poly <- make_circle(5, 5, 3, sf::st_crs(rast))

  expect_error(exact_extract(rast, poly, 'sum', na.rm=TRUE),
               'does not accept additional arguments')

  expect_error(
    exact_extract(rast, poly, cookie = FALSE),
    'Unexpected arguments'
  )
})

test_that('Error is raised for invalid max_cells_in_memory', {
  rast <- make_square_raster(1:100)
  poly <- make_circle(5, 5, 3, sf::st_crs(rast))

  expect_error(exact_extract(rast, poly, 'mean', max_cells_in_memory=-123),
               'Invalid.*max_cells')

  expect_error(
    exact_extract(rast, poly, 'mean', max_cells_in_memory = NA),
    'must be a single numeric')

  expect_error(
    exact_extract(rast, poly, 'mean', max_cells_in_memory = numeric()),
    'must be a single numeric')

  expect_error(
    exact_extract(rast, poly, 'mean', max_cells_in_memory = integer()),
    'must be a single numeric')

  expect_error(
    exact_extract(rast, poly, 'mean', max_cells_in_memory = NULL),
    'must be a single numeric')
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

  expect_error(exact_extract(rast, circles, 'sum', include_area = TRUE),
               'include_area must be FALSE')

  expect_error(exact_extract(rast, circles, 'sum', include_cell = TRUE),
               'include_cell must be FALSE')

  expect_error(exact_extract(rast, circles, 'sum', include_cols = 'fid'),
               'include_cols not supported')
})

test_that('Error is thrown when using include_cols or append_cols with nonexisting columns', {
  rast <- make_square_raster(1:100)

  circles <- st_sf(
    fid = c(2, 9),
    size = c('large', 'small'),
    geometry =  c(
    make_circle(5, 4, 2, sf::st_crs(rast)),
    make_circle(3, 1, 1, sf::st_crs(rast))))

  # append_cols specified but sfc has no attribute columns
  expect_error(
    exact_extract(rast, st_geometry(circles), 'mean', append_cols = 'fid', progress = FALSE),
    'only supported for sf')

  expect_error(
    exact_extract(rast, st_geometry(circles), weighted.mean, append_cols = 'fid', progress = FALSE),
    'only supported for sf')

  # append_cols specified for misspelled column
  expect_error(
    exact_extract(rast, circles, 'mean', append_cols = 'fidd', progress = FALSE),
    'undefined columns'
  )

  expect_error(
    exact_extract(rast, circles, weighted.mean, append_cols = 'fidd', progress = FALSE),
    'undefined columns'
  )

  # include_cols specified for sfc
  expect_error(
    exact_extract(rast, st_geometry(circles), include_cols = 'fidd', progress = FALSE),
    'only supported for sf'
  )

  # include_cols specified for misspelled column
  expect_error(
    exact_extract(rast, circles, include_cols = 'fidd', progress = FALSE),
    'undefined columns'
  )
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

test_that('Warning emitted when value raster is disaggregated', {
  r1 <- make_square_raster(1:100)
  r2 <- make_square_raster(runif(100))
  r1d <- raster::disaggregate(r1, 2)
  r2d <- raster::disaggregate(r2, 2)

  circle <- make_circle(2, 7, 3, sf::st_crs(r1))

  # no warning, values and weights have same resolution
  expect_silent(exact_extract(r1, circle, weights=r2))

  # no warning, values have higher resolution than weights
  expect_silent(exact_extract(r1d, circle, weights=r2))

  # warning, weights have higher resolution than values
  expect_warning(exact_extract(r1, circle, weights=r2d),
                 'value .* disaggregated')
})

test_that('Error raised when value raster is disaggregated and unweighted sum/count requested', {
  r1 <- make_square_raster(1:100)
  r1d <- raster::disaggregate(r1, 2)

  circle <- make_circle(2, 7, 3, sf::st_crs(r1))

  # no error, requested operations either expect disaggregation
  # or are not impacted by it
  expect_silent(exact_extract(r1, circle, c('weighted_sum', 'weighted_mean', 'mean'), weights=r1d))

  # on the other hand, "count" would be messed up by the disaggregation
  expect_error(exact_extract(r1, circle, c('weighted_sum', 'count'), weights=r1d),
               'raster is disaggregated')

  # as would "sum"
  expect_error(exact_extract(r1, circle, c('weighted_sum', 'count'), weights=r1d),
               'raster is disaggregated')

  # no problem if the weights are disaggregated, though
  expect_silent(exact_extract(r1d, circle, c('weighted_sum', 'count'), weights=r1))
})

test_that('We get an error if using stack_apply with incompatible stacks', {
  vals <- stack(replicate(3, make_square_raster(runif(100))))
  names(vals) <- c('a', 'b', 'c')

  weights <- stack(replicate(2, make_square_raster(runif(100))))
  names(weights) <- c('d', 'e')

  circle <- make_circle(2, 7, 3, sf::st_crs(vals))

  expect_error(
    exact_extract(vals, circle, function(v, c, w) 1, weights=weights, stack_apply=TRUE),
    "Can't apply")
})

test_that('Error thrown if summarize_df set where not applicable', {
  rast <- make_square_raster(1:100)
  circle <- make_circle(7.5, 5.5, 4, sf::st_crs(rast))

  expect_error(
    exact_extract(rast, circle, 'mean',  summarize_df = TRUE),
    'can only be used when .* function')

  expect_error(
    exact_extract(rast, circle, summarize_df = TRUE),
    'can only be used when .* function')
})

test_that('Error thrown if stack_apply set where not applicable', {
  rast <- make_square_raster(1:100)
  circle <- make_circle(7.5, 5.5, 4, sf::st_crs(rast))

  expect_error(
    exact_extract(rast, circle, stack_apply = TRUE),
    'can only be used when .* is a summary operation or function'
  )
})

test_that('Error thrown if append_cols set where not applicable', {
  rast <- make_square_raster(1:100)
  circle <- st_sf(make_circle(7.5, 5.5, 4, sf::st_crs(rast)))

  expect_error(
    exact_extract(rast, circle, append_cols = TRUE),
    'can only be used when .* is a summary operation or function'
  )
})

test_that('Error thrown if scalar args have length != 1', {
  rast <- make_square_raster(1:100)
  circle <- make_circle(7.5, 5.5, 4, sf::st_crs(rast))

  flags <- c(
    'coverage_area',
    'force_df',
    'full_colnames',
    'include_area',
    'include_cell',
    'include_xy',
    'progress',
    'stack_apply',
    'summarize_df')

  for (flag in flags) {
    base_args <- list(rast, circle)

    for (bad_value in list(logical(), c(TRUE, TRUE), NA)) {
      args <- base_args
      args[[flag]] <- bad_value

      expect_error(
        do.call(exact_extract, args),
        'must be TRUE or FALSE'
      )
    }
  }
})

test_that('Error thrown if fun is empty', {
  rast <- make_square_raster(1:100)
  circle <- make_circle(7.5, 5.5, 4, sf::st_crs(rast))

  expect_error(
    exact_extract(rast, circle, character()),
    'No summary operations'
  )
})

test_that('Error thrown if fun is incorrect type', {
  rast <- make_square_raster(1:100)
  circle <- make_circle(7.5, 5.5, 4, sf::st_crs(rast))

  expect_error(
    exact_extract(rast, circle, 44),
    'must be a character vector, function')
  expect_error(
    exact_extract(rast, circle, list(function() {}, function() {})),
    'must be a character vector, function')
})

test_that('Error thrown if default values have incorrect type/length', {
  rast <- make_square_raster(1:100)
  circle <- make_circle(7.5, 5.5, 4, sf::st_crs(rast))

  expect_error(
    exact_extract(rast, circle, 'mean', default_value = numeric()),
    'must be a single numeric value'
  )
  expect_error(
    exact_extract(rast, circle, 'mean', default_value = c(3, 8)),
    'must be a single numeric value'
  )
  expect_error(
    exact_extract(rast, circle, 'mean', default_value = NULL),
    'must be a single numeric value'
  )
  expect_error(
    exact_extract(rast, circle, 'mean', default_value = FALSE),
    'must be a single numeric value'
  )

})

