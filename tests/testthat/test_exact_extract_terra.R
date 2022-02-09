# Copyright (c) 2021-2022 ISciences, LLC.
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
context('exact_extract (terra)')

test_that('terra inputs supported (single layer)', {
  ras <- make_square_raster(1:100)

  terra_ras <- terra::rast(ras)

  circ <- make_circle(3, 2, 4, sf::st_crs(ras))

  expect_equal(
    exact_extract(ras, circ),
    exact_extract(terra_ras, circ)
  )

  expect_equal(
    exact_extract(ras, circ, 'mean'),
    exact_extract(terra_ras, circ, 'mean')
  )

  expect_equal(
    exact_extract(ras, circ, weighted.mean),
    exact_extract(terra_ras, circ, weighted.mean)
  )
})

test_that('terra inputs supported (single layer, weighted)', {
  ras <- make_square_raster(1:100)
  ras_w <- sqrt(ras)

  terra_ras <- terra::rast(ras)
  terra_ras_w <- terra::rast(ras_w)

  circ <- make_circle(3, 2, 4, sf::st_crs(ras))

  expect_equal(
    exact_extract(ras, circ, weights = ras_w),
    exact_extract(terra_ras, circ, weights = terra_ras_w)
  )

  expect_equal(
    exact_extract(ras, circ, 'weighted_mean', weights = ras_w),
    exact_extract(terra_ras, circ, 'weighted_mean', weights = terra_ras_w)
  )

  # mixed inputs supported: terra values, raster weights
  expect_equal(
    exact_extract(ras, circ, 'weighted_mean', weights = ras_w),
    exact_extract(terra_ras, circ, 'weighted_mean', weights = ras_w)
  )

  # mixed inputs supported: raster values, terra weights
  expect_equal(
    exact_extract(ras, circ, 'weighted_mean', weights = ras_w),
    exact_extract(ras, circ, 'weighted_mean', weights = terra_ras_w)
  )

  expect_equal(
    exact_extract(ras, circ, weighted.mean),
    exact_extract(terra_ras, circ, weighted.mean)
  )
})


test_that('terra inputs supported (multi-layer)', {
  stk <- raster::stack(list(a = make_square_raster(1:100),
                            b = make_square_raster(101:200)))

  terra_stk <- terra::rast(stk)

  circ <- make_circle(3, 2, 4, sf::st_crs(stk))

  expect_equal(
    exact_extract(stk, circ, 'mean'),
    exact_extract(terra_stk, circ, 'mean')
  )

  expect_equal(
    exact_extract(stk, circ),
    exact_extract(terra_stk, circ)
  )
})

test_that('terra inputs supported (weighted, multi-layer)', {
  stk <- raster::stack(list(a = make_square_raster(1:100),
                            a = make_square_raster(101:200)))
  stk <- terra::rast(stk)
  names(stk) <- c('a', 'a')

  ras <- terra::rast(make_square_raster(runif(100)))
  ras <- terra::disagg(ras, 2)

  circ <- make_circle(3, 2, 4, sf::st_crs(ras))

  expect_error(
    exact_extract(stk, circ, 'mean'),
    'names.*must be unique'
  )
})

test_that('include_* arguments supported for terra inputs', {
  ras <- make_square_raster(1:100)

  terra_ras <- terra::rast(ras)

  circ <- make_circle(3, 2, 4, sf::st_crs(ras))

  expect_equal(
    exact_extract(terra_ras, circ, include_cell = TRUE, include_xy = TRUE),
    exact_extract(ras, circ, include_cell = TRUE, include_xy = TRUE)
  )
})
