# Copyright (c) 2022 ISciences, LLC.
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

context('rasterize_polygons')

test_that('value is assigned to polygon with greatest coverage area', {
  polys <-
    st_as_sf(
      data.frame(id = 1:3,
                 geom = c(
                  'POLYGON ((10 0, 10 5, 5 5, 10 0))',
                  'POLYGON ((0 0, 10 0, 5 5, 1 10, 0 10, 0 0))',
                  'POLYGON ((5 5, 10 5, 10 10, 1 10, 5 5))'
                 )),
      wkt = 'geom'
    )

  rt <- terra::rast(xmin = 0, xmax = 10, ymin = 0, ymax = 10, res = 2)

  r <- rasterize_polygons(polys, rt)

  # lower-right is a tie, so it goes to the first feature encountered
  expect_equal(
    extract(r, cbind(9, 1))[[1]],
    1)

  # center tell is touched by all three, goes to polygon that covers the greatest area
  expect_equal(
    extract(r, cbind(5, 5))[[1]],
    2)
})

test_that('min_coverage excludes cells with small coverage area', {
  rt <- terra::rast(xmin = 0, xmax = 10, ymin = 0, ymax = 10, res = 1)
  circ <- make_circle(5, 5, 3.5, crs=st_crs(rt))
  circ_pieces <- st_sf(st_intersection(circ, st_make_grid(circ, 1)))

  cfrac <- coverage_fraction(rt, circ, crop = FALSE)[[1]]

  # by default, all touched cells are included in output
  r <- rasterize_polygons(circ_pieces, rt)
  expect_equal(
    values(cfrac) > 0,
    !is.na(values(r))
  )

  # min_coverage excludes cells with small coverage area
  r <- rasterize_polygons(circ_pieces, rt, min_coverage = 0.5)
  expect_equal(
    values(cfrac) > 0.5,
    !is.na(values(r))
  )
})

test_that('input type is preserved', {
  rt <- terra::rast(xmin = 0, xmax = 10, ymin = 0, ymax = 10, res = 2)
  rr <- raster::raster(rt)

  circ <- st_sf(make_circle(5, 5, 3.5, crs=st_crs(rt)))

  r <- rasterize_polygons(circ, rt)
  expect_s4_class(r, 'SpatRaster')

  r <- rasterize_polygons(circ, rr)
  expect_s4_class(r, 'RasterLayer')
})

test_that('no error when polygon does not intersect raster', {
  rt <- terra::rast(xmin = 0, xmax = 10, ymin = 0, ymax = 10, res = 2, crs=NA)

  circ <- st_sf(make_circle(500, 500, 3.5, crs=st_crs(rt)))

  r <- rasterize_polygons(circ, rt)

  expect_true(all(is.na(values(r))))
})

test_that('no error when polygon partially intersects raster', {
  rt <- terra::rast(xmin = 0, xmax = 10, ymin = 0, ymax = 10, res = 2, crs=NA)

  circ <- st_sf(make_circle(10, 5, 3.5, crs=st_crs(rt)))

  expect_invisible(
    r <- rasterize_polygons(circ, rt)
  )
})
