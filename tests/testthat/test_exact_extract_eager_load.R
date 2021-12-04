# Copyright (c) 2018-2021 ISciences, LLC.
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

context('exact_extract eager loading')

test_that("message emitted when working area doesn't fit in memory", {
  rast_fname <- system.file(file.path('sao_miguel', 'clc2018_v2020_20u1.tif'),
                            package = 'exactextractr')

  poly_fname <- system.file(file.path('sao_miguel', 'concelhos.gpkg'),
                            package = 'exactextractr')

  r <- terra::rast(rast_fname)
  polys <- st_read(poly_fname, quiet = TRUE)

  # no output when everything fits in memory
  capture.output({
    msg <- capture_messages({
      exact_extract(r, polys, 'mode', progress = TRUE, max_cells_in_memory = 1e7)
    })
  })
  expect_equal(msg, character())

  # message emitted when it doesn't fit
  capture.output({
    expect_message(
      exact_extract(r, polys, 'mode', progress = TRUE, max_cells_in_memory = 1e6),
      'Cannot preload'
    )
  })

  # if progress is disabled, so are hints
  expect_silent(
    exact_extract(r, polys, 'mode', progress = FALSE, max_cells_in_memory = 1e6)
  )

  # get additional warning by blowing out the GDAL block cache
  prevCacheSize <- terra::gdalCache()
  terra::gdalCache(1)

  capture.output({
    expect_message(
      exact_extract(r, polys, 'mode', max_cells_in_memory = 1e6),
      'GDAL block size cache is only 1 MB'
    )
  })

  # get additional warning if we are using a RasterStack
  capture.output({
    expect_message(
      exact_extract(stack(r), polys, 'mode', max_cells_in_memory = 1e6),
      'It is recommended to use a SpatRaster'
    )
  })

  terra::gdalCache(prevCacheSize)
})
