# Copyright (c) 2021-2023 ISciences, LLC.
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

context('block size detection')

test_that('blockSize reports block size in row-col order', {
  landcov_fname <- system.file(file.path('sao_miguel', 'clc2018_v2020_20u1.tif'),
                               package='exactextractr')

  expect_equal(
    .blockSize(raster::raster(landcov_fname)),
    c(2, 3840)
  )

  expect_equal(
    .blockSize(terra::rast(landcov_fname)),
    c(2, 3840)
  )

  # netCDF uses a different code path, so copy our test input to netCDF format
  # and repeat. GDAL doesn't let us control the block size, so hopefully it
  # is stable.
  if ('netCDF' %in% terra::gdal(drivers=TRUE)$name) {
    nc_fname <- tempfile(fileext = '.nc')
    suppressWarnings({
      terra::writeRaster(terra::rast(landcov_fname), nc_fname, gdal=c('FORMAT=NC4', 'COMPRESS=DEFLATE'))
    })

    expect_equal(
      .blockSize(raster::raster(nc_fname)),
      c(1, 3840)
    )

    expect_equal(
      .blockSize(terra::rast(nc_fname)),
      c(1, 3840)
    )

    file.remove(nc_fname)
  }
})
