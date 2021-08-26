# Copyright (c) 2020-2021 ISciences, LLC.
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

context('exact_resample (terra)')

test_that("exact_resample supports SpatRaster arguments", {
  set.seed(123)

  # generate a random raster with a strange extent and resolution
  src <- raster::raster(matrix(runif(10000), nrow=100),
                               xmn=runif(1), xmx=runif(1) + 9,
                               ymn=runif(1), ymx=runif(1) + 9)

  # resample it to a raster with a larger grid with different resolution
  dst <- raster::raster(xmn=0, xmx=10, ymn=0, ymx=10, res=c(1, 2),
                        crs=raster::crs(src))

  dst <- exact_resample(src, dst, 'sum')

  dst_terra <- exact_resample(terra::rast(src), terra::rast(dst), 'sum')

  expect_equal(terra::rast(dst), dst_terra)
})
