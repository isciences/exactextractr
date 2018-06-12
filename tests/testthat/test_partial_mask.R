# Copyright (c) 2018 ISciences, LLC.
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

context('partial_mask')

test_that("Partial mask function works", {
  # This test just verifies a successful journey from R
  # to C++ and back. The correctness of the algorithm
  # is tested at the C++ level.
  square <- sf::st_sfc(sf::st_polygon(
    list(
      matrix(
        c(0.5, 0.5, 2.5, 0.5, 2.5, 2.5, 0.5, 2.5, 0.5, 0.5),
        ncol=2,
        byrow=TRUE))))

  rast <- raster::raster(xmn=0, xmx=3, ymn=0, ymx=3, nrows=3, ncols=3)

  weights <- partial_mask(rast, square)

  expect_equal(weights[[1]],
               rbind(
                 c(0.25, 0.5, 0.25),
                 c(0.50, 1.0, 0.50),
                 c(0.25, 0.5, 0.25)
               ), check.attributes=FALSE)
})
