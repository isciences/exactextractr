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

default_proj <- '+init=epsg:26918' # UTM 18N; avoid wgs84 to keep cartesian calcs in sf

make_rect <- function(xmin, ymin, xmax, ymax, crs) {
  sf::st_sfc(
    sf::st_polygon(
      list(
        matrix(
          c(xmin, ymin,
            xmax, ymin,
            xmax, ymax,
            xmin, ymax,
            xmin, ymin),
          ncol=2,
          byrow=TRUE))),
    crs=crs)
}

make_circle <- function(x, y, r, crs) {
  suppressWarnings(sf::st_buffer(
    sf::st_sfc(
      sf::st_point(c(x, y)),
      crs=crs),
    r))
}

make_square_raster <- function(vals, crs=default_proj) {
  n <- sqrt(length(vals))

  stopifnot(as.integer(n) == n)

  raster::raster(matrix(vals, nrow=n, byrow=TRUE),
                 xmn=0, xmx=n, ymn=0, ymx=n,
                 crs=crs)
}
