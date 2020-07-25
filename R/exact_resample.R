# Copyright (c) 2020 ISciences, LLC.
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

if (!isGeneric("exact_resample")) {
	setGeneric("exact_resample", function(x, y, ...)
		standardGeneric("exact_resample"))
}

#' Resample a raster to a new grid
#'
#' @param x a \code{RasterLayer} to be resampled
#' @param y a \code{RasterLayer} with a grid definition to which \code{x}
#'          should be resampled
#' @param fun a named summary operation to be used for the resampling
#' @return a resampled version of \code{x}
#'
#' @name exact_resample
NULL

#' @import sf
#' @import raster
#' @useDynLib exactextractr
#' @rdname exact_resample
#' @export
setMethod('exact_resample',
          signature(x='RasterLayer', y='RasterLayer'),
          function(x, y, fun) {
            if (sf::st_crs(x) != sf::st_crs(y)) {
              if (is.na(sf::st_crs(x))) {
                warning("No CRS specified for source raster; assuming it has the same CRS as destination raster.")
              } else if (is.na(sf::st_crs(y))) {
                warning("No CRS specified for destination raster; assuming it has the same CRS as source raster.")
              } else {
                stop('Destination raster must have same CRS as source.')
              }
            }

            x <- raster::readStart(x)
            tryCatch({
              CPP_resample(x, y, fun)
            }, finally={
              raster::readStop(x)
            })
          })
