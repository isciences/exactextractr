# Copyright (c) 2018-2019 ISciences, LLC.
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

if (!isGeneric('partial_mask')) {
	setGeneric('partial_mask', function(x, y, ...)
		standardGeneric('partial_mask'))
}

.partial_mask <- function(x, y) {
  lapply(sf::st_as_binary(y), function(wkb) {
    out <- raster::raster(x) # copy input dims, res, etc.

    raster::values(out) <- CPP_weights(as.vector(raster::extent(x)), raster::res(x), wkb)

    return(out)
  })
}

#' Compute the fraction of raster cells covered by a polygon
#'
#' @param     x a RasterLayer
#' @param     y a sf object with polygonal geometries
#' @return    a list with a RasterLayer for each feature in \code{y}.
#'            Values of the raster represent the fraction of each
#'            cell in \code{x} that is covered by \code{y}.
#' @name partial_mask
NULL

#' @import sf
#' @import raster
#' @useDynLib exactextractr
#' @rdname partial_mask
#' @export
setMethod('partial_mask', signature(x='RasterLayer', y='sf'), function(x, y) {
  partial_mask(x, sf::st_geometry(y))
})

#' @rdname partial_mask
#' @export
setMethod('partial_mask', signature(x='RasterLayer', y='sfc_MULTIPOLYGON'), .partial_mask)

#' @rdname partial_mask
#' @export
setMethod('partial_mask', signature(x='RasterLayer', y='sfc_POLYGON'), .partial_mask)

