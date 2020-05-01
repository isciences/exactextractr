# Copyright (c) 2018-2020 ISciences, LLC.
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

if (!isGeneric('coverage_fraction')) {
	setGeneric('coverage_fraction', function(x, y, crop=FALSE, ...)
		standardGeneric('coverage_fraction'))
}

.coverage_fraction <- function(x, y, crop) {
  if(is.na(sf::st_crs(x)) && !is.na(sf::st_crs(y))) {
    warning("No CRS specified for raster; assuming it has the same CRS as the polygons.")
  } else if(is.na(sf::st_crs(y)) && !is.na(sf::st_crs(x))) {
    warning("No CRS specified for polygons; assuming they have the same CRS as the raster.")
  } else if(sf::st_crs(x) != sf::st_crs(y)) {
    y <- sf::st_transform(y, sf::st_crs(x))
    warning("Polygons transformed to raster CRS (EPSG:", sf::st_crs(x)$epsg, ")")
  }

  lapply(sf::st_as_binary(y, EWKB=TRUE), function(wkb) {
    CPP_coverage_fraction(x, wkb, crop)
  })
}

#' Compute the fraction of raster cells covered by a polygon
#'
#' @param     x a (possibly empty) \code{RasterLayer} whose resolution and
#'            extent will be used for the generated \code{RasterLayer}.
#' @param     y a \code{sf} object with polygonal geometries
#' @param     crop if \code{TRUE}, each generated \code{RasterLayer} will be
#'                 cropped to the extent of its associated feature.
#' @return    a list with a \code{RasterLayer} for each feature in \code{y}.
#'            Values of the raster represent the fraction of each
#'            cell in \code{x} that is covered by \code{y}.
#' @examples
#' rast <- raster::raster(matrix(1:100, ncol=10), xmn=0, ymn=0, xmx=10, ymx=10)
#' poly <- sf::st_as_sfc('POLYGON ((2 2, 7 6, 4 9, 2 2))')
#'
#' cov_frac <- coverage_fraction(rast, poly)[[1]]
#' @name coverage_fraction
NULL

#' @import sf
#' @import raster
#' @useDynLib exactextractr
#' @rdname coverage_fraction
#' @export
setMethod('coverage_fraction', signature(x='RasterLayer', y='sf'), function(x, y, crop=FALSE) {
  coverage_fraction(x, sf::st_geometry(y), crop)
})

#' @rdname coverage_fraction
#' @export
setMethod('coverage_fraction', signature(x='RasterLayer', y='sfc_MULTIPOLYGON'), .coverage_fraction)

#' @rdname coverage_fraction
#' @export
setMethod('coverage_fraction', signature(x='RasterLayer', y='sfc_POLYGON'), .coverage_fraction)

