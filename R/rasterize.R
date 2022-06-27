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

setGeneric("rasterize_polygons", function(x, y, ...)
	standardGeneric("rasterize_polygons"))

#' Create a raster approximation of a polygon coverage
#'
#' Returns a raster whose values indicate the index of the polygon covering
#' each cell. Where multiple polygons cover the same cell, the index of the
#' polygon covering the greatest area will be used, with the lowest index
#' returned in the case of ties. Cells that are not covered by any polygon,
#' or whose total covered fraction is less than \code{min_coverage}, will be
#' set to \code{NA}.
#'
#' @param     x a \code{sf} or \code{sfc} object with polygonal geometries
#'
#' @param     y a (possibly empty) \code{RasterLayer} whose resolution and
#'            extent will be used for the generated \code{RasterLayer}.
#' @param     min_coverage minimum fraction of a cell that must be covered by
#'                         polygons to be included in the output
#' @return a \code{RasterLayer} or \code{SpatRaster}, consistent with the type of \code{y}
#' @name rasterize_polygons
NULL

.rasterize_polygons <- function(x, y, min_coverage = 0) {
  num_rows <- terra::nrow(y)
  num_cols <- terra::ncol(y)

  if (min_coverage == 1.0) {
    # account for error in sum of single-precision plots
    min_coverage <- (min_coverage - 1e-6)
  }

  max_coverage_values <- matrix(0.0, nrow = num_rows, ncol = num_cols)
  max_coverage_indexes <- matrix(NA_integer_, nrow = num_rows, ncol = num_cols)
  tot_coverage <- matrix(0.0, nrow = num_rows, ncol = num_cols)

  dest_extent <- .extent(y)
  dest_res <- .res(y)

  geoms <- sf::st_geometry(x)

  for (i in seq_along(geoms)) {
    wkb <- sf::st_as_binary(geoms[[i]], EWKB=TRUE)

    CPP_update_max_coverage(dest_extent, dest_res, max_coverage_values, max_coverage_indexes, tot_coverage, wkb, i)
  }

  max_coverage_indexes[tot_coverage < min_coverage] <- NA

  if (inherits(y, 'SpatRaster')) {
    ret <- terra::rast(y)
    terra::values(ret) <- max_coverage_indexes
  } else {
    ret <- raster::raster(y)
    raster::values(ret) <- max_coverage_indexes
  }

  return(ret)
}

#' @import sf
#' @import raster
#' @useDynLib exactextractr
#' @rdname rasterize_polygons
#' @export
setMethod('rasterize_polygons', signature(x='sf', y='RasterLayer'), .rasterize_polygons)

#' @useDynLib exactextractr
#' @rdname rasterize_polygons
#' @export
setMethod('rasterize_polygons', signature(x='sf', y='SpatRaster'), .rasterize_polygons)
