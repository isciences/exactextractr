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

if (!isGeneric("exact_extract")) {
	setGeneric("exact_extract", function(x, y, ...)
		standardGeneric("exact_extract"))
}

#' Extract values from RasterLayers
#'
#' Extracts the values of cells in a RasterLayer that are covered by a
#' simple feature collection containing polygonal geometries.
#'
#' By default, returns a list with one matrix for each feature in the input
#' feature collection. Each matrix contains two columns, with the first column
#' containing the values of cells that are touched by the feature's polygon,
#' and the second column containing the fraction of the cell (0-1) that is
#' covered by the polygon.
#'
#' If a function \code{fun} is supplied, it will be called for each feature
#' with vectors of cell values and weights as arguments. \code{exact_extract}
#' will then return a list of the return values of \code{fun} instead of a
#' list of matrices.
#'
#' The value of \code{fun} may also set to a string representing a common
#' statistical summary function supported by the exactextract library.
#' Supported statistical functions include:
#' \itemize{
#'  \item{\code{min} - the minimum defined value in any raster cell wholly or
#'                     partially covered by the polygon}
#'  \item{\code{max} - the maximum defined value in any raster cell wholly or
#'                     partially covered by the polygon}
#'  \item{\code{count} - the sum of fractions of raster cells with defined values
#'                       covered by the polygon}
#'  \item{\code{sum}   - the sum of defined raster cell values, multiplied by
#'                       the fraction of the cell that is covered by the polygon}
#'  \item{\code{mean} - the mean cell value, weighted by the fraction of each cell
#'                      that is covered by the polygon}
#'  \item{\code{mode} - the most common cell value, weighted by the fraction of
#'                      each cell that is covered by the polygon. Where multiple
#'                      values occupy the same maximum number of weighted cells,
#'                      the largest value will be returned.}
#'  \item{\code{minority} - the least common cell value, weighted by the fraction
#'                          of each cell that is covered by the polygon. Where
#'                          multiple values occupy the same minimum number of
#'                          weighted cells, the smallest value will be returned.}
#'  \item{\code{variety} - the number of distinct values in cells that are wholly
#'                         or partially covered by the polygon.}
#' }
#'
#' @param     x a RasterLayer
#' @param     y a sf object with polygonal geometries
#' @param     fun an optional function or character vector, as described below
#' @param     ... additional arguments to pass to \code{fun}
#' @import sf
#' @import raster
#' @useDynLib exactextractr
#' @export
setMethod('exact_extract', signature(x='RasterLayer', y='sf'), function(x, y, fun=NULL, ...) {
  exact_extract(x, sf::st_geometry(y), fun, ...)
})

.exact_extract <- function(x, y, fun=NULL, ...) {
  if (is.null(fun)) {
    appfn <- lapply # return list of matrices
  } else {
    appfn <- sapply # return vector of function results
  }

  raster_extent <- as.vector(raster::extent(x))
  raster_res <- raster::res(x)

  if (is.character(fun)) {
    vals <- raster::as.matrix(x)

    appfn(sf::st_as_binary(y), function(wkb) {
      CPP_stat(raster_extent, raster_res, vals, fun, wkb)
    })
  } else {
    appfn(sf::st_as_binary(y), function(wkb) {
      ret <- CPP_exact_extract(raster_extent, raster_res, wkb)

      vals <- raster::getValuesBlock(x,
                                     row=ret$row,
                                     col=ret$col,
                                     nrow=nrow(ret$weights),
                                     ncol=ncol(ret$weights))

      weightvec <- as.vector(t(ret$weights))

      if (!is.null(fun)) {
        return(fun(vals[weightvec > 0], weightvec[weightvec > 0], ...))
      } else {
        return(cbind(vals=vals[weightvec > 0], weights=weightvec[weightvec > 0]))
      }
    })
  }
}

#' @useDynLib exactextractr
#' @export
setMethod('exact_extract', signature(x='RasterLayer', y='sfc_MULTIPOLYGON'), .exact_extract)

#' @useDynLib exactextractr
#' @export
setMethod('exact_extract', signature(x='RasterLayer', y='sfc_POLYGON'), .exact_extract)

