# Copyright (c) 2020-2022 ISciences, LLC.
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

setGeneric("exact_resample", function(x, y, ...)
	standardGeneric("exact_resample"))

#' Resample a raster to a new grid
#'
#' @param x a \code{RasterLayer} or \code{SpatRaster} to be resampled
#' @param y a raster of the same class as \code{x} with a grid definition to
#'          which \code{x} should be resampled
#' @param fun a named summary operation or R function to be used for the resampling
#' @param coverage_area use cell coverage areas instead of coverage fractions
#'                      in \code{fun}
#' @return a resampled version of \code{x}, returned as a \code{RasterLayer} or
#'         \code{SpatRaster}, depending on the values of \code{x} and \code{y}
#'
#' @name exact_resample
NULL

.exact_resample <- function(x, y, fun, coverage_area = FALSE) {
  analysis_crs <- sf::st_crs(x)
  if (sf::st_crs(x) != sf::st_crs(y)) {
    if (is.na(sf::st_crs(x))) {
      warning("No CRS specified for source raster; assuming it has the same CRS as destination raster.")
      analysis_crs <- sf::st_crs(y)
    } else if (is.na(sf::st_crs(y))) {
      warning("No CRS specified for destination raster; assuming it has the same CRS as source raster.")
    } else {
      stop('Destination raster must have same CRS as source.')
    }
  }

  if (is.character(fun)) {
    if (length(fun) != 1) {
      stop("Only a single operation may be used for resampling.")
    }

    if (startsWith(fun, 'weighted')) {
      stop("Weighted operations cannot be used for resampling.")
    }

    if (.numLayers(x) != 1) {
      stop("Raster to be resampled must have a single layer.")
    }

    summary_stat <- fun
    summary_fun <- NULL
  } else {
    if (!is.function(fun)) {
      stop("fun must be a named summary operation or an R function")
    }

    if (.num_expected_args(fun) != 2) {
      stop("exact_extract was called with a function that does not appear to ",
           "be of the form `function(values, coverage_fractions)`.")
    }

    summary_stat <- NULL
    summary_fun <- fun
  }

  .validateFlag(coverage_area, 'coverage_area')

  area_method <- .areaMethod(analysis_crs)

  x <- .startReading(x)
  tryCatch({
    ret <- CPP_resample(x,
                        y,
                        summary_stat,
                        summary_fun,
                        coverage_area,
                        area_method)
    if (inherits(x, 'SpatRaster')) {
      terra::rast(ret)
    } else {
      ret
    }
  }, finally={
    .stopReading(x)
  })
}

#' @import sf
#' @import raster
#' @useDynLib exactextractr
#' @rdname exact_resample
#' @export
setMethod('exact_resample', signature(x='RasterLayer', y='RasterLayer'), .exact_resample)

#' @useDynLib exactextractr
#' @rdname exact_resample
#' @export
setMethod('exact_resample', signature(x='SpatRaster', y='SpatRaster'), .exact_resample)
