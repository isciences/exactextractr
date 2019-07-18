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

if (!isGeneric("exact_extract")) {
	setGeneric("exact_extract", function(x, y, ...)
		standardGeneric("exact_extract"))
}

#' Extract or summarize values from Raster* objects
#'
#' Extracts the values of cells in a RasterLayer that are covered by a
#' simple feature collection containing polygonal geometries. Returns either
#' the result of a summary operation or function applied to the values (if
#' \code{fun} is specified), or the values themselves (if \code{fun} is
#' \code{NULL}.)
#'
#' The value of \code{fun} may be set to a string (or vector of strings)
#' representing summary operations supported by the exactextract library.
#' In this case, \code{exact_extract} will return a vector with the result
#' of the summary operation for each feature in the input.
#'
#' The following summary operations are supported:
#'
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
#'  \item{\code{majority} - synonym for \code{mode}}
#'  \item{\code{minority} - the least common cell value, weighted by the fraction
#'                          of each cell that is covered by the polygon. Where
#'                          multiple values occupy the same minimum number of
#'                          weighted cells, the smallest value will be returned.}
#'  \item{\code{variety} - the number of distinct values in cells that are wholly
#'                         or partially covered by the polygon.}
#' }
#'
#' Alternatively, an R function may be provided as \code{fun}. The function will be
#' called for each feature with  with vectors of cell values and weights as arguments.
#' \code{exact_extract} will then return a vector of the return values of \code{fun}.
#'
#' If \code{fun} is not specified, \code{exact_extract} will return a list with
#' one data frame for each feature in the input feature collection. The data
#' frame contains a column with values from each layer in the input `Raster*`,
#' and a final column indicating the fraction of the cell that is covered by the
#' polygon.
#'
#' @param     x a RasterLayer
#' @param     y a sf object with polygonal geometries
#' @param     include_xy if \code{TRUE}, augment the returned data frame with
#'                        columns for cell center coordinates (\code{x} and
#'                        \code{y}) or pass them to \code{fun}
#' @param     fun an optional function or character vector, as described above
#' @param     max_cells_in_memory the maximum number of raster cells to load at
#'                                a given time when using a named summary operation
#'                                for \code{fun} (as opposed to a function defined using
#'                                R code). If a polygon covers more than \code{max_cells_in_memory}
#'                                raster cells, it will be processed in multiple chunks.
#' @param     progress show progress bar
#' @param     ... additional arguments to pass to \code{fun}
#' @name exact_extract
NULL

#' @import sf
#' @import raster
#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='sf'), function(x, y, fun=NULL, ..., include_xy=FALSE, progress=TRUE, max_cells_in_memory=30000000) {
  exact_extract(x, sf::st_geometry(y), fun=fun, ..., include_xy=include_xy, progress=progress, max_cells_in_memory=max_cells_in_memory)
})

# Return the number of standard (non-...) arguments in a supplied function that
# do not have a default value. This is used to fail if the summary function
# provided by the user cannot accept arguments of values and weights.
.num_expected_args <- function(fun) {
  a <- formals(args(fun))
  a <- a[names(a) != '...']
  sum(sapply(a, nchar) == 0)
}

.exact_extract <- function(x, y, fun=NULL, ..., include_xy=FALSE, progress=TRUE, max_cells_in_memory=30000000) {
  if(is.na(sf::st_crs(x)) && !is.na(sf::st_crs(y))) {
    warning("No CRS specified for raster; assuming it has the same CRS as the polygons.")
  } else if(is.na(sf::st_crs(y)) && !is.na(sf::st_crs(x))) {
    warning("No CRS specified for polygons; assuming they have the same CRS as the raster.")
  } else if(sf::st_crs(x) != sf::st_crs(y)) {
    old_crs <- sf::st_crs(y)
    y <- sf::st_transform(y, sf::st_crs(x))
    warning("Polygons transformed from EPSG:", old_crs$epsg, " to EPSG:", sf::st_crs(x)$epsg)
  }

  if (!is.null(fun) && !is.character(fun) && .num_expected_args(fun) < 2) {
    stop("exact_extract was called with a function that does not appear to ",
         "be of the form `function(values, coverage_fractions, ...)`")
  }

  if (is.character(fun) && length(list(...)) > 0) {
    stop("exact_extract was called with a named summary operation that",
         "does not accept additional arguments ...")
  }

  raster_extent <- as.vector(raster::extent(x))
  raster_res <- raster::res(x)

  if (progress && length(y) > 1) {
    n <- length(y)
    pb <- utils::txtProgressBar(min = 0, max = n, initial=0, style=3)
    update_progress <- function() {
      i <- 1 + utils::getTxtProgressBar(pb)
      utils::setTxtProgressBar(pb, i)
      if (i == n) {
        close(pb)
      }
    }
  } else {
    update_progress <- function() {}
  }

  tryCatch({
    x <- readStart(x)

    if (is.character(fun)) {
      if (raster::nlayers(x) > 1) stop("Predefined summary operations only available for single-layer rasters. Please define a summary function using R code.")

      results <- sapply(sf::st_as_binary(y), function(wkb) {
        update_progress()
        CPP_stats(x, wkb, fun, max_cells_in_memory)
      })

      if (length(fun) > 1) {
        # Return a data frame with a column for each stat
        results <- t(results)
        dimnames(results) <- list(NULL, fun)
        return(as.data.frame(results))
      } else {
        # Just return a vector of stat results
        return(results)
      }
    } else {
      if (is.null(fun)) {
        appfn <- lapply # return list of data frames
      } else {
        appfn <- sapply
      }

      appfn(sf::st_as_binary(y), function(wkb) {
        ret <- CPP_exact_extract(raster_extent, raster_res, wkb)

        vals <- raster::getValuesBlock(x,
                                       row=ret$row,
                                       col=ret$col,
                                       nrow=nrow(ret$weights),
                                       ncol=ncol(ret$weights))

        if(is.matrix(vals)) {
          vals <- as.data.frame(vals)
        } else {
          vals <- data.frame(value=vals)
        }

        if (include_xy) {
            x_coords <- raster::xFromCol(x, col=ret$col:(ret$col+ncol(ret$weights) - 1))
            y_coords <- raster::yFromRow(x, row=ret$row:(ret$row+nrow(ret$weights) - 1))

            vals$x <- rep.int(x_coords, times=nrow(ret$weights))
            vals$y <- rep(y_coords, each=ncol(ret$weights))
        }

        cov_fracs <- as.vector(t(ret$weights))
        vals <- vals[cov_fracs > 0, , drop=FALSE]
        cov_fracs <- cov_fracs[cov_fracs > 0]

        update_progress()

        if (is.null(fun)) {
          vals$coverage_fraction <- cov_fracs
          return(vals)
        } else {
          if (ncol(vals) == 1) {
            return(fun(vals[,1], cov_fracs, ...))
          } else {
            return(fun(vals, cov_fracs, ...))
          }
        }
      })
    }
  }, finally={
    readStop(x)
  })
}

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='sfc_MULTIPOLYGON'), .exact_extract)

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='sfc_POLYGON'), .exact_extract)

