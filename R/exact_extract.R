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
#' If \code{fun} is not specified, \code{exact_extract} will return a list
#' with one matrix for each feature in the input feature collection. The matrix
#' contains a column of cell values from each layer in the input `Raster*`, and
#' a final column indicating the fraction of the cell that is covered by the
#' polygon.
#'
#' @param     x a RasterLayer
#' @param     y a sf object with polygonal geometries
#' @param     include_xy if \code{TRUE}, augmente the returned matrix with
#'                        columns for cell center coordinates (\code{x} and
#'                        \code{y}) or pass them to \code{fun}
#' @param     fun an optional function or character vector, as described above
#' @param     ... additional arguments to pass to \code{fun}
#' @name exact_extract
NULL

#' @import sf
#' @import raster
#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='sf'), function(x, y, fun=NULL, ..., include_xy=FALSE) {
  exact_extract(x, sf::st_geometry(y), fun=fun, ..., include_xy=include_xy)
})

# Return the number of standard (non-...) arguments in a supplied function that
# do not have a default value. This is used to fail if the summary function
# provided by the user cannot accept arguments of values and weights.
.num_expected_args <- function(fun) {
  a <- formals(args(fun))
  a <- a[names(a) != '...']
  sum(sapply(a, nchar) == 0)
}

.exact_extract <- function(x, y, fun=NULL, ..., include_xy=FALSE) {
  if(sf::st_crs(x) != sf::st_crs(y)) {
    stop("Raster and polygons must be in the same coordinate reference system.")
  }

  if (is.null(fun)) {
    appfn <- lapply # return list of matrices
  } else {
    appfn <- sapply

    if (!is.character(fun) && .num_expected_args(fun) < 2) {
      stop("exact_extract was called with a function that does not appear to ",
           "be of the form `function(values, coverage_fractions, ...)`")
    }
  }

  raster_extent <- as.vector(raster::extent(x))
  raster_res <- raster::res(x)

  ret <- tryCatch({
    x <- readStart(x)

    if (is.character(fun)) {
      if (raster::nlayers(x) > 1) stop("Predefined summary operations only available for single-layer rasters. Please define a summary function using R code.")

      appfn(sf::st_as_binary(y), function(wkb) {
        CPP_stats(x, wkb, fun)
      })
    } else {
      appfn(sf::st_as_binary(y), function(wkb) {
        ret <- CPP_exact_extract(raster_extent, raster_res, wkb)

        vals <- raster::getValuesBlock(x,
                                       row=ret$row,
                                       col=ret$col,
                                       nrow=nrow(ret$weights),
                                       ncol=ncol(ret$weights))

        if (!is.matrix(vals)) {
          vals <- matrix(vals)
          dimnames(vals)[[2]] <- list('values')
        }

        if (include_xy) {
            x_coords <- raster::xFromCol(x, col=ret$col:(ret$col+ncol(ret$weights) - 1))
            y_coords <- raster::yFromRow(x, row=ret$row:(ret$row+nrow(ret$weights) - 1))

            vals <- cbind(vals,
                          x=rep.int(x_coords, times=nrow(ret$weights)),
                          y=rep(y_coords, each=ncol(ret$weights)))
        }

        weightvec <- as.vector(t(ret$weights))

        if (!is.null(fun)) {
          return(fun(vals[weightvec > 0,, drop=FALSE], weightvec[weightvec > 0], ...))
        } else {
          return(cbind(vals[weightvec > 0,, drop=FALSE], weights=weightvec[weightvec > 0]))
        }
      })
    }
  }, finally={
    readStop(x)
  })

  if (is.matrix(ret)) {
    t(ret)
  } else {
    ret
  }
}

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='sfc_MULTIPOLYGON'), .exact_extract)

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='sfc_POLYGON'), .exact_extract)

