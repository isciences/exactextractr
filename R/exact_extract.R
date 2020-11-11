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

if (!isGeneric("exact_extract")) {
	setGeneric("exact_extract", function(x, y, ...)
		standardGeneric("exact_extract"))
}

#' Extract or summarize values from Raster* objects
#'
#' Extracts the values of cells in a Raster* that are covered by a
#' simple feature collection containing polygonal geometries, as well as the
#' fraction of each cell that is covered by the polygon. Returns either
#' the result of a summary operation or function applied to the values
#' and coverage fractions (if \code{fun} is specified), or a data frame
#' containing the values and coverage fractions themselves (if \code{fun}
#' is \code{NULL}.)
#'
#' The value of \code{fun} may be set to a string (or vector of strings)
#' representing summary operations supported by the exactextract library.
#' If the input raster has a single layer and a single summary operation
#' is specified, \code{exact_extract} will return a vector with the result
#' of the summary operation for each feature in the input. If the input
#' raster has multiple layers, or if multiple summary operations are specified,
#' \code{exact_extract} will return a data frame with a row for each feature
#' and a column for each summary operation / layer combination. (The
#' \code{force_df} can be used to always return a data frame instead of a vector.)
#' In all of the summary operations, \code{NA} values in the raster are ignored
#' (i.e., \code{na.rm = TRUE}.)
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
#'  \item{\code{median} - the median cell value, weighted by the fraction of each
#'                        cell that is covered by the polygon}
#'  \item{\code{quantile} - arbitrary quantile(s) of cell values, specified in
#'                          \code{quantiles}, weighted by the fraction of each
#'                          cell that is covered by the polygon}
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
#'  \item{\code{variance} - the population variance of cell values, weighted by the
#'                          fraction of each cell that is covered by the polygon.}
#'  \item{\code{stdev} - the population standard deviation of cell values, weighted
#'                       by the fraction of each cell that is covered by the polygon.}
#'  \item{\code{coefficient_of_variation} - the population coefficient of variation of
#'                       cell values, weighted by the fraction of each cell that is
#'                       covered by the polygon.}
#'  \item{\code{weighted_mean} - the mean cell value, weighted by the product of
#'                               the fraction of each cell covered by the polygon
#'                               and the value of a second weighting raster provided
#'                               as \code{weights}}
#'  \item{\code{weighted_sum} - the sum of defined raster cell values, multiplied by
#'                              the fraction of each cell that is covered by the polygon
#'                               and the value of a second weighting raster provided
#'                               as \code{weights}}
#' }
#'
#' Alternatively, an R function may be provided as \code{fun}. The function will be
#' called for each feature with with vectors of cell values and weights as arguments.
#' \code{exact_extract} will then return a vector of the return values of \code{fun}.
#'
#' If \code{fun} is not specified, \code{exact_extract} will return a list with
#' one data frame for each feature in the input feature collection. The data
#' frame will contain a column with values from each layer in the input `Raster*`,
#' and a final column indicating the fraction of the cell that is covered by the
#' polygon.
#'
#' @param     x a \code{RasterLayer}, \code{RasterStack}, or \code{RasterBrick}
#' @param     y a sf object with polygonal geometries
#' @param     fun an optional function or character vector, as described below
#' @param     weights  a weighting raster to be used with the \code{weighted_mean}
#'                     and \code{weighted_sum} summary operations.
#' @param     quantiles   quantiles to be computed when \code{fun == 'quantile'}
#' @param     append_cols when \code{fun} is not \code{NULL}, an optional
#'                        character vector of columns from \code{y} to be
#'                        included in returned data frame.
#' @param     force_df always return a data frame instead of a vector, even if
#'                     \code{x} has only one layer and \code{fun} has length 1
#' @param     full_colnames include the names of \code{x} in the names of the
#'                          returned data frame, even if \code{x} has only one
#'                          layer. This is useful when the results of multiple
#'                          calls to \code{exact_extract} are combined with
#'                          \code{cbind}.
#' @param     include_cell if \code{TRUE}, and \code{fun} is \code{NULL}, augment
#'                       the returned data frame for each feature with a column
#'                       for the cell index (\code{cell}). If \code{TRUE} and
#'                       \code{fun} is not \code{NULL}, add \code{cell} to the
#'                       data frame passed to \code{fun} for each feature.
#' @param     include_cols an optional character vector of column names in
#'                         \code{y} to be added to the data frame for each
#'                         feature that is either returned (when \code{fun} is
#'                         \code{NULL}) or passed to \code{fun}.
#' @param     include_xy if \code{TRUE}, and \code{fun} is \code{NULL}, augment
#'                       the returned data frame for each feature with columns
#'                       for cell center coordinates (\code{x} and \code{y}). If
#'                       \code{TRUE} and \code{fun} is not \code{NULL}, add
#'                       \code{x} and {y} to the data frame passed to \code{fun}
#'                       for each feature.
#' @param     stack_apply   if \code{TRUE}, apply \code{fun} to each layer of
#'                          \code{x} independently. If \code{FALSE}, apply \code{fun}
#'                          to all layers of \code{x} simultaneously.
#' @param     max_cells_in_memory the maximum number of raster cells to load at
#'                                a given time when using a named summary operation
#'                                for \code{fun} (as opposed to a function defined using
#'                                R code). If a polygon covers more than \code{max_cells_in_memory}
#'                                raster cells, it will be processed in multiple chunks.
#' @param     progress if \code{TRUE}, display a progress bar during processing
#' @param     ... additional arguments to pass to \code{fun}
#' @return a vector or list of data frames, depending on the type of \code{x} and the
#'         value of \code{fun} (see Details)
#' @examples
#' rast <- raster::raster(matrix(1:100, ncol=10), xmn=0, ymn=0, xmx=10, ymx=10)
#' poly <- sf::st_as_sfc('POLYGON ((2 2, 7 6, 4 9, 2 2))')
#'
#' # named summary operation on RasterLayer, returns vector
#' exact_extract(rast, poly, 'mean')
#'
#' # two named summary operations on RasterLayer, returns data frame
#' exact_extract(rast, poly, c('min', 'max'))
#'
#' # named summary operation on RasterStack, returns data frame
#' stk <- raster::stack(list(a=rast, b=sqrt(rast)))
#' exact_extract(stk, poly, 'mean')
#'
#' # named weighted summary operation, returns vector
#' weights <- raster::raster(matrix(runif(100), ncol=10), xmn=0, ymn=0, xmx=10, ymx=10)
#' exact_extract(rast, poly, 'weighted_mean', weights=weights)
#'
#' # custom summary function, returns vector
#' exact_extract(rast, poly, function(value, cov_frac) length(value[cov_frac > 0.9]))
#'
#' @name exact_extract
NULL

#' @import sf
#' @import raster
#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='sf'),
          function(x, y, fun=NULL, ...,
                   include_xy=FALSE,
                   progress=TRUE,
                   max_cells_in_memory=30000000,
                   include_cell=FALSE,
                   force_df=FALSE,
                   full_colnames=FALSE,
                   stack_apply=FALSE,
                   append_cols=NULL,
                   include_cols=NULL,
                   quantiles=NULL) {
  .exact_extract(x, y, fun=fun, ...,
                include_xy=include_xy,
                progress=progress,
                max_cells_in_memory=max_cells_in_memory,
                include_cell=include_cell,
                force_df=force_df,
                full_colnames=full_colnames,
                stack_apply=stack_apply,
                append_cols=append_cols,
                include_cols=include_cols,
                quantiles=quantiles)
})

# Return the number of standard (non-...) arguments in a supplied function that
# do not have a default value. This is used to fail if the summary function
# provided by the user cannot accept arguments of values and weights.
.num_expected_args <- function(fun) {
  a <- formals(args(fun))
  a <- a[names(a) != '...']
  sum(sapply(a, nchar) == 0)
}

emptyVector <- function(rast) {
  switch(substr(raster::dataType(rast), 1, 3),
         LOG=logical(),
         INT=integer(),
         numeric())
}

.exact_extract <- function(x, y, fun=NULL, ...,
                           weights=NULL,
                           include_xy=FALSE,
                           progress=TRUE,
                           max_cells_in_memory=30000000,
                           include_cell=FALSE,
                           force_df=FALSE,
                           full_colnames=FALSE,
                           stack_apply=FALSE,
                           append_cols=NULL,
                           include_cols=NULL,
                           quantiles=NULL) {
  if(!is.null(append_cols)) {
    if (!inherits(y, 'sf')) {
      stop(sprintf('append_cols only supported for sf arguments (received %s)',
                   paste(class(y), collapse = ' ')))
    }

    force_df <- TRUE
  }

  if(sf::st_geometry_type(y, by_geometry = FALSE) == 'GEOMETRY') {
    if (!all(sf::st_dimension(y) == 2)) {
      stop("Features in sfc_GEOMETRY must be polygonal")
    }
    y <- sf::st_cast(y, 'MULTIPOLYGON')
  }

  if(!is.null(weights)) {
    if (!startsWith(class(weights), 'Raster')) {
      stop("Weights must be a Raster object.")
    }

    if (!is.character(fun)) {
      stop("Weighting raster can only be used with named summary operations.")
    }

    if (!any(startsWith(fun, "weighted"))) {
      warning("Weights provided but no requested operations use them.")
    }

    if (!is.na(sf::st_crs(x))) {
      if (is.na(sf::st_crs(weights))) {
        warning("No CRS specified for weighting raster; assuming it has the same CRS as the value raster.")
      } else if (sf::st_crs(x) != sf::st_crs(weights)) {
        stop("Weighting raster does not have the same CRS as value raster.")
      }
    }
  }

  if(is.na(sf::st_crs(x)) && !is.na(sf::st_crs(y))) {
    warning("No CRS specified for raster; assuming it has the same CRS as the polygons.")
  } else if(is.na(sf::st_crs(y)) && !is.na(sf::st_crs(x))) {
    warning("No CRS specified for polygons; assuming they have the same CRS as the raster.")
  } else if(sf::st_crs(x) != sf::st_crs(y)) {
    y <- sf::st_transform(y, sf::st_crs(x))
    warning("Polygons transformed to raster CRS (EPSG:", sf::st_crs(x)$epsg, ")")
  }

  if (!is.null(fun) && !is.character(fun) && .num_expected_args(fun) < 2) {
    stop("exact_extract was called with a function that does not appear to ",
         "be of the form `function(values, coverage_fractions, ...)`")
  }

  if (is.character(fun)) {
    if (length(list(...)) > 0) {
      stop("exact_extract was called with a named summary operation that ",
           "does not accept additional arguments ...")
    }

    if (include_xy) {
      stop("include_xy must be FALSE for named summary operations")
    }

    if (include_cell) {
      stop("include_cell must be FALSE for named summary operations")
    }

    if (!is.null(include_cols)) {
      stop("include_cols not supported for named_summary operations (see argument append_cols)")
    }
  }

  geoms <- sf::st_geometry(y)

  if (progress && length(geoms) > 1) {
    n <- length(geoms)
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
    x <- raster::readStart(x)
    if (!is.null(weights)) {
      weights <- raster::readStart(weights)
    }

    if (is.character(fun)) {
      # Compute all stats in C++
      results <- sapply(sf::st_as_binary(geoms, EWKB=TRUE), function(wkb) {
        ret <- CPP_stats(x, weights, wkb, fun, max_cells_in_memory, quantiles)
        update_progress()
        return(ret)
      })

      if (length(fun) == 1 &&
          raster::nlayers(x) == 1 &&
          (is.null(quantiles) || length(quantiles) == 1) &&
          !force_df) {
        # Just return a vector of stat results
        return(as.vector(results))
      } else {
        # Return a data frame with a column for each stat
        colnames <- .cppStatColNames(x, fun, full_colnames, quantiles)

        if (is.matrix(results)) {
          results <- t(results)
        } else {
          results <- matrix(results, nrow=length(results))
        }

        dimnames(results) <- list(NULL, colnames)
        ret <- as.data.frame(results)

        if (!is.null(append_cols)) {
          ret <- cbind(sf::st_drop_geometry(y[, append_cols]), ret)
        }

        return(ret)
      }
    } else {
      ret <- lapply(seq_along(geoms), function(feature_num) {
        wkb <- sf::st_as_binary(geoms[[feature_num]], EWKB=TRUE)

        ret <- CPP_exact_extract(x, wkb)

        if (length(ret$weights) > 0) {
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
        } else {
          # Polygon does not intersect raster.
          # Construct a zero-row data frame with correct column names/types.
          vals <- do.call(data.frame, lapply(seq_len(raster::nlayers(x)),
                                                     function(i) emptyVector(x[[i]])))
          if (raster::nlayers(x) == 1) {
            names(vals) <- 'value'
          } else {
            names(vals) <- names(x)
          }
        }

        if (include_xy) {
           vals <- .appendXY(vals, x, ret$row, nrow(ret$weights), ret$col, ncol(ret$weights))
        }

        if (include_cell) {
          vals <- .appendCell(vals, x, ret$row, nrow(ret$weights), ret$col, ncol(ret$weights))
        }

        if (!is.null(include_cols)) {
          # use vals as first argument to cbind, then rearrange names so that
          # include_cols come first
          vals <- cbind(vals,
                        sf::st_drop_geometry(y[feature_num, include_cols]),
                        row.names = NULL)[, c(include_cols, names(vals))]
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
            # Only one layer, nothing appended (cells or XY)
            return(fun(vals[,1], cov_fracs, ...))
          } else {
            if (stack_apply) {
              # Pass each layer in stack to callback individually
              nlay <- raster::nlayers(x)
              appended_cols <- seq_len(ncol(vals))[-seq_len(nlay)]

              if (length(appended_cols) == 0) {
                result <- lapply(seq_len(nlay), function(z)
                  fun(vals[, z], cov_fracs, ...))
              } else {
                result <- lapply(seq_len(nlay), function(z)
                  fun(cbind(data.frame(value=vals[, z]),
                            vals[, c(appended_cols)]), cov_fracs, ...))
              }

              names(result) <- paste('fun', names(x), sep='.')
              return(do.call(data.frame, result))
            } else {
              # Pass all layers to callback, to be handled together
              return(fun(vals, cov_fracs, ...))
            }
          }
        }
      })

      if (!is.null(fun)) {
        if (all(sapply(ret, is.data.frame))) {
          if (requireNamespace('dplyr', quietly = TRUE)) {
            ret <- dplyr::bind_rows(ret) # handle column name mismatches
          } else {
            ret <- do.call(rbind, ret)
          }
        } else {
          ret <- simplify2array(ret)

          if (force_df) {
            ret <- data.frame(result = ret)
          }
        }
      }

      if (!is.null(append_cols)) {
        ret <- cbind(sf::st_drop_geometry(y[, append_cols]), ret)
      }

      return(ret)
    }
  }, finally={
    raster::readStop(x)
    if (!is.null(weights)) {
      raster::readStop(weights)
    }
  })
}

.appendXY <- function(vals_df, rast, first_row, nrow, first_col, ncol) {
  if (nrow(vals_df) == 0) {
    vals_df$x <- numeric()
    vals_df$y <- numeric()
  } else {
    x_coords <- raster::xFromCol(rast, col=seq(first_col, first_col + ncol - 1))
    y_coords <- raster::yFromRow(rast, row=seq(first_row, first_row + nrow - 1))

    vals_df$x <- rep(x_coords, times=nrow)
    vals_df$y <- rep(y_coords, each=ncol)
  }

  return(vals_df)
}

.appendCell <- function(vals_df, rast, first_row, nrow, first_col, ncol) {
  if (nrow(vals_df) == 0) {
    vals_df$cell <- numeric()
  } else {
    rows <- rep(seq(first_row, first_row + nrow - 1), each = ncol)
    cols <- rep(seq(first_col, first_col + ncol - 1), times = nrow)

    vals_df$cell <- raster::cellFromRowCol(rast, row=rows, col=cols)
  }

  return(vals_df)
}

.cppStatColNames <- function(rast, stat_names, full_colnames, quantiles) {
  quantile_index = which(stat_names == 'quantile')
  if (length(quantile_index) != 0) {
    stat_names <- c(stat_names[seq_along(stat_names) < quantile_index],
                    sprintf('q%02d', as.integer(100 * quantiles)),
                    stat_names[seq_along(stat_names) > quantile_index])
  }

  if (raster::nlayers(rast) > 1 || full_colnames) {
    z <- expand.grid(names(rast), stat_names, stringsAsFactors=TRUE)
    mapply(paste, z[[2]], z[[1]], MoreArgs=list(sep='.'))
  } else {
    stat_names
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

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='sfc_GEOMETRY'), .exact_extract)
