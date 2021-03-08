# Copyright (c) 2018-2021 ISciences, LLC.
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
#'                     and \code{weighted_sum} summary operations. When
#'                     \code{weights} is set to \code{'mean'}, the cell areas
#'                     of \code{x} will be calculated and used as weights.
#' @param     coverage_area  if \code{TRUE}, output pixel \code{coverage_area}
#'                           instead of \code{coverage_fraction}
#' @param     default_value  an optional value to use instead of \code{NA} in \code{x}
#' @param     default_weight an optional value to use instead of \code{NA} in
#'                           \code{weights}
#' @param     quantiles   quantiles to be computed when \code{fun == 'quantile'}
#' @param     append_cols when \code{fun} is not \code{NULL}, an optional
#'                        character vector of columns from \code{y} to be
#'                        included in returned data frame.
#' @param     force_df always return a data frame instead of a vector, even if
#'                     \code{x} has only one layer and \code{fun} has length 1
#' @param     summarize_df  pass values, coverage fraction/area, and weights to
#'                          \code{fun} as a single data frame instead of
#'                          separate arguments.
#' @param     full_colnames include the names of \code{x} in the names of the
#'                          returned data frame, even if \code{x} has only one
#'                          layer. This is useful when the results of multiple
#'                          calls to \code{exact_extract} are combined with
#'                          \code{cbind}.
#' @param     include_area if \code{TRUE}, and \code{fun} is \code{NULL}, augment
#'                       the returned data frame for each feature with a column
#'                       for the cell area. If the units of the raster CRS are
#'                       degrees, the area in square meters will be calculated
#'                       based on a spherical approximation of Earth. Otherwise,
#'                       a Cartesian area will be calculated (and will be the
#'                       same for all pixels.) If \code{TRUE} and \code{fun} is
#'                       not \code{NULL}, add \code{area} to the data frame passed
#'                       to \code{fun} for each feature.
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
                           summarize_df=FALSE,
                           full_colnames=FALSE,
                           stack_apply=FALSE,
                           append_cols=NULL,
                           include_area=FALSE,
                           include_cols=NULL,
                           coverage_area=FALSE,
                           quantiles=NULL,
                           default_value=NA_real_,
                           default_weight=NA_real_) {
  area_weights <- is.character(weights) && length(weights) == 1 && weights == 'area'
  if (area_weights) {
    weights <- NULL
  }

  .validateFlag(coverage_area, 'coverage_area')
  .validateFlag(force_df, 'force_df')
  .validateFlag(full_colnames, 'full_colnames')
  .validateFlag(include_area, 'include_area')
  .validateFlag(include_cell, 'include_cell')
  .validateFlag(include_xy, 'include_xy')
  .validateFlag(progress, 'progress')
  .validateFlag(stack_apply, 'stack_apply')
  .validateFlag(summarize_df, 'summarize_df')
  .validateNumericScalar(max_cells_in_memory, 'max_cells_in_memory')
  .validateNumericScalarOrNA(default_value, 'default_value')
  .validateNumericScalarOrNA(default_weight, 'default_weight')

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
    if (!inherits(weights, 'BasicRaster')) {
      stop("Weights must be a Raster object or \"area\".")
    }

    if (is.character(fun) && !any(startsWith(fun, "weighted"))) {
      warning("Weights provided but no requested operations use them.")
    }

    if (inherits(weights, 'BasicRaster') && !is.na(sf::st_crs(x))) {
      if (is.na(sf::st_crs(weights))) {
        warning("No CRS specified for weighting raster; assuming it has the same CRS as the value raster.")
      } else if (sf::st_crs(x) != sf::st_crs(weights)) {
        stop("Weighting raster does not have the same CRS as value raster.")
      }
    }
  }

  analysis_crs <- sf::st_crs(x)
  if(is.na(sf::st_crs(x)) && !is.na(sf::st_crs(y))) {
    warning("No CRS specified for raster; assuming it has the same CRS as the polygons.")
    analysis_crs <- sf::st_crs(y)
  } else if(is.na(sf::st_crs(y)) && !is.na(sf::st_crs(x))) {
    warning("No CRS specified for polygons; assuming they have the same CRS as the raster.")
  } else if(sf::st_crs(x) != sf::st_crs(y)) {
    y <- sf::st_transform(y, sf::st_crs(x))
    warning("Polygons transformed to raster CRS (EPSG:", sf::st_crs(x)$epsg, ")")
  }
  if(area_weights || include_area || coverage_area) {
    area_method <- .areaMethod(analysis_crs)
  } else {
    area_method <- NULL
  }

  if (is.character(fun)) {
    if (length(fun) == 0) {
      stop("No summary operations provided.")
    }

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

    if (include_area) {
      stop("include_area must be FALSE for named summary operations")
    }

    if (!is.null(include_cols)) {
      stop("include_cols not supported for named_summary operations (see argument append_cols)")
    }

    if (summarize_df) {
      stop("summarize_df can only be used when `fun` is an R function")
    }
  } else if (is.function(fun)) {
    if (summarize_df) {
      if (.num_expected_args(fun) < 1) {
        stop("exact_extract was called with a function that does not appear to ",
             "be of the form `function(df, ...`).")
      }
    } else if (is.null(weights)) {
      if (.num_expected_args(fun) < 2) {
        stop("exact_extract was called with a function that does not appear to ",
             "be of the form `function(values, coverage_fractions, ...)`. If ",
             "the summary function should accept a single data frame argument, ",
             "set `summarize_df = TRUE`.")
      }
    } else if (.num_expected_args(fun) < 3) {
        stop("exact_extract was called with a function that does not appear to ",
             "be of the form `function(values, coverage_fractions, weights, ...)`.",
             "If the summary function should accept a single data frame argument, ",
             "set `summarize_df = TRUE`.")
    }
  } else if (is.null(fun)) {
    if (length(list(...)) > 0) {
      stop("Unexpected arguments: ", paste(names(list(...)), collapse = ','))
    }

    if (summarize_df) {
      stop("summarize_df can only be used when `fun` is an R function")
    }
  } else {
    stop("fun must be a character vector, function, or NULL")
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
    if (inherits(weights, 'BasicRaster')) {
      weights <- raster::readStart(weights)
    }

    if (is.character(fun)) {
      # Compute all stats in C++.
      # CPP_stats returns a matrix, which gets turned into a column by sapply
      # Results has one column per feature and one row per stat/raster layer
      if((is.null(weights) && !area_weights) && any(.isWeighted(fun))) {
        stop("Weighted stat requested but no weights provided.")
      }

      results <- sapply(sf::st_as_binary(geoms, EWKB=TRUE), function(wkb) {
        ret <- CPP_stats(x, weights, wkb, default_value, default_weight, coverage_area, area_method, fun, max_cells_in_memory, quantiles)
        update_progress()
        return(ret)
      })

      if (!is.matrix(results) && !force_df) {
        # Single stat? Return a vector unless asked otherwise via force_df.
        return(results)
      } else {
        # Return a data frame with a column for each stat
        colnames <- .resultColNames(names(x), names(weights), fun, full_colnames, quantiles)

        if (is.matrix(results)) {
          results <- t(results)
        } else {
          results <- matrix(results, nrow=length(results))
        }

        dimnames(results) <- list(NULL, colnames)
        ret <- as.data.frame(results)

        # drop duplicated columns (occurs when an unweighted stat is
        # requested alongside a weighted stat, with a stack of weights)
        ret <- ret[, unique(names(ret)), drop = FALSE]

        if (!is.null(append_cols)) {
          ret <- cbind(sf::st_drop_geometry(y[, append_cols]), ret)
        }

        return(ret)
      }
    } else {
      num_values <- raster::nlayers(x)
      value_names <- names(x)
      if (inherits(weights, 'BasicRaster')) {
        num_weights <- raster::nlayers(weights)
        weight_names <- names(weights)
      } else if (area_weights) {
        num_weights <- 1
        weight_names <- 'area'
        weights <- NULL
      } else {
        num_weights <- 0
        weight_names <- NULL
      }

      if (stack_apply || (num_values == 1 && num_weights <= 1)) {
        apply_layerwise <- TRUE

        if (num_values > 1 && num_weights > 1 && num_values != num_weights) {
          stop(sprintf("Can't apply function layerwise with stacks of %d value layers and %d layers", num_values, num_weights))
        }

        result_names <- .resultColNames(value_names, weight_names, fun, full_colnames)
        num_results <- max(num_weights, num_values)
        ind <- .valueWeightIndexes(num_values, num_weights)
      } else {
        apply_layerwise <- FALSE
      }

      # For R summary functions and data frame output, we avoid using
      # input layer names if there is only one value/weight layer
      if (length(value_names) == 1) {
        value_names <- 'value'
      }
      if (length(weight_names) == 1) {
        weight_names <- 'weight'
      }

      ret <- lapply(seq_along(geoms), function(feature_num) {
        wkb <- sf::st_as_binary(geoms[[feature_num]], EWKB=TRUE)

        if (!is.null(include_cols)) {
          include_col_values <- sf::st_drop_geometry(y[feature_num, include_cols])
        } else {
          include_col_values <- NULL
        }

        # only raise a disaggregation warning for the first feature
        warn_on_disaggregate <- feature_num == 1

        col_list <- CPP_exact_extract(x,
                                      weights,
                                      wkb,
                                      default_value,
                                      default_weight,
                                      include_xy,
                                      include_cell,
                                      include_area,
                                      area_weights,
                                      coverage_area,
                                      area_method,
                                      include_col_values,
                                      value_names,
                                      weight_names,
                                      warn_on_disaggregate)
        if (coverage_area) {
          coverage_col <- 'coverage_area'
        } else {
          coverage_col <- 'coverage_fraction'
        }

        if (!is.null(include_cols)) {
          # Replicate the include_cols vectors to be as long as the other columns,
          # so we can use quickDf
          nrow <- length(col_list[[coverage_col]])
          col_list[include_cols] <- lapply(col_list[include_cols], rep, nrow)
        }
        df <- .quickDf(col_list)

        update_progress()

        # No summary function? Return the whole data frame.
        if (is.null(fun)) {
          return(df)
        }

        # Summary function accepts a data frame? Pass it along.
        if (summarize_df && !apply_layerwise) {
          return(fun(df, ...))
        }

        # Break the data frame into components that we can pass
        # separately to the summary functions.
        included_cols_df <- df[, !(names(df) %in% c(value_names, weight_names, coverage_col)), drop = FALSE]
        vals_df <- df[, value_names, drop = FALSE]
        weights_df <- df[, weight_names, drop = FALSE]
        cov_fracs <- df[[coverage_col]]

        if (apply_layerwise) {
          result <- lapply(seq_len(num_results), function(i) {
            if (summarize_df) {
              # Pack everything into a single data frame and pass it to `fun`
              arg_df <- cbind(value = vals_df[[ind$values[i]]],
                              included_cols_df)
              if (num_weights > 0) {
                arg_df$weight <- weights_df[[ind$weights[i]]]
              }
              arg_df[[coverage_col]] <- df[[coverage_col]]

              return(fun(arg_df, ...))
            } else {
              # Pull values and weights out into vectors, unless we have
              # included columns (x/y/cell, etc.) in which case values
              # remain a data frame. Retained for backwards compat.

              vx <- vals_df[, ind$values[i]]
              if (ncol(included_cols_df) > 0) {
                vx <- cbind(data.frame(value = vx), included_cols_df)
              }
              if (num_weights == 0) {
                return(fun(vx, cov_fracs, ...))
              } else {
                return(fun(vx, cov_fracs, weights_df[, ind$weights[i]], ...))
              }
            }
          })

          if (num_results == 1) {
            return(result[[1]])
          }

          names(result) <- result_names

          return(.quickDf(result))
        } else {
          # Pass all layers to callback, to be handled together
          # Included columns (x/y/cell) are passed with the values.
          # Pass single-column data frames as vectors.
          vals_df <- .singleColumnToVector(cbind(vals_df, included_cols_df))
          weights_df <- .singleColumnToVector(weights_df)

          if (num_weights == 0) {
            return(fun(vals_df, cov_fracs, ...))
          } else {
            return(fun(vals_df, cov_fracs, weights_df, ...))
          }
        }
      })

      if (!is.null(fun)) {
        if (all(sapply(ret, is.data.frame))) {
        # function returned a data frame for each polygon? rbind them
          if (requireNamespace('dplyr', quietly = TRUE)) {
            ret <- dplyr::bind_rows(ret) # handle column name mismatches
          } else {
            ret <- do.call(rbind, ret)
          }
        } else {
          # function returned something else; combine the somethings into
          # an array
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
    if (inherits(weights, 'BasicRaster')) {
      raster::readStop(weights)
    }
  })
}

# faster replacement for as.data.frame when input is a named list
# with equal-length columns
# from Advanced R, sec. 24.4.2
.quickDf <- function(lst) {
  class(lst) <- 'data.frame'
  attr(lst, 'row.names') <- .set_row_names(length(lst[[1]]))
  lst
}

.singleColumnToVector <- function(df) {
  if (ncol(df) == 1) {
    df[, 1]
  } else {
    df
  }
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

#' @import sf
#' @import raster
#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='sf'),
          .exact_extract)

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='SpatialPolygonsDataFrame'),
          function(x, y, ...) .exact_extract(x, sf::st_as_sf(y), ...))

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='SpatialPolygons'),
          function(x, y, ...) .exact_extract(x, sf::st_as_sf(y), ...))

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
