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
#' Extracts the values of cells in `Raster*` that are covered by polygons in a
#' simple feature collection (`sf` or `sfc`) or `SpatialPolygonsDataFrame`.
#' Returns either a summary of the extracted values or the extracted values
#' themselves.
#'
#' @details
#' `exact_extract` extracts the values of cells in a `Raster*` that are covered
#' by polygons in a simple feature collection (`sf` or `sfc`) or
#' `SpatialPolygonDataFrame`, as well as the fraction or area of each cell that
#' is covered by the polygon. The function can either return these values
#' directly to the caller, or can return the result of a predefined summary
#' operation or user-defined R function applied to the values. These three
#' approaches are described in the subsections below.
#'
#' ## Returning extracting values directly
#'
#' If `fun` is not specified, `exact_extract` will return a list with
#' one data frame for each feature in the input feature collection. The data
#' frame will contain a column with cell values from each layer in the input
#' `Raster*` (and optional weighting `Raster*`) and a column indicating
#' the fraction or area of the cell that is covered by the polygon.
#'
#' If the input rasters have only one layer, the corresponding columns in the
#' data frame will be named `values` or `weights`. When the input rasters have
#' more than one layer, the columns will be named according to `names(x)` and
#' `names(weights)`. The column containing pixel coverage will be called
#' `coverage_fraction` when `coverage_area = FALSE`, or `coverage_area` when
#' `coverage_area = TRUE`.
#'
#' If the output data frames are to be combined (e.g., with `rbind`, it may be
#' useful to include identifying column(s) from the input features in the
#' returned data frames using `include_cols`. Additional columns can be added
#' to the returned data frames with the `include_area`, `include_cell`, and
#' `include_xy` arguments.
#'
#' ## Predefined summary operations
#'
#' Often the individual pixel values are not needed; only one or more summary
#' statistics (e.g., mean, sum) is required for each polygon. Common summary
#' statistics can be calculated by `exact_extract` directly using a predefined
#' summary operation. Where possible, this approach is advantageous because it
#' allows the package to calculate the statistics incrementally, avoiding the
#' need to store all pixel values in memory at the same time. This allows the
#' package to process arbitrarily large data with a small amount of memory. (The
#' `max_pixels_in_memory` argument can be used to fine-tune the amount of memory
#' made available to `exact_extract`.)
#'
#' To summarize pixel values using a predefined summary option, `fun` should be
#' set to a character vector of one or more operation names. If the input raster
#' has a single layer and a single summary operation is specified,
#' `exact_extract` will return a vector with the result of the summary operation
#' for each feature in the input. If the input raster has multiple layers, or if
#' multiple summary operations are specified, `exact_extract` will return a data
#' frame with a row for each feature and a column for each summary operation /
#' layer combination. (The `force_df` option can be used to always return a data
#' frame instead of a vector.)
#'
#' The following summary operations are supported:
#'
#'  * `min` - the minimum defined (non-`NA`) value in any raster cell wholly or
#'            partially covered by the polygon
#'  * `max` - the maximum defined (non-`NA`) value in any raster cell wholly or
#'            partially covered by the polygon
#'  * `count` - the sum of fractions of raster cells with defined non-`NA`
#'              values covered by the polygon
#'  * `sum`   - the sum of defined (non-`NA`) raster cell values, multiplied by
#'              the fraction of the cell that is covered by the polygon
#'  * `mean` - the mean cell value, weighted by the fraction of each cell
#'             that is covered by the polygon
#'  * `median` - the median cell value, weighted by the fraction of each cell
#'               that is covered by the polygon
#'  * `quantile` - arbitrary quantile(s) of cell values, specified in
#'                 `quantiles`, weighted by the fraction of each cell that is
#'                  covered by the polygon
#'  * `mode` - the most common cell value, weighted by the fraction of
#'             each cell that is covered by the polygon. Where multiple
#'             values occupy the same maximum number of weighted cells,
#'             the largest value will be returned.
#'  * `majority` - synonym for `mode`
#'  * `minority` - the least common cell value, weighted by the fraction
#'                 of each cell that is covered by the polygon. Where
#'                 multiple values occupy the same minimum number of
#'                 weighted cells, the smallest value will be returned.
#'  * `variety` - the number of distinct values in cells that are wholly or
#'                partially covered by the polygon.
#'  * `variance` - the population variance of cell values, weighted by the
#'                 fraction of each cell that is covered by the polygon.
#'  * `stdev` - the population standard deviation of cell values, weighted by
#'              the fraction of each cell that is covered by the polygon.
#'  * `coefficient_of_variation` - the population coefficient of variation of
#'                                 cell values, weighted by the fraction of each
#'                                 cell that is covered by the polygon.
#'  * `weighted_mean` - the mean cell value, weighted by the product of
#'                      the fraction of each cell covered by the polygon
#'                      and the value of a second weighting raster provided
#'                      as `weights`
#'  * `weighted_sum` - the sum of defined raster cell values, multiplied by
#'                     the fraction of each cell that is covered by the polygon
#'                     and the value of a second weighting raster provided
#'                     as `weights`
#'
#' In all of the summary operations, `NA` values in the the primary raster (`x`)
#' raster are ignored (i.e., `na.rm = TRUE`.) If `NA` values occur in the
#' weighting raster, the result of the weighted operation will be `NA`. `NA`
#' values in both `x` and `weights` can be replaced on-the-fly using the
#' `default_value` and `default_weight` arguments.
#'
#' ## User-defined summary functions
#'
#' If no predefined summary operation is suitable, a user-defined R function may
#' be provided as `fun`. The function will be called once for each feature and
#' must return either a single value or a data frame. The results of the
#' function for each feature will be combined and returned by `exact_extract`.
#'
#' The simplest way to write a summary function is to set
#' argument `summarize_df = TRUE`. (For backwards compatibility, this is not the
#' default.) In this mode, the summary function takes the signature
#' `function(df, ...)` where `df` is the same data frame that would be returned
#' by `exact_extract` with `fun = NULL`.
#'
#' With `summarize_df = FALSE`, the function must have the signature
#' `function(values, coverage_fractions, ...)` when weights are not used, and
#' `function(values, coverage_fractions, weights, ...)` when weights are used.
#' If the value and weight rasters are `RasterLayers`, the function arguments
#' will be vectors; if either is a `RasterStack`, the function arguments will
#' be data frames, with column names taken from the names of the value/weight
#' rasters. Values brought in through the `include_xy`, `include_area`,
#' `include_cell`, and `include_cols` arguments will be added to the `values`
#' data frame. For most applications, it is simpler to set `summarize_df = TRUE`
#' and work with all inputs in a single data frame.
#'
#' @param     x a `RasterLayer`, `RasterStack`, or `RasterBrick`
#' @param     y a `sf`, `sfc`, `SpatialPolygonsDataFrame`, or `SpatialPolygons`
#'            object with polygonal geometries
#' @param     fun an optional function or character vector, as described below
#' @param     weights  a weighting raster to be used with the `weighted_mean`
#'                     and `weighted_sum` summary operations, or a user-defined
#'                     summary function. When `weights` is set to `'area'`, the
#'                     cell areas of `x` will be calculated and used as weights.
#' @param     coverage_area  if `TRUE`, output pixel `coverage_area`
#'                           instead of `coverage_fraction`
#' @param     default_value  an optional value to use instead of `NA` in `x`
#' @param     default_weight an optional value to use instead of `NA` in `weights`
#' @param     quantiles   quantiles to be computed when `fun = 'quantile'`
#' @param     append_cols when `fun` is not `NULL`, an optional character vector
#'                        of columns from `y` to be included in returned data frame.
#' @param     force_df always return a data frame instead of a vector, even if
#'                     `x` has only one layer and `fun` has length 1
#' @param     summarize_df  pass values, coverage fraction/area, and weights to
#'                          `fun` as a single data frame instead of
#'                          separate arguments.
#' @param     full_colnames include the names of `x` and `weights` in
#'                          the names of the data frame for each feature, even if
#'                          `x` or `weights` has only one layer.
#'                          This is useful when the results of multiple
#'                          calls to `exact_extract` are combined with
#'                          `cbind`.
#' @param     include_area if `TRUE`, and `fun` is `NULL`, augment
#'                       the data frame for each feature with a column
#'                       for the cell area. If the units of the raster CRS are
#'                       degrees, the area in square meters will be calculated
#'                       based on a spherical approximation of Earth. Otherwise,
#'                       a Cartesian area will be calculated (and will be the
#'                       same for all pixels.) If `TRUE` and `fun` is
#'                       not `NULL`, add `area` to the data frame passed
#'                       to `fun` for each feature.
#' @param     include_cell if `TRUE`, and `fun` is `NULL`, augment
#'                       the data frame for each feature with a column
#'                       for the cell index (`cell`). If `TRUE` and
#'                       `fun` is not `NULL`, add `cell` to the
#'                       data frame passed to `fun` for each feature.
#' @param     include_cols an optional character vector of column names in
#'                         `y` to be added to the data frame for each
#'                         feature that is either returned (when `fun` is
#'                         `NULL`) or passed to `fun`.
#' @param     include_xy if `TRUE`, and `fun` is `NULL`, augment
#'                       the returned data frame for each feature with columns
#'                       for cell center coordinates (`x` and `y`). If
#'                       `TRUE` and `fun` is not `NULL`, add
#'                       `x` and `y` to the data frame passed to `fun`
#'                       for each feature.
#' @param     stack_apply   if `TRUE`, apply `fun` independently to
#'                          each layer or `x` (and its corresponding layer
#'                          of `weights`, if provided.) The number of
#'                          layers in `x` and `weights` must equal
#'                          each other or `1`, in which case the
#'                          single layer raster will be recycled.
#'                          If `FALSE`, apply `fun` to all layers of
#'                          `x` (and `weights`) simultaneously.
#' @param     max_cells_in_memory the maximum number of raster cells to load at
#'                                a given time when using a named summary operation
#'                                for `fun` (as opposed to a function defined using
#'                                R code). If a polygon covers more than `max_cells_in_memory`
#'                                raster cells, it will be processed in multiple chunks.
#' @param     progress if `TRUE`, display a progress bar during processing
#' @param     ... additional arguments to pass to `fun`
#' @return a vector, data frame, or list of data frames, depending on the type
#'         of `x` and the value of `fun` (see Details)
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
#' @md
NULL

.exact_extract <- function(x, y, fun=NULL, ...,
                           weights=NULL,
                           append_cols=NULL,
                           coverage_area=FALSE,
                           default_value=NA_real_,
                           default_weight=NA_real_,
                           include_area=FALSE,
                           include_cell=FALSE,
                           include_cols=NULL,
                           include_xy=FALSE,
                           force_df=FALSE,
                           full_colnames=FALSE,
                           stack_apply=FALSE,
                           summarize_df=FALSE,
                           quantiles=NULL,
                           progress=TRUE,
                           max_cells_in_memory=30000000
                           ) {
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

  if(!is.null(include_cols)) {
    if (!inherits(y, 'sf')) {
      stop(sprintf('include_cols only supported for sf arguments (received %s)',
                   paste(class(y), collapse = ' ')))
    }
  }

  if(sf::st_geometry_type(y, by_geometry = FALSE) == 'GEOMETRY') {
    if (!all(sf::st_dimension(y) == 2)) {
      stop("Features in sfc_GEOMETRY must be polygonal")
    }
  }

  if(!is.null(weights)) {
    if (!(.isRaster(weights) || weights == 'area')) {
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
    if (stack_apply) {
      stop("stack_apply can only be used when `fun` is a summary operation or function")
    }
    if (!is.null(append_cols)) {
      stop("append_cols can only be used when `fun` is a summary operation or function. See `include_cols`.")
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
    x <- .startReading(x)
    weights <- .startReading(weights)

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
      num_values <- .numLayers(x)
      value_names <- names(x)
      if (.isRaster(weights)) {
        num_weights <- .numLayers(weights)
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

      # if we have columns to append, iterate over the results and
      # append the columns to each, coercing into a data frame if
      # necessary. We need to iterate over the results before binding
      # them together because we don't know the legnth of each result
      # or even if their lengths are constant.
      if (!is.null(append_cols)) {
        for (i in seq_along(ret)) {
          if (!is.data.frame(ret[[i]])) {
            ret[[i]] = data.frame(result = ret[[i]])
          }

          if (nrow(ret[[i]]) >= 1) {
            append_row <- i
          } else {
            # cbinding row 0 cleanly brings in zero-length columns
            # of the correct names and types, in the correct positions
            append_row <- 0
          }

          ret[[i]] <- cbind(sf::st_drop_geometry(y[append_row, append_cols]),
                            ret[[i]],
                            row.names = NULL)
        }
      }

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

      return(ret)
    }
  }, finally={
    .stopReading(x)
    .stopReading(weights)
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

# CRAN version of sf does not explicitly declare this class. If we do not do it
# ourselves, documentation is not generated for the `sfc_GEOMETRYCOLLECTION`
# overload.
setOldClass('sfc_GEOMETRYCOLLECTION')

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='Raster', y='sfc_GEOMETRYCOLLECTION'), .exact_extract)

#' @import sf
#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='SpatRaster', y='sf'),
          .exact_extract)

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='SpatRaster', y='SpatialPolygonsDataFrame'),
          function(x, y, ...) .exact_extract(x, sf::st_as_sf(y), ...))

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='SpatRaster', y='SpatialPolygons'),
          function(x, y, ...) .exact_extract(x, sf::st_as_sf(y), ...))

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='SpatRaster', y='sfc_MULTIPOLYGON'), .exact_extract)

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='SpatRaster', y='sfc_POLYGON'), .exact_extract)

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='SpatRaster', y='sfc_GEOMETRY'), .exact_extract)

#' @useDynLib exactextractr
#' @rdname exact_extract
#' @export
setMethod('exact_extract', signature(x='SpatRaster', y='sfc_GEOMETRYCOLLECTION'), .exact_extract)
