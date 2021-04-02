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

#' Return column names to be used for summary operations
#'
#' @param value_names names of value raster layers
#' @param weight_names names of weighting raster layers
#' @param fun functions or names of summary operations
#' @param full_colnames return a complete column name even when there is no
#'                      ambiguity?
#' @param quantiles quantiles to use when \code{stat_names} contains \code{quantile}
#' @return character vector of column names
#' @keywords internal
.resultColNames <- function(value_names, weight_names, fun, full_colnames, quantiles=numeric()) {
  if (inherits(fun, 'standardGeneric')) {
    stat_names <- fun@generic[1]
  } else if (is.function(fun)) {
    stat_names <- 'fun'
  } else {
    stat_names <- fun
  }

  quantile_index = which(stat_names == 'quantile')
  if (length(quantile_index) != 0) {
    stat_names <- c(stat_names[seq_along(stat_names) < quantile_index],
                    sprintf('q%02d', as.integer(100 * quantiles)),
                    stat_names[seq_along(stat_names) > quantile_index])
  }

  ind <- .valueWeightIndexes(length(value_names), length(weight_names))
  vn <- value_names[ind$values]
  wn <- weight_names[ind$weights]

  if (length(vn) > 1 || full_colnames) {
    # determine all combinations of index and stat
    z <- expand.grid(index=seq_along(vn),
                     stat=stat_names, stringsAsFactors=FALSE)
    z$value <- vn[z$index]
    if (is.null(wn)) {
      z$weights <- NA
    } else {
      z$weights <- wn[z$index]
    }

    # construct column names for each index, stat
    # add weight layer name only if layer is ambiguously weighted
    mapply(function(stat, value, weight) {
      ret <- stat
      if (full_colnames || length(value_names) > 1) {
        ret <- paste(ret, value, sep='.')
      }
      if (.includeWeightInColName(stat) && ((full_colnames & length(weight_names) > 0)
                                            || length(weight_names) > 1)) {
        ret <- paste(ret, weight, sep='.')
      }

      return(ret)
    }, z$stat, z$value, z$weight, USE.NAMES = FALSE)
  } else {
    stat_names
  }
}

.includeWeightInColName <- function(fun) {
  .isWeighted(fun) | fun == 'fun'
}

.isWeighted <- function(stat_name) {
  stat_name %in% c('weighted_mean', 'weighted_sum')
}

#' Compute indexes for the value and weight layers that should be
#' processed together
#'
#' @param num_values number of layers in value raster
#' @param num_weights number of layers in weighting raster
#' @return list with \code{values} and \code{weights} elements
#'         providing layer indexes
#' @keywords internal
.valueWeightIndexes <- function(num_values, num_weights) {
  if (num_weights == 0) {
    vi <- seq_len(num_values)
    wi <- NA
  } else if (num_values == num_weights) {
    # process in parallel
    vi <- seq_len(num_values)
    wi <- seq_len(num_weights)
  } else if (num_values == 1 && num_weights > 1) {
    # recycle values
    vi <- rep.int(1, num_weights)
    wi <- seq_len(num_weights)
  } else if (num_values > 1 && num_weights == 1) {
    # recycle weights
    vi <- seq_len(num_values)
    wi <- rep.int(1, num_values)
  }

  list(values = vi, weights = wi)
}

.areaMethod <- function(crs_obj) {
  if (!(is.na(crs_obj)) && crs_obj$units_gdal == 'degree') {
    return('spherical')
  } else {
    return('cartesian')
  }
}

.validateFlag <- function(value, name) {
  if(!(is.logical(value) && length(value) == 1 && !is.na(value))) {
    stop(name, ' must be TRUE or FALSE')
  }
}

.validateNumericScalar <- function(value, name) {
  if (!(is.numeric(value) && length(value) == 1 && !is.na(value))) {
    stop(name, ' must be a single numeric value')
  }
}

.validateNumericScalarOrNA <- function(value, name) {
  if (!(is.numeric(value) && length(value) == 1)) {
    stop(name, ' must be a single numeric value')
  }
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

# Return the number of standard (non-...) arguments in a supplied function that
# do not have a default value. This is used to fail if the summary function
# provided by the user cannot accept arguments of values and weights.
.num_expected_args <- function(fun) {
  a <- formals(args(fun))
  a <- a[names(a) != '...']
  sum(sapply(a, nchar) == 0)
}

.startReading <- function(r) {
  if(inherits(r, 'BasicRaster')) {
    return(raster::readStart(r))
  } else if (inherits(r, 'SpatRaster')) {
    terra::readStart(r)
  }

  return(r)
}

.stopReading <- function(r) {
  if(inherits(r, 'BasicRaster')) {
    return(raster::readStop(r))
  } else if (inherits(r, 'SpatRaster')) {
    terra::readStop(r)
  }

  return(r)
}

.numLayers <- function(r) {
  if(inherits(r, 'BasicRaster')) {
    return(raster::nlayers(r))
  } else if (inherits(r, 'SpatRaster')) {
    return(terra::nlyr(r))
  } else {
    stop('Unknown type: ', class(r))
  }
}

.isRaster <- function(r) {
  inherits(r, 'BasicRaster') | inherits(r, 'SpatRaster')
}

.xFromCol <- function(r, col) {
  if (inherits(r, 'BasicRaster')) {
    raster::xFromCol(r, col)
  } else if (inherits(r, 'SpatRaster')) {
    terra::xFromCol(r, col)
  } else {
    stop('Unknown type: ', class(r))
  }
}

.yFromRow <- function(r, row) {
  if (inherits(r, 'BasicRaster')) {
    raster::yFromRow(r, row)
  } else if (inherits(r, 'SpatRaster')) {
    terra::yFromRow(r, row)
  } else {
    stop('Unknown type: ', class(r))
  }
}

.cellFromRowCol <- function(r, row, col) {
  if (inherits(r, 'BasicRaster')) {
    raster::cellFromRowCol(r, row, col)
  } else if (inherits(r, 'SpatRaster')) {
    terra::cellFromRowCol(r, row, col)
  } else {
    stop('Unknown type: ', class(r))
  }
}

.extent <- function(r) {
  if (inherits(r, 'BasicRaster')) {
    ex <- r@extent
    c(ex@xmin, ex@ymin, ex@xmax, ex@ymax)
  } else if (inherits(r, 'SpatRaster')) {
    ex <- terra::ext(r)
    c(ex$xmin, ex$ymin, ex$xmax, ex$ymax)
  } else {
    stop('Unknown type: ', class(r))
  }
}

.res <- function(r) {
  if (inherits(r, 'BasicRaster')) {
    raster::res(r)
  } else if (inherits(r, 'SpatRaster')) {
    terra::res(r)
  } else {
    stop('Unknown type: ', class(r))
  }
}

.getValuesBlock <- function(r, row, nrows, col, ncols) {
  if (inherits(r, 'BasicRaster')) {
    raster::getValuesBlock(r, row, nrows, col, ncols, format = 'm')
  } else if (inherits(r, 'SpatRaster')) {
    terra::readValues(r, row, nrows, col, ncols, mat = TRUE)
  } else {
    stop('Unknown type: ', class(r))
  }
}
