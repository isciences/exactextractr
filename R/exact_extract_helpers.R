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
