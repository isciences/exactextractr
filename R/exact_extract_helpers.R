# Copyright (c) 2018-2022 ISciences, LLC.
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

.resultColNames <- function(...) {
  .resultColumns(...)$colname
}

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
.resultColumns <- function(value_names, weight_names, fun, full_colnames, quantiles=numeric(), unique_values=numeric(), colname_fun = NULL) {
  if (inherits(fun, 'standardGeneric')) {
    stat_names <- fun@generic[1]
  } else if (is.function(fun)) {
    stat_names <- 'fun'
  } else if (is.null(fun)) {
    stat_names <- ''
  } else {
    stat_names <- fun
  }

  quantile_index = which(stat_names == 'quantile')
  if (length(quantile_index) != 0) {
    stat_names <- splice(stat_names, quantile_index, rep('quantile', length(quantiles)))
  }

  for (stat in c('frac', 'weighted_frac')) {
    frac_index = which(stat_names == stat)
    if (length(frac_index) != 0) {
      stat_names <- splice(stat_names, frac_index, rep(stat, length(unique_values)))
    }
  }

  ind <- .valueWeightIndexes(length(value_names), length(weight_names))
  vn <- value_names[ind$values]
  wn <- weight_names[ind$weights]

  # determine all combinations of index and stat
  z <- expand.grid(index=seq_along(vn),
                   stat_name=stat_names, stringsAsFactors=FALSE)

  z$values <- vn[z$index]
  z$base_value <- NA_real_

  for (stat in c('frac', 'weighted_frac')) {
    ifrac <- which(z$stat_name == stat)
    z$base_value[ifrac] <- rep(unique_values, each = length(ifrac) / length(unique_values))
  }

  iquantile <- which(z$stat_name == 'quantile')
  z$base_value[iquantile] <- rep(quantiles, each = length(iquantile) / length(quantiles))

  if (is.null(wn)) {
    z$weights <- NA
  } else {
    z$weights <- wn[z$index]
  }
  z$weights[!.includeWeightInColName(z$stat_name)] <- NA

  if (is.null(colname_fun)) {
    colname_fun <- function(...) {
      .makeColname(full_colnames = full_colnames, ...)
    }
  }

  z$colname <- mapply(colname_fun,
                      fun_name = z$stat_name,
                      values = z$values,
                      weights = z$weights,
                      fun_value = z$base_value,
                      MoreArgs = list(nvalues = length(value_names),
                                      nweights = length(weight_names)),
                      USE.NAMES = FALSE)

  return(z)
}

.makeColname <- function(fun_name, values, weights, fun_value, full_colnames, nvalues, nweights) {
  # construct column names for each index, stat
  # add weight layer name only if layer is ambiguously weighted
  if (fun_name == 'quantile') {
    fun_component <- sprintf('q%02d', as.integer(100 * fun_value))
  } else if (fun_name %in% c('frac', 'weighted_frac')) {
    fun_component <- sprintf('%s_%s', fun_name, fun_value)
  } else {
    fun_component <- fun_name
  }

  ret <- fun_component
  if (full_colnames || nvalues > 1) {
    ret <- paste(ret, values, sep='.')
  }
  if ((!is.na(weights)) && ((full_colnames & nweights > 0) || nweights > 1)) {
    ret <- paste(ret, weights, sep='.')
  }

  return(ret)
}

.includeWeightInColName <- function(fun) {
  .isWeighted(fun) | fun == 'fun'
}

.isWeighted <- function(stat_name) {
  stat_name %in% c('weighted_mean', 'weighted_sum', 'weighted_frac')
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

.validateUniqueNames <- function(x) {
  nm <- names(x)
  if (!is.null(nm)) {
    if (length(nm) != length(unique(nm))) {
      stop('names of input rasters must be unique')
    }
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

.crs <- function(r) {
  if(inherits(r, 'BasicRaster')) {
    if (utils::packageVersion('raster') < numeric_version('3.5')) {
      return(raster::crs(r))
    } else {
      return(terra::crs(r))
    }
  } else if (inherits(r, 'SpatRaster')) {
    return(terra::crs(r))
  } else {
    stop('Unknown type: ', class(r))
  }
}

.setValues <- function(r, x) {
  if(inherits(r, 'BasicRaster')) {
    if (utils::packageVersion('raster') < numeric_version('3.5')) {
      raster::values(r) <- x
    } else {
      terra::values(r) <- x
    }
  } else if (inherits(r, 'SpatRaster')) {
    raster::values(r) <- x
  } else {
    stop('Unknown type: ', class(r))
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

.colFromX <- function(r, x) {
  if (inherits(r, 'BasicRaster')) {
    raster::colFromX(r, x)
  } else if (inherits(r, 'SpatRaster')) {
    terra::colFromX(r, x)
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

.rowFromY <- function(r, y) {
  if (inherits(r, 'BasicRaster')) {
    raster::rowFromY(r, y)
  } else if (inherits(r, 'SpatRaster')) {
    terra::rowFromY(r, y)
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

.createProgress <- function(progress, n) {
  if (progress && n > 1) {
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

  return(update_progress)
}

.isInMemory <- function(r) {
  if (inherits(r, 'BasicRaster')) {
    return(raster::inMemory(r))
  } else if (inherits(r, 'SpatRaster')) {
    return(terra::inMemory(r)[1])
  } else {
    stop('Unknown type: ', class(r))
  }
}

.netCDFBlockSize <- function(fname, varname) {
  nc <- NULL
  sz <- NA

  tryCatch({
    nc <- ncdf4::nc_open(fname)
    sz <- nc$var[[varname]]$chunksizes

    dim_index <- nc$var[[varname]]$dimids + 1L
    dim_names <- sapply(dim_index, function(i) nc$dim[[i]]$name)

    if (all(is.na(sz))) {
      # file is not compressed
      sz <- rep.int(1, length(dim_names))
    }

    names(sz) <- dim_names
  }, finally = {
    if (!is.null(nc)) {
      ncdf4::nc_close(nc)
    }
  })

  # flip dimensions 1 and 2 so we return row/col
  return(sz[c(2, 1, seq_along(sz)[-(1:2)])])
}

.blockSize <- function(r) {
  # set default return value in case file is uncompressed and has
  # has no block size, or we simply can't figure it out
  ret <- c(1, 1)

  if (inherits(r, 'BasicRaster')) {
    if (r[[1]]@file@driver == 'netcdf') {
      ret <- .netCDFBlockSize(r[[1]]@file@name, attr(r[[1]]@data, 'zvar'))
    } else if (r@file@driver == 'gdal') {
      ret <- c(r@file@blockrows, r@file@blockcols)
    }
  } else if (inherits(r, 'SpatRaster')) {
    ret <- terra::fileBlocksize(r)[1, ]
  }

  unname(ret)
}

.eagerLoad <- function(r, geoms, max_cells_in_memory, message_on_fail) {
  if (is.null(r)) {
    return(NULL)
  }

  cells_required <- .numCells(r, geoms)

  if (cells_required <= max_cells_in_memory) {
    box <- sf::st_bbox(geoms)
    geom_ext <- terra::ext(box[c('xmin', 'xmax', 'ymin', 'ymax')])

    if (!inherits(r, 'SpatRaster')) {
      # current CRAN version of terra (1.4-22) does not preserve
      # names on conversion (https://github.com/rspatial/terra/issues/430)
      nm <- names(r)
      r <- terra::rast(r)
      names(r) <- nm
    }

    overlap_ext <- terra::intersect(terra::ext(r), geom_ext)
    if (is.null(overlap_ext)) {
      # Extents do not overlap, and terra::crop will throw an error
      # if we try to crop. Return the input raster as-is; nothing will be
      # read from it anyway.
      return(r)
    }

    r <- terra::crop(r, geom_ext, snap = 'out')
  } else if (message_on_fail) {
    message('Cannot preload entire working area of ', cells_required,
            ' cells with max_cells_in_memory = ', max_cells_in_memory, '.',
            ' Raster values will be read for each feature individually.')
  }

  return(r)
}

.numCells <- function(r, g) {
  if (is.null(r)) {
    return(0)
  }

  box <- sf::st_bbox(g)

  top <- .rowFromY(r, box['ymax'])
  bottom <- .rowFromY(r, box['ymin'])
  left <- .colFromX(r, box['xmin'])
  right <- .colFromX(r, box['xmax'])

  if (is.na(top) && is.na(bottom)) {
    return(0L)
  }
  if (is.na(left) && is.na(right)) {
    return(0L)
  }
  if (is.na(top)) {
    top <- 1
  }
  if (is.na(bottom)) {
    bottom <- nrow(r)
  }
  if (is.na(left)) {
    left <- 1
  }
  if (is.na(right)) {
    right <- ncol(r)
  }

  return( (bottom - top + 1) * (right - left + 1) * .numLayers(r) )
}

splice <- function(x, i, replacement) {
  c(x[seq_along(x) < i],
    replacement,
    x[seq_along(x) > i])
}
