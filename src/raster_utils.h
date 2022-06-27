// Copyright (c) 2018-2022 ISciences, LLC.
// All rights reserved.
//
// This software is licensed under the Apache License, Version 2.0 (the "License").
// You may not use this file except in compliance with the License. You may
// obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#pragma once

// [[Rcpp::plugins("cpp14")]]
#include <Rcpp.h>

#include <memory.h>

#include "exactextract/src/grid.h"
#include "exactextract/src/raster.h"
#include "exactextract/src/raster_area.h"

// Construct a grid corresponding to 'rast'
exactextract::Grid<exactextract::bounded_extent> make_grid(const Rcpp::S4 & rast);
exactextract::Grid<exactextract::bounded_extent> make_grid(const Rcpp::NumericVector & extent, const Rcpp::NumericVector & res);

// Return the number of layers in 'rast'
int get_nlayers(Rcpp::S4 & rast);

// Return a column number in 'rast' for each cell in 'grid'
Rcpp::IntegerVector cols_for_x(Rcpp::S4 & rast, exactextract::Grid<exactextract::bounded_extent> grid);

// Return a row number in 'rast' for each cell in 'grid'
Rcpp::IntegerVector rows_for_y(Rcpp::S4 & rast, exactextract::Grid<exactextract::bounded_extent> grid);

// Return a vector or x values in 'rast' for each cell in 'grid'
Rcpp::NumericVector get_x_values(Rcpp::S4 & rast, exactextract::Grid<exactextract::bounded_extent> grid);

// Return a vector of y values in 'rast' for each cell in 'grid'
Rcpp::NumericVector get_y_values(Rcpp::S4 & rast, exactextract::Grid<exactextract::bounded_extent> grid);

// Return a vector of cell numbers in 'rast' for each cell in 'grid'
Rcpp::NumericVector get_cell_numbers(Rcpp::S4 & rast, exactextract::Grid<exactextract::bounded_extent> grid);

// Construct a row-major vector of the values in 'r'
template<typename T>
Rcpp::NumericVector as_vector(const exactextract::AbstractRaster<T> & r) {
  // Convert Raster to a vector, using row-major storage (consistent with
  // the Raster package and distinct from the R representation of matrices)
  Rcpp::NumericVector ret = Rcpp::no_init(r.rows() * r.cols());

  size_t k = 0;
  for (size_t i = 0; i < r.rows(); i++) {
    for (size_t j = 0; j < r.cols(); j++) {
      ret[k++] = r(i, j);
    }
  }

  return ret;
}

template<typename T>
bool requires_stored_values(T s) {
  return s == "mode" ||
    s == "majority" ||
    s == "minority" ||
    s == "variety" ||
    s == "median" ||
    s == "quantile" ||
    s == "frac" ||
    s == "weighted_frac";
}

template<typename T, typename G>
std::unique_ptr<exactextract::AbstractRaster<double>> get_area_raster(T method, const G& grid) {
    if (method == "cartesian") {
      return std::unique_ptr<exactextract::AbstractRaster<double>>(
        static_cast<exactextract::AbstractRaster<double>*>(
          new exactextract::CartesianAreaRaster<double>(grid)));
    } else if (method == "spherical") {
      return std::unique_ptr<exactextract::AbstractRaster<double>>(
        static_cast<exactextract::AbstractRaster<double>*>(
          new exactextract::SphericalAreaRaster<double>(grid)));
    } else {
      Rcpp::stop("Unknown area method: " + method);
    }
}
