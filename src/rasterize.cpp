// Copyright (c) 2022 ISciences, LLC.
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

// [[Rcpp::plugins("cpp14")]]
#include <Rcpp.h>

#include "geos_r.h"
#include "raster_utils.h"

#include "exactextract/src/raster_cell_intersection.h"

// [[Rcpp::export]]
void CPP_update_max_coverage(Rcpp::NumericVector & extent,
		Rcpp::NumericVector & res,
		Rcpp::NumericMatrix & max_coverage,
		Rcpp::IntegerMatrix & max_coverage_index,
		Rcpp::NumericMatrix & tot_coverage,
		const Rcpp::RawVector & wkb,
		int index)
{
  GEOSAutoHandle geos;

  auto grid = make_grid(extent, res);

  auto coverage_fraction = exactextract::raster_cell_intersection(grid, geos.handle, read_wkb(geos.handle, wkb).get());

  auto ix = grid.row_offset(coverage_fraction.grid());
  auto jx = grid.col_offset(coverage_fraction.grid());

  for (size_t i = 0; i < coverage_fraction.rows(); i++) {
    for (size_t j = 0; j < coverage_fraction.cols(); j++) {
      auto cov = coverage_fraction(i, j);
      if (cov > 0) {
        tot_coverage(i + ix, j + jx) += cov;
        if (cov > max_coverage(i + ix, j + jx)) {
          max_coverage(i + ix, j + jx) = cov;
          max_coverage_index(i + ix, j + jx) = index;
        }
      }
    }
  }
}


