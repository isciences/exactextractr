// Copyright (c) 2018-2020 ISciences, LLC.
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
#include <memory>

#include <Rcpp.h>

#include "geos_r.h"
#include "raster_utils.h"

#include "s4_raster_source.h"

#include "exactextract/src/geos_utils.h"
#include "exactextract/src/grid.h"
#include "exactextract/src/matrix.h"
#include "exactextract/src/raster_cell_intersection.h"
#include "exactextract/src/raster_source.h"
#include "exactextract/src/raster_stats.h"

using exactextract::Box;
using exactextract::Matrix;
using exactextract::Grid;
using exactextract::subdivide;
using exactextract::bounded_extent;
using exactextract::Raster;
using exactextract::RasterView;
using exactextract::raster_cell_intersection;
using exactextract::RasterStats;
using exactextract::RasterSource;


// [[Rcpp::export]]
Rcpp::List CPP_exact_extract(Rcpp::S4 & rast, const Rcpp::RawVector & wkb) {
  GEOSAutoHandle geos;

  auto grid = make_grid(rast);
  auto coverage_fractions = raster_cell_intersection(grid, geos.handle, read_wkb(geos.handle, wkb).get());

  size_t nrow = coverage_fractions.rows();
  size_t ncol = coverage_fractions.cols();

  Rcpp::NumericMatrix weights = Rcpp::no_init(nrow, ncol);
  for (size_t i = 0; i < nrow; i++) {
    for (size_t j = 0; j < ncol; j++) {
      weights(i, j) = coverage_fractions(i ,j);
    }
  }

  if (nrow > 0) {
    size_t row_us = 1 + coverage_fractions.grid().row_offset(grid);
    if (row_us > static_cast<size_t>(std::numeric_limits<int>::max())) {

    }
  }

  int row = NA_INTEGER;
  int col = NA_INTEGER;
  if (nrow > 0) {
    size_t row_us = (1 + coverage_fractions.grid().row_offset(grid));
    if (row_us > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("Cannot represent row offset as an R integer");
    }
    row = static_cast<int>(row_us);
  }
  if (ncol > 0) {
    size_t col_us = (1 + coverage_fractions.grid().col_offset(grid));
    if (col_us > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("Cannot represent column offset as an R integer");
    }
    col = static_cast<int>(col_us);
  }

  return Rcpp::List::create(
    Rcpp::Named("row")     = row,
    Rcpp::Named("col")     = col,
    Rcpp::Named("weights") = weights
  );
}


// [[Rcpp::export]]
Rcpp::NumericMatrix CPP_stats(Rcpp::S4 & rast,
                              Rcpp::Nullable<Rcpp::S4> weights,
                              const Rcpp::RawVector & wkb,
                              const Rcpp::StringVector & stats,
                              int max_cells_in_memory,
                              const Rcpp::Nullable<Rcpp::NumericVector> & quantiles) {
  try {
    GEOSAutoHandle geos;

    if (max_cells_in_memory < 1) {
      Rcpp::stop("Invalid value for max_cells_in_memory: %d", max_cells_in_memory);
    }

    int nlayers = get_nlayers(rast);

    S4RasterSource rsrc(rast);

    std::unique_ptr<S4RasterSource> rweights;
    bool weighted = false;
    if (weights.isNotNull()) {
      Rcpp::S4 weights_s4 = weights.get();

      if (get_nlayers(weights_s4) != 1) {
        Rcpp::stop("Weighting raster must have only a single layer.");
      }

      rweights = std::make_unique<S4RasterSource>(weights_s4);
      weighted = true;
    }

    bool store_values = false;
    int stat_result_size = 0;
    for (const auto & stat : stats) {
      // explicit construction of std::string seems necessary to avoid ambiguous overload error
      if (stat == std::string("mode") ||
          stat == std::string("majority") ||
          stat == std::string("minority") ||
          stat == std::string("variety") ||
          stat == std::string("median") ||
          stat == std::string("quantile")) {
        store_values = true;
      }

      if (!weighted &&
          (stat == std::string("weighted_mean") ||
          stat == std::string("weighted_sum"))) {
        Rcpp::stop("Weighted stat requested but no weights provided.");
      }

      if (stat == std::string("quantile")) {
        int num_quantiles = 0;

        if (quantiles.isNotNull()) {
          Rcpp::NumericVector qvec = quantiles.get();
          num_quantiles = qvec.size();
        }

        if (num_quantiles == 0) {
          Rcpp::stop("Quantiles not specified.");
        }

        stat_result_size += num_quantiles;
      } else {
        stat_result_size += 1;
      }
    }

    std::vector<RasterStats<double>> raster_stats;
    raster_stats.reserve(nlayers);
    for (int i = 0; i < nlayers; i++) {
      raster_stats.emplace_back(store_values);
    }

    Rcpp::NumericMatrix stat_results = Rcpp::no_init(nlayers, stat_result_size);

    auto geom = read_wkb(geos.handle, wkb);
    auto bbox = exactextract::geos_get_box(geos.handle, geom.get());

    auto grid = weighted ? rsrc.grid().common_grid(rweights->grid()) : rsrc.grid();

    if (bbox.intersects(grid.extent())) {
      auto cropped_grid = grid.crop(bbox);

      for (const auto &subgrid : subdivide(cropped_grid, max_cells_in_memory)) {
        auto coverage_fraction = raster_cell_intersection(subgrid, geos.handle, geom.get());
        auto& cov_grid = coverage_fraction.grid();

        if (!cov_grid.empty()) {
          if (weighted) {
            auto weights = rweights->read_box(cov_grid.extent(), 0);

            for (int i = 0; i < nlayers; i++) {
              auto values = rsrc.read_box(cov_grid.extent(), i);
              raster_stats[i].process(coverage_fraction, *values, *weights);
            }
          } else {
            for (int i = 0; i < nlayers; i++) {
              auto values = rsrc.read_box(cov_grid.extent(), i);
              raster_stats[i].process(coverage_fraction, *values);
            }
          }
        }
      }
    }

    for (int j = 0; j < nlayers; j++) {
      const auto& rs = raster_stats[j];

      int i = 0;
      for(const auto& stat : stats) {
        if (stat == std::string("mean")) stat_results(j, i++) = rs.mean();

        else if (stat == std::string("sum")) stat_results(j, i++) = rs.sum();
        else if (stat == std::string("count")) stat_results(j, i++) = rs.count();

        else if (stat == std::string("min")) stat_results(j, i++) = rs.min().value_or(NA_REAL);
        else if (stat == std::string("max")) stat_results(j, i++) = rs.max().value_or(NA_REAL);

        else if (stat == std::string("median")) stat_results(j, i++) = rs.quantile(0.5).value_or(NA_REAL);

        else if (stat == std::string("mode")) stat_results(j, i++) = rs.mode().value_or(NA_REAL);
        else if (stat == std::string("majority")) stat_results(j, i++) = rs.mode().value_or(NA_REAL);
        else if (stat == std::string("minority")) stat_results(j, i++) = rs.minority().value_or(NA_REAL);

        else if (stat == std::string("variety")) stat_results(j, i++) = rs.variety();
        else if (stat == std::string("weighted_mean")) stat_results(j, i++) = rs.weighted_mean();
        else if (stat == std::string("weighted_sum")) stat_results(j, i++) = rs.weighted_sum();

        else if (stat == std::string("variance")) stat_results(j, i++) = rs.variance();
        else if (stat == std::string("stdev")) stat_results(j, i++) = rs.stdev();
        else if (stat == std::string("coefficient_of_variation")) stat_results(j, i++) = rs.coefficient_of_variation();

        else if (stat == std::string("quantile")) {
          Rcpp::NumericVector qvec = quantiles.get();
          for (double q : qvec) {
            stat_results(j, i++) = rs.quantile(q).value_or(NA_REAL);
          }
        }

        else Rcpp::stop("Unknown stat: " + stat);
      }
    }

    return stat_results;
  } catch (std::exception & e) {
    // throw predictable exception class
    Rcpp::stop(e.what());
  }
}
