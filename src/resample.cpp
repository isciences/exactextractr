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

#include <Rcpp.h>

#include "s4_raster_source.h"
#include "raster_utils.h"

#include "exactextract/src/box.h"
#include "exactextract/src/raster_stats.h"

using exactextract::Box;
using exactextract::RasterStats;
using exactextract::RasterView;

// TODO merge with nearly-identical code in exact_extract.cpp
static double get_stat_value(const RasterStats<double> & stats, const std::string & stat_name) {
  if (stat_name == "mean") return stats.mean();

  else if (stat_name == "sum") return stats.sum();
  else if (stat_name == "count") return stats.count();

  else if (stat_name == "min") return stats.min().value_or(NA_REAL);
  else if (stat_name == "max") return stats.max().value_or(NA_REAL);

  else if (stat_name == "mode") return stats.mode().value_or(NA_REAL);
  else if (stat_name == "majority") return stats.mode().value_or(NA_REAL);
  else if (stat_name == "minority") return stats.minority().value_or(NA_REAL);

  else if (stat_name == "variety") return stats.variety();
  else if (stat_name == "weighted_mean") return stats.weighted_mean();
  else if (stat_name == "weighted_sum") return stats.weighted_sum();

  else if (stat_name == "variance") return stats.variance();
  else if (stat_name == "stdev") return stats.stdev();
  else if (stat_name == "coefficient_of_variation") return stats.coefficient_of_variation();

  else Rcpp::stop("Unknown stat: " + stat_name);
}

// [[Rcpp::export]]
Rcpp::S4 CPP_resample(Rcpp::S4 & rast_in,
                      Rcpp::S4 & rast_out,
                      Rcpp::Nullable<Rcpp::CharacterVector> p_stat,
                      Rcpp::Nullable<Rcpp::Function> p_fun,
                      bool coverage_area,
                      std::string area_method) {
  try {
    Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
    Rcpp::Environment xx = Rcpp::Environment::namespace_env("exactextractr");
    Rcpp::Function rasterFn = raster["raster"];
    Rcpp::Function valuesFn = xx[".setValues"];
    Rcpp::Function nlyrFn = xx[".numLayers"];

    int numLayers = Rcpp::as<int>(nlyrFn(rast_in));

    S4RasterSource rsrc(rast_in);

    Rcpp::S4 out = rasterFn(rast_out);

    auto grid_in = make_grid(rast_in);
    auto grid_out = make_grid(rast_out);

    std::string stat_name;
    bool store_values = false;
    bool r_summary_function = false;

    if (p_stat.isNotNull()) {
      Rcpp::CharacterVector stat = p_stat.get();
      stat_name = Rcpp::as<std::string>(stat[0]);
      store_values = requires_stored_values(stat_name);
    } else {
      r_summary_function = true;
    }

    Rcpp::NumericMatrix values_out = Rcpp::no_init(grid_out.rows(), grid_out.cols());

    std::vector<std::unique_ptr<NumericVectorRaster>> values(numLayers);

    for (size_t row = 0; row < grid_out.rows(); row++) {
      // Read enough source raster data to process an entire destination row at
      // a time, since getValuesBlock calls have a lot of overhead.
      auto y = grid_out.y_for_row(row);
      auto ymin = y - grid_out.dy();
      auto ymax = y + grid_out.dy();

      Box row_box{ grid_out.xmin(), ymin, grid_out.xmax(), ymax };

      for (int i = 0; i < numLayers; i++) {
        values[i] = std::unique_ptr<NumericVectorRaster>(
          static_cast<NumericVectorRaster*>(rsrc.read_box(row_box, i).release()));
      }

      for (size_t col = 0; col < grid_out.cols(); col++) {
        Box cell = grid_cell(grid_out, row, col);
        auto coverage_fraction = raster_cell_intersection(grid_in, cell);

        auto& cov_grid = coverage_fraction.grid();
        if (coverage_area) {
          auto areas = get_area_raster(area_method, cov_grid);
          for (size_t i = 0; i < coverage_fraction.rows(); i++) {
            for (size_t j = 0; j < coverage_fraction.cols(); j++) {
              coverage_fraction(i, j) = coverage_fraction(i, j) * (*areas)(i, j);
            }
          }
        }

        if (r_summary_function) {
          Rcpp::Function summary_fun = p_fun.get();

          // Transform values to same grid as coverage fractions
          auto coverage_vec = as_vector(coverage_fraction);

          Rcpp::NumericVector result;

          if (numLayers == 1) {
            RasterView<double> rt(*(values[0]), cov_grid);
            auto value_vec = as_vector(rt);
            result = summary_fun(value_vec, coverage_vec);
          } else {
            Rcpp::NumericMatrix value_mat = Rcpp::no_init(coverage_vec.size(), numLayers);
            for (int i = 0; i < numLayers; i++) {
              RasterView<double> rt(*(values[i]), cov_grid);
              auto value_vec = as_vector(rt);
              value_mat(Rcpp::_, i) = value_vec;
            }

            result = summary_fun(value_mat, coverage_vec);
          }

          if (result.size() != 1) {
            Rcpp::stop("Summary function must return a single value");
          }

          values_out(row, col) = result[0];
        } else {
          RasterStats<double> stats{store_values};

          if (!cov_grid.empty()) {
            stats.process(coverage_fraction, *(values[0]));
          }

          values_out(row, col) = get_stat_value(stats, stat_name);
        }
      }
    }

    out = valuesFn(out, values_out);
    return out;
  } catch (std::exception & e) {
    // throw predictable exception class
#ifdef __SUNPRO_CC
    // Rcpp::stop crashes CRAN Solaris build
    // https://github.com/RcppCore/Rcpp/issues/1159
    Rf_error(e.what());
    return R_NilValue;
#else
    Rcpp::stop(e.what());
#endif
  }
}
