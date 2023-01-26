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

// [[Rcpp::plugins("cpp14")]]
#include <Rcpp.h>

#include <memory>

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

#include <vector>
#include <string>

// [[Rcpp::export]]
Rcpp::List CPP_exact_extract(Rcpp::S4 & rast,
                             const Rcpp::NumericVector & rast_ext,
                             const Rcpp::NumericVector & rast_res,
                             Rcpp::Nullable<Rcpp::S4> & rast_uncropped,
                             Rcpp::Nullable<Rcpp::S4> & weights,
                             const Rcpp::RawVector & wkb,
                             double default_value,
                             double default_weight,
                             bool include_xy,
                             bool include_cell_number,
                             bool include_area,
                             bool area_weights,
                             bool coverage_areas,
                             Rcpp::Nullable<Rcpp::CharacterVector> & p_area_method,
                             Rcpp::Nullable<Rcpp::List> & include_cols,
                             Rcpp::CharacterVector & src_names,
                             Rcpp::Nullable<Rcpp::CharacterVector> & p_weights_names,
                             bool warn_on_disaggregate,
                             double grid_compat_tol) {
  try {
    GEOSAutoHandle geos;
    Rcpp::Function names("names");

    auto grid = make_grid(rast);
    auto weights_grid = exactextract::Grid<bounded_extent>::make_empty();
    auto common_grid = grid;

    S4RasterSource rsrc(rast, rast_ext, rast_res, default_value);
    int src_nlayers = get_nlayers(rast);

    std::unique_ptr<S4RasterSource> rweights;
    int weights_nlayers = 0;
    Rcpp::CharacterVector weights_names;

    std::string area_method;
    if (p_area_method.isNotNull()) {
      Rcpp::CharacterVector amethod = p_area_method.get();
      area_method = amethod[0];
    }

    if (weights.isNotNull()) {
      Rcpp::S4 weights_s4 = weights.get();
      weights_nlayers = get_nlayers(weights_s4);
      weights_grid = make_grid(weights_s4);

      common_grid = grid.common_grid(weights_grid, grid_compat_tol);

      rweights = std::make_unique<S4RasterSource>(weights_s4, default_weight);
      weights_names = p_weights_names.get();

      if (warn_on_disaggregate && (common_grid.dx() < grid.dx() || common_grid.dy() < grid.dy())) {
        Rcpp::warning("value raster implicitly disaggregated to match higher resolution of weights");
      }
    } else if (area_weights) {
      weights_nlayers = 1;
      weights_names = p_weights_names.get();
    }

    auto geom = read_wkb(geos.handle, wkb);

    auto bbox = exactextract::geos_get_box(geos.handle, geom.get());

    common_grid = common_grid.crop(bbox);

    auto coverage_fractions = raster_cell_intersection(common_grid, geos.handle, geom.get());
    auto& cov_grid = coverage_fractions.grid();

    // We can't construct an Rcpp::DataFrame with a variable number of columns
    // because of a bug in Rcpp (see https://github.com/RcppCore/Rcpp/pull/1099)
    // Instead, we build up a list, and then create a data frame from the list.
    // Profiling shows that this ends up in as.data.frame, which is slow.
    // Once Rcpp 1.0.6 is released, this can be reworked to be simpler and faster.
    Rcpp::List cols;

    Rcpp::NumericVector coverage_vec = as_vector(coverage_fractions);
    if (coverage_areas) {
        auto areas = get_area_raster(area_method, cov_grid);
        Rcpp::NumericVector area_vec = as_vector(*areas);
        coverage_vec = coverage_vec * area_vec;
    }
    Rcpp::LogicalVector covered = coverage_vec > 0;

    if (include_cols.isNotNull()) {
      Rcpp::List include_cols_list = include_cols.get();
      Rcpp::CharacterVector include_names = include_cols_list.attr("names");

      for (int i = 0; i < include_names.size(); i++) {
        std::string name(include_names[i]);
        cols[name] = include_cols_list[name];
      }
    }

    for (int i = 0; i < src_nlayers; i++) {
      auto values = rsrc.read_box(cov_grid.extent(), i);
      const NumericVectorRaster* r = static_cast<NumericVectorRaster*>(values.get());

      // TODO Perhaps extend this to preserve types (integer, logical.)
      // A bit challenging, since we don't know the type of the raster
      // until we call getValuesBlock, and even then we can't be sure
      // that the returned type is correct (as in a RasterStack with mixed types.)
      // Since R integers are only 32-bit, we are not going to lose data by
      // converting everything to numeric, although we pay a storage penalty.
      Rcpp::NumericVector value_vec = r->vec();
      if (grid.dx() != cov_grid.dx() || grid.dy() != cov_grid.dy() ||
          value_vec.size() != covered.size()) {
        // Transform values to same grid as coverage fractions
        RasterView<double> rt(*r, cov_grid);
        value_vec = as_vector(rt);
      }

      value_vec = value_vec[covered];
      cols[std::string(src_names[i])] = value_vec;
    }

    for (int i = 0; i < weights_nlayers; i++) {
      Rcpp::NumericVector weight_vec;
      if (area_weights) {
        auto weights = get_area_raster(area_method, cov_grid);
        weight_vec = as_vector(*weights);
      } else {
        auto values = rweights->read_box(cov_grid.extent(), i);
        const NumericVectorRaster* r = static_cast<NumericVectorRaster*>(values.get());
        weight_vec = r->vec();

        if (weights_grid.dx() != cov_grid.dx() || weights_grid.dy() != cov_grid.dy() ||
            weight_vec.size() != covered.size()) {
          // Transform weights to same grid as coverage fractions
          RasterView<double> rt (*r, cov_grid);

          weight_vec = as_vector(rt);
        }
      }

      weight_vec = weight_vec[covered];

      std::string colname(weights_names[i]);
      if (cols.containsElementNamed(colname.c_str())) {
        // append ".1" to the column name, to match the behavior of
        // data.frame(). We're safe to just add ".1" (instead of incrementing
        // .1 to .2, for example) because duplicated names within the values
        // or weight stack will already have been made unique by raster::stack()
        colname = colname + ".1"; // append .1
      }
      cols[colname] = weight_vec;
    }

    if (include_xy) {
      // Include xy values from whichever input raster has the higher resolution
      if (weights_nlayers > 0 && (weights_grid.dx() < grid.dx() || weights_grid.dy() < grid.dy())) {
        Rcpp::S4 weights_s4 = weights.get();
        cols["x"] = get_x_values(weights_s4, cov_grid)[covered];
        cols["y"] = get_y_values(weights_s4, cov_grid)[covered];
      } else {
        cols["x"] = get_x_values(rast, cov_grid)[covered];
        cols["y"] = get_y_values(rast, cov_grid)[covered];
      }
    }

    if (include_cell_number) {
      if (rast_uncropped.isNull()) {
        cols["cell"] = get_cell_numbers(rast, cov_grid)[covered];
      } else {
        Rcpp::S4 tmp = rast_uncropped.get();
        cols["cell"] = get_cell_numbers(tmp, cov_grid)[covered];
      }
    }

    if (include_area) {
      auto area_rast = get_area_raster(area_method, cov_grid);
      Rcpp::NumericVector area_vec = as_vector(*area_rast);
      cols["area"] = area_vec[covered];
    }

    if (coverage_areas) {
      cols["coverage_area"] = coverage_vec[covered];
    } else {
      cols["coverage_fraction"] = coverage_vec[covered];
    }

    return cols;
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

enum class WeightingMethod {
  NONE,
  RASTER,
  AREA
};


static int get_num_stats(const Rcpp::StringVector & stats,
                  const std::size_t num_quantiles,
                  const std::size_t num_unique_values) {
  int num_stats = 0;
  for (const auto & stat : stats) {
      if (stat == std::string("frac") ||
          stat == std::string("weighted_frac")) {
        num_stats += static_cast<int>(num_unique_values);
      } else if (stat == std::string("quantile")) {
        num_stats += static_cast<int>(num_quantiles);
      } else {
        num_stats += 1;
      }
  }

  return num_stats;
}

// Return a matrix with one row per stat and one row per raster layer
// [[Rcpp::export]]
Rcpp::NumericMatrix CPP_stats(Rcpp::S4 & rast,
                              const Rcpp::NumericVector & rast_ext,
                              const Rcpp::NumericVector & rast_res,
                              Rcpp::Nullable<Rcpp::S4> weights,
                              const Rcpp::RawVector & wkb,
                              double default_value,
                              double default_weight,
                              bool coverage_areas,
                              Rcpp::Nullable<Rcpp::CharacterVector> & p_area_method,
                              const Rcpp::StringVector & stats,
                              int max_cells_in_memory,
                              double grid_compat_tol,
                              const Rcpp::Nullable<Rcpp::NumericVector> & quantiles) {
  try {
    GEOSAutoHandle geos;

    if (max_cells_in_memory < 1) {
      Rcpp::stop("Invalid value for max_cells_in_memory: %d", max_cells_in_memory);
    }

    int nlayers = get_nlayers(rast);

    S4RasterSource rsrc(rast, rast_ext, rast_res, default_value);

    std::unique_ptr<S4RasterSource> rweights;
    std::string area_method;
    if (p_area_method.isNotNull()) {
      area_method = ((Rcpp::CharacterVector) p_area_method.get())[0];
    }

    WeightingMethod weighting = WeightingMethod::NONE;
    int nweights = 0;
    if (weights.isNotNull()) {
      weighting = WeightingMethod::RASTER;
      Rcpp::S4 weights_s4 = weights.get();
      nweights = get_nlayers(weights_s4);

      if (nlayers > 1 && nweights > 1 && nlayers != nweights) {
        Rcpp::stop("Incompatible number of layers in value and weighting rasters");
      }

      rweights = std::make_unique<S4RasterSource>(weights_s4, default_weight);
    } else if (p_area_method.isNotNull()) {
      weighting = WeightingMethod::AREA;
      nweights = 1;
    }

    auto geom = read_wkb(geos.handle, wkb);
    auto bbox = exactextract::geos_get_box(geos.handle, geom.get());

    auto grid = weighting == WeightingMethod::RASTER ?
      rsrc.grid().common_grid(rweights->grid(), grid_compat_tol) : rsrc.grid();

    bool disaggregated = (grid.dx() < rsrc.grid().dx() || grid.dy() < rsrc.grid().dy());

    bool store_values = false;
    bool calc_value_set = false;
    std::size_t num_unique_values = 0;
    std::size_t num_quantiles = 0;

    for (const auto & stat : stats) {
      store_values = store_values || requires_stored_values(stat);

      if (disaggregated && (stat == std::string("count") ||
                            stat == std::string("sum"))) {
        Rcpp::stop("Cannot compute 'count' or 'sum' when value raster is disaggregated to resolution of weights.");
      }

      if (stat == std::string("frac") ||
          stat == std::string("weighted_frac")) {
        calc_value_set = true;
      }

      if (stat == std::string("quantile")) {
        if (quantiles.isNotNull()) {
          Rcpp::NumericVector qvec = quantiles.get();
          num_quantiles = qvec.size();
        }

        if (num_quantiles == 0) {
          Rcpp::stop("Quantiles not specified.");
        }
      }
    }

    int nresults = std::max(nlayers, nweights);

    std::vector<RasterStats<double>> raster_stats;
    raster_stats.reserve(nresults);
    for (int i = 0; i < nresults; i++) {
      raster_stats.emplace_back(store_values);
    }

    if (bbox.intersects(grid.extent())) {
      auto cropped_grid = grid.crop(bbox);

      for (const auto &subgrid : subdivide(cropped_grid, max_cells_in_memory)) {
        auto coverage_fraction = raster_cell_intersection(subgrid, geos.handle, geom.get());
        if (coverage_areas) {
          auto areas = get_area_raster(area_method, subgrid);
          for (size_t i = 0; i < coverage_fraction.rows(); i++) {
            for (size_t j = 0; j < coverage_fraction.cols(); j++) {
              coverage_fraction(i, j) = coverage_fraction(i, j) * (*areas)(i, j);
            }
          }
        }

        auto& cov_grid = coverage_fraction.grid();

        if (!cov_grid.empty()) {
          if (weighting == WeightingMethod::NONE) {
            for (int i = 0; i < nlayers; i++) {
              auto values = rsrc.read_box(cov_grid.extent(), i);
              raster_stats[i].process(coverage_fraction, *values);
            }
          } else if (weighting == WeightingMethod::AREA) {
            auto weights = get_area_raster(area_method, cov_grid);

            for (int i = 0; i < nlayers; i++) {
              auto values = rsrc.read_box(cov_grid.extent(), i);
              raster_stats[i].process(coverage_fraction, *values, *weights);
            }
          } else {
            if (nlayers > nweights) {
              // recycle weights
              auto weights = rweights->read_box(cov_grid.extent(), 0);

              for (int i = 0; i < nlayers; i++) {
                auto values = rsrc.read_box(cov_grid.extent(), i);
                raster_stats[i].process(coverage_fraction, *values, *weights);
              }
            } else if (nweights > nlayers) {
              // recycle values
              auto values = rsrc.read_box(cov_grid.extent(), 0);

              for (int i = 0; i < nweights; i++) {
                auto weights = rweights->read_box(cov_grid.extent(), i);
                raster_stats[i].process(coverage_fraction, *values, *weights);
              }
            } else {
              // process values and weights in parallel
              for (int i = 0; i < nlayers; i++) {
                auto values = rsrc.read_box(cov_grid.extent(), i);
                auto weights = rweights->read_box(cov_grid.extent(), i);

                raster_stats[i].process(coverage_fraction, *values, *weights);
              }
            }
          }
        }
      }
    }


    std::set<double> value_set;
    if (calc_value_set) {
      for (const auto& rs : raster_stats) {
        for (const auto& value : rs) {
          value_set.insert(value);
        }
      }
      num_unique_values = value_set.size();
    }

    auto stat_result_size = get_num_stats(stats, num_quantiles, num_unique_values);
    Rcpp::NumericMatrix stat_results = Rcpp::no_init(nresults, stat_result_size);

    if (calc_value_set) {
      Rcpp::NumericVector unique_values(value_set.begin(), value_set.end());
      stat_results.attr("unique_values") = unique_values;
    }

    // rows (j) represent layers
    // cols (i) represent stats
    for (int j = 0; j < nresults; j++) {
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

        else if (stat == std::string("weighted_variance")) stat_results(j, i++) = rs.weighted_variance();
        else if (stat == std::string("weighted_stdev")) stat_results(j, i++) = rs.weighted_stdev();

        else if (stat == std::string("quantile")) {
          Rcpp::NumericVector qvec = quantiles.get();
          for (double q : qvec) {
            stat_results(j, i++) = rs.quantile(q).value_or(NA_REAL);
          }
        }

        else if (stat == std::string("frac")) {
          for (double v : value_set) {
            stat_results(j, i++) = rs.frac(v).value_or(0);
          }
        }

        else if (stat == std::string("weighted_frac")) {
          for (double v : value_set) {
            stat_results(j, i++) = rs.weighted_frac(v).value_or(0);
          }
        }

        else Rcpp::stop("Unknown stat: " + stat);
      }
    }

    // FIXME set dimnames here

    return stat_results;
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
