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
#include <geos_c.h>

#include "exactextract/src/geos_utils.h"
#include "exactextract/src/grid.h"
#include "exactextract/src/matrix.h"
#include "exactextract/src/raster_cell_intersection.h"
#include "exactextract/src/raster_source.h"
#include "exactextract/src/raster_stats.h"

using geom_ptr= std::unique_ptr<GEOSGeometry, std::function<void(GEOSGeometry*)>>;
using wkb_reader_ptr = std::unique_ptr<GEOSWKBReader, std::function<void(GEOSWKBReader*)>>;

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

static Grid<bounded_extent> make_grid(const Rcpp::S4 & rast) {
  Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");

  Rcpp::S4 extent = rast.slot("extent");

  Rcpp::Function resFn = raster["res"];

  Rcpp::NumericVector res = resFn(rast);

  return {{
    extent.slot("xmin"),
    extent.slot("ymin"),
    extent.slot("xmax"),
    extent.slot("ymax"),
    },
    res[0],
    res[1]
  };
}

// Construct a Raster using an R vector for storage
class NumericVectorRaster : public exactextract::AbstractRaster<double> {
public:
  NumericVectorRaster(const Rcpp::NumericVector & vec, const Grid<bounded_extent> & g) :
    AbstractRaster<double>(g),
    m_vec(vec)
  {}

  double operator()(size_t row, size_t col) const final {
    return m_vec[row*cols() + col];
  }

private:
  const Rcpp::NumericVector m_vec;
};

// Read raster values from an R raster object
class S4RasterSource {
public:
  S4RasterSource(Rcpp::S4 rast) :
    m_grid(Grid<bounded_extent>::make_empty()),
    m_rast(rast),
    m_last_box(0, 0, 0, 0)
  {
    m_grid = make_grid(rast);
  }

  const Grid<bounded_extent> &grid() const {
    return m_grid;
  }

  std::unique_ptr<exactextract::AbstractRaster<double>> read_box(const exactextract::Box & box, int layer) {
    auto cropped_grid = m_grid.crop(box);

    if (!(box == m_last_box)) {
      m_last_box = box;

      Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
      Rcpp::Function getValuesBlockFn = raster["getValuesBlock"];

      if (cropped_grid.empty()) {
        m_rast_values = Rcpp::no_init(0, 0);
      } else {
        // Instead of reading only values for the requested band, we read values
        // for all requested bands and then cache them to return from subsequent
        // calls to read_box. There are two reasons for this:
        // 1) Use of the 'lyrs` argument to getValuesBlock does not work for
        //    a single band in the CRAN version of raster
        // 2) Reading everything once avoids the very significant per-call
        //    overhead of getValuesBlock, which largely comes from operations
        //    like setting names on the result vector.
        m_rast_values = getValuesBlockFn(m_rast,
                                         1 + cropped_grid.row_offset(m_grid),
                                         cropped_grid.rows(),
                                         1 + cropped_grid.col_offset(m_grid),
                                         cropped_grid.cols(),
                                         Rcpp::Named("format", "m"));
      }
    }

    if (cropped_grid.empty()) {
      return std::make_unique<NumericVectorRaster>(m_rast_values, cropped_grid);
    } else {
      return std::make_unique<NumericVectorRaster>(m_rast_values(Rcpp::_, layer),
                                                   cropped_grid);
    }
  }

private:
  Grid<bounded_extent> m_grid;
  Rcpp::S4 m_rast;
  Rcpp::NumericMatrix m_rast_values;
  exactextract::Box m_last_box;
};

// GEOS warning handler
static void geos_warn(const char* fmt, ...) {
  char buf[BUFSIZ] = { '\0' };

  va_list msg;
  va_start(msg, fmt);
  vsnprintf(buf, BUFSIZ*sizeof(char), fmt, msg);
  va_end(msg);

  Rcpp::Function warning("warning");
  warning(buf);
}

// GEOS error handler
static void geos_error(const char* fmt, ...) {
  char buf[BUFSIZ] = { '\0' };

  va_list msg;
  va_start(msg, fmt);
  vsnprintf(buf, BUFSIZ*sizeof(char), fmt, msg);
  va_end(msg);

  Rcpp::stop(buf);
}

// GEOSContextHandle wrapper to ensure finishGEOS is called.
struct GEOSAutoHandle {
  GEOSAutoHandle() {
    handle = initGEOS_r(geos_warn, geos_error);
  }

  ~GEOSAutoHandle() {
    finishGEOS_r(handle);
  }

  GEOSContextHandle_t handle;
};

// Return a smart pointer to a Geometry, given WKB input
static geom_ptr read_wkb(const GEOSContextHandle_t & context, const Rcpp::RawVector & wkb) {
  wkb_reader_ptr wkb_reader{ GEOSWKBReader_create_r(context), [context](GEOSWKBReader* r) { GEOSWKBReader_destroy_r(context, r); } };

  geom_ptr geom{ GEOSWKBReader_read_r(context,
                                      wkb_reader.get(),
                                      std::addressof(wkb[0]),
                                      wkb.size()),
                 [context](GEOSGeometry* g) { GEOSGeom_destroy_r(context, g); } };

  if (geom.get() == nullptr) {
    Rcpp::stop("Failed to parse WKB geometry");
  }

  return geom;
}

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
Rcpp::S4 CPP_coverage_fraction(Rcpp::S4 & rast, const Rcpp::RawVector & wkb, bool crop)
{
  GEOSAutoHandle geos;
  Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
  Rcpp::Function rasterFn = raster["raster"];
  Rcpp::Function crsFn = raster["crs"];

  auto grid = make_grid(rast);
  auto coverage_fraction = raster_cell_intersection(grid, geos.handle, read_wkb(geos.handle, wkb).get());

  if (crop) {
    grid = coverage_fraction.grid();
  }
  RasterView<float> coverage_view(coverage_fraction, grid);

  Rcpp::NumericMatrix weights{static_cast<int>(grid.rows()),
                              static_cast<int>(grid.cols())};

  for (size_t i = 0; i < grid.rows(); i++) {
    for (size_t j = 0; j < grid.cols(); j++) {
      weights(i, j) = coverage_view(i, j);
      if (!crop && std::isnan(weights(i, j))) {
        weights(i, j) = 0;
      }
    }
  }

  return rasterFn(weights,
                  Rcpp::Named("xmn")=grid.xmin(),
                  Rcpp::Named("xmx")=grid.xmax(),
                  Rcpp::Named("ymn")=grid.ymin(),
                  Rcpp::Named("ymx")=grid.ymax(),
                  Rcpp::Named("crs")=crsFn(rast));
}

static int get_nlayers(Rcpp::S4 & rast) {
  Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
  Rcpp::Function nlayersFn = raster["nlayers"];

  Rcpp::NumericVector nlayersVec = nlayersFn(rast);

  return static_cast<int>(nlayersVec[0]);
}

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

// [[Rcpp::export]]
Rcpp::S4 CPP_resample(Rcpp::S4 & rast_in,
                      Rcpp::S4 & rast_out,
                      const Rcpp::StringVector & stat) {
  Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
  Rcpp::Function rasterFn = raster["raster"];
  Rcpp::Function valuesFn = raster["values<-"];

  if (stat.size() != 1) {
    Rcpp::stop("Only a single operation may be used for resampling.");
  }

  S4RasterSource rsrc(rast_in);

  Rcpp::S4 out = rasterFn(rast_out);

  auto grid_in = make_grid(rast_in);
  auto grid_out = make_grid(rast_out);

  std::string stat_name = Rcpp::as<std::string>(stat[0]);

  Rcpp::NumericMatrix values_out = Rcpp::no_init(grid_out.rows(), grid_out.cols());

  for (size_t row = 0; row < grid_out.rows(); row++) {
    // Read enough source raster data to process an entire destination row at
    // a time, since getValuesBlock calls have a lot of overhead.
    auto y = grid_out.y_for_row(row);
    auto ymin = y - grid_out.dy();
    auto ymax = y + grid_out.dy();

    Box row_box{ grid_out.xmin(), ymin, grid_out.xmax(), ymax };
    auto values = rsrc.read_box(row_box, 0);

    for (size_t col = 0; col < grid_out.cols(); col++) {
      RasterStats<double> stats;

      Box cell = grid_cell(grid_out, row, col);
      auto coverage_fraction = raster_cell_intersection(grid_in, cell);

      auto& cov_grid = coverage_fraction.grid();

      if (!cov_grid.empty()) {
        stats.process(coverage_fraction, *values);
      }

      values_out(row, col) = get_stat_value(stats, stat_name);
    }
  }

  out = valuesFn(out, values_out);
  return out;
}
