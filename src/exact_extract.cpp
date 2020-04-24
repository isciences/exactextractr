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

using exactextract::Matrix;
using exactextract::Grid;
using exactextract::subdivide;
using exactextract::bounded_extent;
using exactextract::Raster;
using exactextract::RasterView;
using exactextract::raster_cell_intersection;
using exactextract::RasterStats;
using exactextract::RasterSource;

// Construct a Raster using an R matrix for storage
class NumericMatrixRaster : public exactextract::AbstractRaster<double> {
public:
  NumericMatrixRaster(const Rcpp::NumericMatrix & mat, const Grid<bounded_extent> & g) :
    AbstractRaster<double>(g),
    m_mat(mat)
  {}

  double operator()(size_t row, size_t col) const final {
    return m_mat(row, col);
  }

private:
  const Rcpp::NumericMatrix m_mat;
};

// Read raster values from an R raster object
class S4RasterSource : public RasterSource {
public:
  S4RasterSource(Rcpp::S4 rast) : m_grid(Grid<bounded_extent>::make_empty()), m_rast(rast) {
    Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
    Rcpp::Function extentFn = raster["extent"];
    Rcpp::Function resFn = raster["res"];

    Rcpp::S4 extent = extentFn(rast);
    Rcpp::NumericVector res = resFn(rast);

    m_grid = {{
      extent.slot("xmin"),
      extent.slot("ymin"),
      extent.slot("xmax"),
      extent.slot("ymax"),
      },
      res[0],
      res[1]
    };
  }

  const Grid<bounded_extent> &grid() const override {
    return m_grid;
  }

  std::unique_ptr<exactextract::AbstractRaster<double>> read_box(const exactextract::Box & box) override {
    Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
    Rcpp::Function getValuesBlockFn = raster["getValuesBlock"];

    auto cropped_grid = m_grid.crop(box);

    Rcpp::NumericMatrix rast_values;

    if (cropped_grid.empty()) {
      rast_values = Rcpp::no_init(0, 0);
    } else {
      rast_values = getValuesBlockFn(m_rast,
                                     1 + cropped_grid.row_offset(m_grid),
                                     cropped_grid.rows(),
                                     1 + cropped_grid.col_offset(m_grid),
                                     cropped_grid.cols(),
                                     "matrix");
    }

    return std::make_unique<NumericMatrixRaster>(rast_values, cropped_grid);
  }

private:
  Grid<bounded_extent> m_grid;
  Rcpp::S4 m_rast;
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

static Grid<bounded_extent> make_grid(const Rcpp::S4 & rast) {
  Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");

  Rcpp::Function extentFn = raster["extent"];
  Rcpp::Function resFn = raster["res"];

  Rcpp::S4 extent = extentFn(rast);
  Rcpp::NumericVector res = resFn(rast);

  return {{ extent.slot("xmin"), extent.slot("ymin"), extent.slot("xmax"), extent.slot("ymax") },
            res[0], res[1] };
}

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
    if (row_us > std::numeric_limits<int>::max()) {

    }
  }

  int row = NA_INTEGER;
  int col = NA_INTEGER;
  if (nrow > 0) {
    size_t row_us = (1 + coverage_fractions.grid().row_offset(grid));
    if (row_us > std::numeric_limits<int>::max()) {
      throw std::runtime_error("Cannot represent row offset as an R integer");
    }
    row = static_cast<int>(row_us);
  }
  if (ncol > 0) {
    size_t col_us = (1 + coverage_fractions.grid().col_offset(grid));
    if (col_us > std::numeric_limits<int>::max()) {
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

static Rcpp::S4 layer(Rcpp::S4 & rast, int n) {
  Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
  Rcpp::Function subsetFn = raster["subset"];

  return subsetFn(rast, n+1);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix CPP_stats(Rcpp::S4 & rast,
                              Rcpp::Nullable<Rcpp::S4> weights,
                              const Rcpp::RawVector & wkb,
                              const Rcpp::StringVector & stats,
                              int max_cells_in_memory) {
  try {
    GEOSAutoHandle geos;

    if (max_cells_in_memory < 1) {
      Rcpp::stop("Invalid value for max_cells_in_memory: %d", max_cells_in_memory);
    }

    int nlayers = get_nlayers(rast);

    std::vector<S4RasterSource> rsrc;
    rsrc.reserve(nlayers);
    for (int i = 0; i < nlayers; i++) {
      rsrc.emplace_back(layer(rast, i));
    }

    std::unique_ptr<S4RasterSource> rweights;
    bool weighted = false;
    if (weights.isNotNull()) {
      Rcpp::S4 weights_s4 = weights.get();

      if (get_nlayers(weights_s4) != 1) {
        Rcpp::stop("Weighting raster must have only a single layer.");
      }

      rweights = std::make_unique<S4RasterSource>(layer(weights_s4, 0));
      weighted = true;
    }

    bool store_values = false;

    for (const auto & stat : stats) {
      // explicit construction of std::string seems necessary to avoid ambiguous overload error
      if (stat == std::string("mode") ||
          stat == std::string("majority") ||
          stat == std::string("minority") ||
          stat == std::string("variety")) {
        store_values = true;
      }

      if (!weighted &&
          (stat == std::string("weighted_mean") ||
          stat == std::string("weighted_sum"))) {
        Rcpp::stop("Weighted stat requested but no weights provided.");
      }
    }

    std::vector<RasterStats<double>> raster_stats;
    raster_stats.reserve(nlayers);
    for (int i = 0; i < nlayers; i++) {
      raster_stats.emplace_back(store_values);
    }

    Rcpp::NumericMatrix stat_results = Rcpp::no_init(nlayers, stats.size());

    auto geom = read_wkb(geos.handle, wkb);
    auto bbox = exactextract::geos_get_box(geos.handle, geom.get());

    auto grid = weighted ? rsrc[0].grid().common_grid(rweights->grid()) : rsrc[0].grid();

    if (bbox.intersects(grid.extent())) {
      auto cropped_grid = grid.crop(bbox);

      for (const auto &subgrid : subdivide(cropped_grid, max_cells_in_memory)) {
        auto coverage_fraction = raster_cell_intersection(subgrid, geos.handle, geom.get());
        auto& cov_grid = coverage_fraction.grid();

        if (!cov_grid.empty()) {
          if (weighted) {
            auto weights = rweights->read_box(cov_grid.extent());

            for (int i = 0; i < nlayers; i++) {
              auto values = rsrc[i].read_box(cov_grid.extent());
              raster_stats[i].process(coverage_fraction, *values, *weights);
            }
          } else {
            for (int i = 0; i < nlayers; i++) {
              auto values = rsrc[i].read_box(cov_grid.extent());
              raster_stats[i].process(coverage_fraction, *values);
            }
          }
        }
      }
    }

    for (int j = 0; j < nlayers; j++) {
      for(int i = 0; i < stats.size(); i++) {
        if (stats[i] == std::string("mean")) stat_results(j, i) = raster_stats[j].mean();

        else if (stats[i] == std::string("sum")) stat_results(j, i) = raster_stats[j].sum();
        else if (stats[i] == std::string("count")) stat_results(j, i) = raster_stats[j].count();

        else if (stats[i] == std::string("min")) stat_results(j, i) = raster_stats[j].min().value_or(NA_REAL);
        else if (stats[i] == std::string("max")) stat_results(j, i) = raster_stats[j].max().value_or(NA_REAL);

        else if (stats[i] == std::string("mode")) stat_results(j, i) = raster_stats[j].mode().value_or(NA_REAL);
        else if (stats[i] == std::string("majority")) stat_results(j, i) = raster_stats[j].mode().value_or(NA_REAL);
        else if (stats[i] == std::string("minority")) stat_results(j, i) = raster_stats[j].minority().value_or(NA_REAL);

        else if (stats[i] == std::string("variety")) stat_results(j, i) = raster_stats[j].variety();
        else if (stats[i] == std::string("weighted_mean")) stat_results(j, i) = raster_stats[j].weighted_mean();
        else if (stats[i] == std::string("weighted_sum")) stat_results(j, i) = raster_stats[j].weighted_sum();

        else Rcpp::stop("Unknown stat: " + stats[i]);
      }
    }

    return stat_results;
  } catch (std::exception & e) {
    // throw predictible exception class
    Rcpp::stop(e.what());
  }
}
