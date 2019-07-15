// Copyright (c) 2018-2019 ISciences, LLC.
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

#include "exactextract/src/grid.h"
#include "exactextract/src/matrix.h"
#include "exactextract/src/raster_cell_intersection.h"
#include "exactextract/src/raster_stats.h"

using geom_ptr= std::unique_ptr<GEOSGeometry, std::function<void(GEOSGeometry*)>>;
using wkb_reader_ptr = std::unique_ptr<GEOSWKBReader, std::function<void(GEOSWKBReader*)>>;

using exactextract::Matrix;
using exactextract::Grid;
using exactextract::bounded_extent;
using exactextract::Raster;
using exactextract::RasterView;
using exactextract::raster_cell_intersection;
using exactextract::RasterStats;

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
  const Rcpp::NumericMatrix& m_mat;
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

// Create a Grid from vectors representing the spatial extent and resolution
static Grid<bounded_extent> make_grid(const Rcpp::NumericVector & extent, const Rcpp::NumericVector & res) {
  double xmin = extent[0];
  double xmax = extent[1];
  double ymin = extent[2];
  double ymax = extent[3];

  double dx = res[0];
  double dy = res[1];

  return {{xmin, ymin, xmax, ymax}, dx, dy};
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
Rcpp::List CPP_exact_extract(const Rcpp::NumericVector & extent,
                             const Rcpp::NumericVector & res,
                             const Rcpp::RawVector & wkb) {
  GEOSAutoHandle geos;

  auto grid = make_grid(extent, res);
  auto coverage_fractions = raster_cell_intersection(grid, geos.handle, read_wkb(geos.handle, wkb).get());

  size_t nrow = coverage_fractions.rows();
  size_t ncol = coverage_fractions.cols();

  Rcpp::NumericMatrix weights = Rcpp::no_init(nrow, ncol);
  for (size_t i = 0; i < nrow; i++) {
    for (size_t j = 0; j < ncol; j++) {
      weights(i, j) = coverage_fractions(i ,j);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("row")     = 1 + coverage_fractions.grid().row_offset(grid),
    Rcpp::Named("col")     = 1 + coverage_fractions.grid().col_offset(grid),
    Rcpp::Named("weights") = weights
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix CPP_coverage_fraction(const Rcpp::NumericVector & extent,
                                          const Rcpp::NumericVector & res,
                                          const Rcpp::RawVector & wkb)
{
  GEOSAutoHandle geos;

  auto grid = make_grid(extent, res);
  auto coverage_fraction = raster_cell_intersection(grid, geos.handle, read_wkb(geos.handle, wkb).release());

  RasterView<float> coverage_view(coverage_fraction, grid);

  Rcpp::NumericMatrix weights{static_cast<int>(grid.rows()),
                              static_cast<int>(grid.cols())};

  for (size_t i = 0; i < grid.rows(); i++) {
    for (size_t j = 0; j < grid.cols(); j++) {
      weights(i, j) = coverage_view(i, j);
    }
  }

  return weights;
}

// [[Rcpp::export]]
SEXP CPP_stats(Rcpp::S4 & rast, const Rcpp::RawVector & wkb, const Rcpp::StringVector & stats) {
  GEOSAutoHandle geos;

  Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
  Rcpp::Function getValuesBlockFn = raster["getValuesBlock"];
  Rcpp::Function extentFn = raster["extent"];
  Rcpp::Function resFn = raster["res"];

  Rcpp::S4 extent = extentFn(rast);
  Rcpp::NumericVector res = resFn(rast);

  Grid<bounded_extent> grid {{
    extent.slot("xmin"),
    extent.slot("ymin"),
    extent.slot("xmax"),
    extent.slot("ymax"),
    },
    res[0],
    res[1]
  };

  bool store_values = false;
  for (const auto & stat : stats) {
    // explicit construction of std::string seems necessary to avoid ambiguous overload error
    if (stat == std::string("mode") ||
        stat == std::string("majority") ||
        stat == std::string("minority") ||
        stat == std::string("variety")) {
      store_values = true;
    }
  }

  RasterStats<double> raster_stats(store_values);
  Rcpp::NumericVector stat_results = Rcpp::no_init(stats.size());

  auto coverage_fraction = raster_cell_intersection(grid, geos.handle, read_wkb(geos.handle, wkb).get());
  auto& subgrid = coverage_fraction.grid();

  if (!subgrid.empty()) {
    Rcpp::NumericMatrix rast_values = getValuesBlockFn(rast,
                                                       1 + subgrid.row_offset(grid),
                                                       subgrid.rows(),
                                                       1 + subgrid.col_offset(grid),
                                                       subgrid.cols(),
                                                       "matrix");
    NumericMatrixRaster values(rast_values, subgrid);

    raster_stats.process(coverage_fraction, values);
  }

  int i = 0;
  for (const auto & stat : stats) {
    if (stat == std::string("mean")) stat_results[i] = raster_stats.mean();

    else if (stat == std::string("sum")) stat_results[i] = raster_stats.sum();
    else if (stat == std::string("count")) stat_results[i] = raster_stats.count();

    else if (stat == std::string("min")) stat_results[i] = raster_stats.min();
    else if (stat == std::string("max")) stat_results[i] = raster_stats.max();

    else if (stat == std::string("mode")) stat_results[i] = raster_stats.mode();
    else if (stat == std::string("majority")) stat_results[i] = raster_stats.mode();
    else if (stat == std::string("minority")) stat_results[i] = raster_stats.minority();

    else if (stat == std::string("variety")) stat_results[i] = raster_stats.variety();

    else Rcpp::stop("Unknown stat: " + stat);

    i++;
  }

  return stat_results;
}
