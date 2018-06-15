// Copyright (c) 2018 ISciences, LLC.
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

// [[Rcpp::plugins(cpp14)]]
#include <memory>

#include <Rcpp.h>
#include <geos_c.h>

#include "exactextract/src/extent.h"
#include "exactextract/src/raster_cell_intersection.h"
#include "exactextract/src/raster_stats.h"

#include "matrix_wrapper.h"

using geom_ptr= std::unique_ptr<GEOSGeometry, decltype(&GEOSGeom_destroy)>;
using wkb_reader_ptr = std::unique_ptr<GEOSWKBReader, decltype(&GEOSWKBReader_destroy)>;

using exactextract::Extent;
using exactextract::RasterCellIntersection;
using exactextract::RasterStats;

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

// Create an Extent from vectors representing the spatial extent and resolution
static Extent make_extent(const Rcpp::NumericVector & extent, const Rcpp::NumericVector & res) {
  double xmin = extent[0];
  double xmax = extent[1];
  double ymin = extent[2];
  double ymax = extent[3];

  double dx = res[0];
  double dy = res[1];

  return {xmin, ymin, xmax, ymax, dx, dy};
}

// Return a smart pointer to a Geometry, given WKB input
static geom_ptr read_wkb(const Rcpp::RawVector & wkb) {
  wkb_reader_ptr wkb_reader{ GEOSWKBReader_create(), GEOSWKBReader_destroy };

  geom_ptr geom{ GEOSWKBReader_read(wkb_reader.get(),
                                    std::addressof(wkb[0]),
                                    wkb.size()),
                 GEOSGeom_destroy };

  if (geom.get() == nullptr) {
    Rcpp::stop("Failed to parse WKB geometry");
  }

  return geom;
}

// Compute and return a single statistic, given a Matrix 'rast_values'
// that covers 'extent' at resolution 'res'.
// This is set up as a templated function so that NumericMatrix and
// IntegerMatrix can both be accepted. At the moment, we don't do
// anything intelligent in terms of return type, and all of the stats
// are returned as doubles, even when an integer would be more
// appropriate.
template<typename T>
double single_stat(const Rcpp::NumericVector & extent,
                   const Rcpp::NumericVector & res,
                   const T & rast_values,
                   const std::string & stat,
                   const Rcpp::RawVector & wkb)
{
  Extent ex = make_extent(extent, res);

  initGEOS(geos_warn, geos_error);

  const RasterCellIntersection rci(ex, read_wkb(wkb).get());

  auto mat = wrap(rast_values);
  RasterStats<decltype(mat)> stats{rci, mat, false};

  if (stat == "mean") return stats.mean();
  if (stat == "sum") return stats.sum();
  if (stat == "count") return stats.count();

  if (stat == "min") return stats.min();
  if (stat == "max") return stats.max();

  if (stat == "mode") return stats.mode();
  if (stat == "minority") return stats.minority();

  if (stat == "variety") return stats.variety();

  Rcpp::stop("Unknown stat: " + stat);
}

// [[Rcpp::export]]
Rcpp::List CPP_exact_extract(const Rcpp::NumericVector & extent,
                             const Rcpp::NumericVector & res,
                             const Rcpp::RawVector & wkb) {
  Extent ex = make_extent(extent, res);

  initGEOS(geos_warn, geos_error);

  const exactextract::RasterCellIntersection rci(ex, read_wkb(wkb).get());

  size_t nrow = rci.rows();
  size_t ncol = rci.cols();

  Rcpp::NumericMatrix weights = Rcpp::no_init(nrow, ncol);
  for (size_t i = 0; i < nrow; i++) {
    for (size_t j = 0; j < ncol; j++) {
      weights(i, j) = rci.get_local(i ,j);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("row")     = 1+rci.min_row(),
    Rcpp::Named("col")     = 1+rci.min_col(),
    Rcpp::Named("weights") = weights
  );
}

// [[Rcpp::export]]
double CPP_stat(const Rcpp::NumericVector & extent,
                const Rcpp::NumericVector & res,
                const Rcpp::NumericMatrix & rast_values,
                const std::string & stat,
                const Rcpp::RawVector & wkb)
{
  return single_stat(extent, res, rast_values, stat, wkb);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix CPP_weights(const Rcpp::NumericVector & extent,
                                const Rcpp::NumericVector & res,
                                const Rcpp::RawVector & wkb)
{
  Extent ex = make_extent(extent, res);

  size_t nrow = ex.rows();
  size_t ncol = ex.cols();

  initGEOS(geos_warn, geos_error);

  exactextract::RasterCellIntersection rci(ex, read_wkb(wkb).get());

  Rcpp::NumericMatrix weights(nrow, ncol);
  for (size_t i = 0; i < nrow; i++) {
    for (size_t j = 0; j < ncol; j++) {
      weights(i, j) = rci.get(i, j);
    }
  }

  return weights;
}
