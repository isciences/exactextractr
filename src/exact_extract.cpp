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
  for (size_t i = rci.min_row(); i < rci.max_row(); i++) {
    for (size_t j = rci.min_col(); j < rci.max_col(); j++) {
      weights(i, j) = rci.get(i, j);
    }
  }

  return weights;
}

// [[Rcpp::export]]
SEXP CPP_stats(Rcpp::S4 & rast, const Rcpp::RawVector & wkb, const Rcpp::StringVector & stats) {
  initGEOS(geos_warn, geos_error);

  Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
  Rcpp::Function getValuesBlockFn = raster["getValuesBlock"];
  Rcpp::Function extentFn = raster["extent"];
  Rcpp::Function resFn = raster["res"];

  Rcpp::S4 extent = extentFn(rast);
  Rcpp::NumericVector res = resFn(rast);

  Extent ex {
    extent.slot("xmin"),
    extent.slot("ymin"),
    extent.slot("xmax"),
    extent.slot("ymax"),
    res[0],
    res[1]
  };

  exactextract::RasterCellIntersection rci(ex, read_wkb(wkb).get());

  Rcpp::NumericVector stat_results = Rcpp::no_init(stats.size());

  Rcpp::NumericMatrix rast_values = getValuesBlockFn(rast,
                                                     1 + rci.min_row(),
                                                     rci.rows(),
                                                     1 + rci.min_col(),
                                                     rci.cols(),
                                                     "matrix");

  auto mat = wrap(rast_values);
  RasterStats<decltype(mat)> raster_stats{rci, mat, true};

  int i = 0;
  for (const auto & stat : stats) {
    if (stat == std::string("mean")) stat_results[i] = raster_stats.mean();

    else if (stat == std::string("sum")) stat_results[i] = raster_stats.sum();
    else if (stat == std::string("count")) stat_results[i] = raster_stats.count();

    else if (stat == std::string("min")) stat_results[i] = raster_stats.min();
    else if (stat == std::string("max")) stat_results[i] = raster_stats.max();

    else if (stat == std::string("mode")) stat_results[i] = raster_stats.mode();
    else if (stat == std::string("minority")) stat_results[i] = raster_stats.minority();

    else if (stat == std::string("variety")) stat_results[i] = raster_stats.variety();

    else Rcpp::stop("Unknown stat: " + stat);

    i++;
  }

  return stat_results;
}
