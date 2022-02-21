// Copyright (c) 2018-2021 ISciences, LLC.
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

#include "geos_r.h"

#include "raster_utils.h"

#include "exactextract/src/raster.h"
#include "exactextract/src/raster_cell_intersection.h"

using exactextract::RasterView;
using exactextract::raster_cell_intersection;

// [[Rcpp::export]]
Rcpp::S4 CPP_coverage_fraction(Rcpp::S4 & rast,
                               const Rcpp::RawVector & wkb,
                               bool crop)
{
  try {
    GEOSAutoHandle geos;
    Rcpp::Environment xx = Rcpp::Environment::namespace_env("exactextractr");

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

    if (rast.inherits("SpatRaster")) {
      Rcpp::Environment terra = Rcpp::Environment::namespace_env("terra");
      Rcpp::Function rastFn = terra["rast"];
      Rcpp::Function extFn = terra["ext"];
      Rcpp::Function crsFn = terra["crs"];

      Rcpp::S4 ext = extFn(grid.xmin(), grid.xmax(), grid.ymin(), grid.ymax());

      return rastFn(weights,
                    Rcpp::Named("ext") = ext,
                    Rcpp::Named("crs") = crsFn(rast));
    } else {
      Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
      Rcpp::Function rasterFn = raster["raster"];
      Rcpp::Function crsFn = xx[".crs"];

      return rasterFn(weights,
                      Rcpp::Named("xmn")=grid.xmin(),
                      Rcpp::Named("xmx")=grid.xmax(),
                      Rcpp::Named("ymn")=grid.ymin(),
                      Rcpp::Named("ymx")=grid.ymax(),
                      Rcpp::Named("crs")=crsFn(rast));
    }
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
