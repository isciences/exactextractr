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

#include "raster_utils.h"

exactextract::Grid<exactextract::bounded_extent> make_grid(const Rcpp::S4 & rast) {
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

int get_nlayers(Rcpp::S4 & rast) {
  Rcpp::Environment raster = Rcpp::Environment::namespace_env("raster");
  Rcpp::Function nlayersFn = raster["nlayers"];

  Rcpp::NumericVector nlayersVec = nlayersFn(rast);

  return static_cast<int>(nlayersVec[0]);
}
