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
#pragma once

#include <memory>

#include <Rcpp.h>

#include "exactextract/src/grid.h"

#include "numeric_vector_raster.h"
#include "raster_utils.h"

// Read raster values from an R raster object
class S4RasterSource {
public:
  explicit S4RasterSource(SEXP rast, double default_value = std::numeric_limits<double>::quiet_NaN()) :
    m_grid(exactextract::Grid<exactextract::bounded_extent>::make_empty()),
    m_rast(rast),
    m_last_box(std::numeric_limits<double>::quiet_NaN(),
               std::numeric_limits<double>::quiet_NaN(),
               std::numeric_limits<double>::quiet_NaN(),
               std::numeric_limits<double>::quiet_NaN()),
    m_default_value(default_value)
  {
    m_grid = make_grid(rast);
  }

  S4RasterSource(SEXP rast,
                 const Rcpp::NumericVector & ext,
                 const Rcpp::NumericVector & res,
                 double default_value = std::numeric_limits<double>::quiet_NaN()) :
    m_grid(exactextract::Grid<exactextract::bounded_extent>::make_empty()),
    m_rast(rast),
    m_last_box(std::numeric_limits<double>::quiet_NaN(),
               std::numeric_limits<double>::quiet_NaN(),
               std::numeric_limits<double>::quiet_NaN(),
               std::numeric_limits<double>::quiet_NaN()),
    m_default_value(default_value)
  {
    m_grid = make_grid(ext, res);
  }

  const exactextract::Grid<exactextract::bounded_extent> &grid() const {
    return m_grid;
  }

  std::unique_ptr<exactextract::AbstractRaster<double>> read_box(const exactextract::Box & box, int layer) {
    auto cropped_grid = m_grid.crop(box);

    if (!(box == m_last_box)) {
      m_last_box = box;

      Rcpp::Environment xx = Rcpp::Environment::namespace_env("exactextractr");
      Rcpp::Function getValuesBlockFn = xx[".getValuesBlock"];

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
                                         cropped_grid.cols());

        if (!std::isnan(m_default_value)) {
          for (double& x : m_rast_values) {
            if (std::isnan(x)) {
              x = m_default_value;
            }
          }
        }
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
  exactextract::Grid<exactextract::bounded_extent> m_grid;
  SEXP m_rast;
  Rcpp::NumericMatrix m_rast_values;
  exactextract::Box m_last_box;
  double m_default_value;
};
