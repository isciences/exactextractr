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
#pragma once

#include <Rcpp.h>

#include "exactextract/src/raster.h"

// Construct a Raster using an R vector for storage
// This class uses row-major storage, consistent with the return value of
// raster::getValuesBlock, but inconsistent with the representation of
// matrices in R.
class NumericVectorRaster : public exactextract::AbstractRaster<double> {
public:
  NumericVectorRaster(const Rcpp::NumericVector & vec, const exactextract::Grid<exactextract::bounded_extent> & g) :
    AbstractRaster<double>(g),
    m_vec(vec)
  {}

  double operator()(size_t row, size_t col) const final {
    return m_vec[row*cols() + col];
  }

  const Rcpp::NumericVector vec() const {
    return m_vec;
  }

private:
  const Rcpp::NumericVector m_vec;
};

