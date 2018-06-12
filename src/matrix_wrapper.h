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

// This file provides a wrapper around the Rcpp NumericMatrix and
// IntegerMatrix types, so that their elements can be accessed by
// code expecting a type of the form Container<T>

template<typename T>
struct MatrixWrapper {};

template<>
struct MatrixWrapper<double> {
  MatrixWrapper(const Rcpp::NumericMatrix& mat) : m_mat{mat} {}

  using value_type= double;

  double operator()(size_t i, size_t j) const {
    return m_mat(i, j);
  }

  const Rcpp::NumericMatrix& m_mat;
};

template<>
struct MatrixWrapper<int> {
  MatrixWrapper(const Rcpp::IntegerMatrix& mat) : m_mat{mat} {}

  using value_type= int;

  double operator()(size_t i, size_t j) const {
    return m_mat(i, j);
  }

  const Rcpp::IntegerMatrix& m_mat;
};

MatrixWrapper<double> wrap(const Rcpp::NumericMatrix & m) {
  return { m };
}

MatrixWrapper<int> wrap(const Rcpp::IntegerMatrix & m) {
  return { m };
}

