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

#pragma once

#include <memory>

#include <Rcpp.h>

#include <geos_c.h>

using geom_ptr = std::unique_ptr<GEOSGeometry, std::function<void(GEOSGeometry*)>>;
using wkb_reader_ptr = std::unique_ptr<GEOSWKBReader, std::function<void(GEOSWKBReader*)>>;

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
static inline geom_ptr read_wkb(const GEOSContextHandle_t & context, const Rcpp::RawVector & wkb) {
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
