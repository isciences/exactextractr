if(getRversion() < "3.3.0") {
  stop("Your version of R is too old. This package requires R-3.3.0 or newer on Windows.")
}

VERSION <- commandArgs(TRUE)
if(!file.exists(sprintf("../windows/gdal2-%s/include/gdal/gdal.h", VERSION))){
  download.file(sprintf("https://github.com/rwinlib/gdal2/archive/v%s.zip", VERSION), "lib.zip", quiet = FALSE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}
