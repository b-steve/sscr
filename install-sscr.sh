#!/bin/bash
cd ~/GitHub/sscr
rm -rfv man
rm -rfv inst/tmb/bin
rm -fv NAMESPACE
rm -fv src/*.o src/RcppExports.cpp R/RcppExports.R
rm -rfv package-build
R --slave -e "library(roxygen2); roxygenise('.')"
R --slave -e "library(Rcpp); compileAttributes()"
R CMD build .
mkdir -p package-build
mv sscr_*.tar.gz package-build/
R CMD check package-build/sscr_*.tar.gz --no-tests
R CMD INSTALL --install-tests package-build/sscr_*.tar.gz
R --slave -e "library(sscr); compile.sscr()"
