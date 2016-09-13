#!/bin/bash
cd ~/GitHub/sscr
rm -rfv man
rm -fv NAMESPACE
##R --slave -e "library(Rcpp); compileAttributes()"
R --slave -e "library(roxygen2); roxygenise('.')"
##rm -rfv ..Rcheck/ ..pdf
##rm -rfv src/*.o src/*.so src/*.rds
##rm -rfv src-i386/ src-x64/
R CMD build .
mkdir -p package-build
mv sscr_*.tar.gz package-build/
R CMD check package-build/sscr_*.tar.gz
R CMD INSTALL --install-tests package-build/sscr_*.tar.gz
