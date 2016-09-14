#!/bin/bash
cd ~/GitHub/sscr
rm -rfv man
rm -rfv inst/tmb/bin
rm -fv NAMESPACE
R --slave -e "library(roxygen2); roxygenise('.')"
R CMD build .
mkdir -p package-build
mv sscr_*.tar.gz package-build/
R CMD check package-build/sscr_*.tar.gz
R CMD INSTALL --install-tests package-build/sscr_*.tar.gz
R --slave -e "library(sscr); compile.sscr()"
