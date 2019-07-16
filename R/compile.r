#' Compiling TMB C++ templates
#'
#' Compiles the \code{sscr} TMB templates into a shared object
#' file. This must be done a single time following installation or
#' updating of the package.
#'
#' @param template The TMB template(s) to compile. If \code{NULL}, all
#'     will be compiled, which is probably what you want to do.
#' 
#' @export
compile.sscr <- function(template = NULL){
    wd <- getwd()
    dir <- paste(system.file(package = "sscr"), "/tmb/src", sep = "")
    setwd(dir)
    if (!dir.exists("../bin")){
        dir.create("../bin")
    }
    files <- strsplit(list.files(), "[.]")
    base <- sapply(files, function(x) x[1])
    ext <- sapply(files, function(x) x[2])
    if (is.null(template)){
        v <- base[ext == "cpp"]
    } else {
        v <- template
    }
    for (i in v){
        compile(paste(i, ".cpp", sep = ""))
        unlink(paste(i, ".o", sep = ""))
        file.rename(paste(i, ".so", sep = ""),
                    paste("../bin/", i, ".so", sep = ""))
    }
    setwd(wd)
}

#' @import TMB Rcpp
#' @importFrom mvtnorm rmvnorm
#' @importFrom minqa bobyqa
#' @importFrom numDeriv grad
#' @importFrom spatstat crossdist
#' @importFrom stats dist dpois nlm nlminb optimHess plogis qlogis rbinom rnorm rpois runif sd
#' @useDynLib sscr
NULL

## Data documentation.

#' Test data
#'
#' Data used for unit tests.
#'
#' @name test.data
#' @format A list.
#' @usage test.data
#' @docType data
#' @keywords datasets
NULL
