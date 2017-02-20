#' Compiling TMB C++ templates
#'
#' Compiles the \code{sscr} TMB templates into a shared object
#' file. This must be done a single time following installation or
#' updating of the package.
#'
#' @export
compile.sscr <- function(){
    wd <- getwd()
    dir <- paste(system.file(package = "sscr"), "/tmb/src", sep = "")
    setwd(dir)
    if (!dir.exists("../bin")){
        dir.create("../bin")
    }
    files <- strsplit(list.files(), "[.]")
    base <- sapply(files, function(x) x[1])
    ext <- sapply(files, function(x) x[2]) 
    for (i in base[ext == "cpp"]){
        compile(paste(i, ".cpp", sep = ""))
        unlink(paste(i, ".o", sep = ""))
        file.rename(paste(i, ".so", sep = ""),
                    paste("../bin/", i, ".so", sep = ""))
    }
    setwd(wd)
}

#' @import TMB
#' @importFrom mvtnorm rmvnorm
#' @importFrom spatstat crossdist
#' @importFrom stats dist nlminb rbinom rpois runif sd
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
