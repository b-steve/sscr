#' Fitting an SCR model with second-order spatial dependence
#'
#' Fits an SSCR model. Estimation is by maximum likelihood. The
#' second-order spatial dependence is modelled via trap-level random
#' effects for each detected individual. The likelihood function is
#' calculated by integrating over these random effects using the
#' Laplace approximation.
#'
#' @param capt Capture history object.
#' @param traps Traps object.
#' @param mask Mask object.
#' @param resp Response distribution for capture history elements.
#' @param re.structure Covariance structure of the random effects. The
#'     current options are (1) \code{"none"} for no random effects
#'     (regular SCR), (2) \code{"independent"}, for independent random
#'     effects (equivalent to counts of detections being
#'     overdispersed), (3) \code{"constant"}, for independent random
#'     effects that are restricted to being the same at all traps
#'     (equivalent to having a random effect on \code{lambda0} for
#'     each individual), and (4) \code{"exponential"}, for random
#'     effects with an exponential covariance structure (spatially
#'     structured random effects).
#' @param trace Logical. If \code{TRUE}, parameter values for each
#'     step of the optimisation algorithm are printed.
#' @param test Logical. If \code{TRUE}, the likelihood is calculated
#'     at parameter start values. If \code{FALSE}, a model is fitted.
#' 
#' @export
fit.sscr <- function(capt, traps, mask, resp = "binom", re.structure = "none", trace = FALSE, test = FALSE){
    ## Loading DLLs.
    dll.dir <- paste(system.file(package = "sscr"), "/tmb/bin/", sep = "")
    for (i in paste(dll.dir, list.files(dll.dir), sep = "")){
        dyn.load(i)
    }
    ## Calculating mask distances.
    mask.dists <- crossdist(mask[, 1], mask[, 2],
                            traps[, 1], traps[, 2])
    ## Extracting mask area.
    mask.area <- attr(mask, "area")
    ## Number of mask points.
    n.mask <- nrow(mask)
    ## Distances between traps.
    
    ## Number of traps.
    n.traps <- nrow(traps)
    ## Packaging the data up into a list.
    survey.data <- list(capt = capt,
                        mask.dists = mask.dists,
                        mask.area = mask.area,
                        n.mask = n.mask,
                        n.traps = n.traps,
                        trace = trace)
    ## Optimisation object constructor function.
    make.obj <- switch(re.structure,
                       none = make.obj.none,
                       independent = make.obj.independent,
                       make.obj.error)
    ## Making optimisation object.
    opt.obj <- make.obj(survey.data, resp)
    ## Fitting model or testing likelihood.
    if (test){
        fit <- opt.obj$fn(opt.obj$par)
    } else {
        fit <- nlminb(opt.obj$par, opt.obj$fn, opt.obj$gr)
    }
    fit
}

