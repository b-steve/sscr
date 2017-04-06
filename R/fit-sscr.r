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
#' @param resp Response distribution for capture history
#'     elements. Either \code{"binom"} for a binomial distribution, or
#'     \code{"pois"} for a Poisson distribution.
#' @param resp.pars A named vector of known, fixed parameters for the
#'     response distribution. If \code{resp} is \code{"binom"}, then
#'     this must have a single element named \code{"size"} giving the
#'     fixed number of trials; if this argument is not provided, then
#'     the default is 1.
#' @param detfn Detection function, given by a character string. Use
#'     \code{"hn"} for halfnormal and \code{"hr"} for hazard rate.
#' @param cov.structure Covariance structure of the random
#'     effects. The current options are (1) \code{"none"} for no
#'     random effects (regular SCR), (2) \code{"independent"}, for
#'     independent random effects (equivalent to counts of detections
#'     being overdispersed), (3) \code{"exponential"}, for random
#'     effects with an exponential covariance structure, (4)
#'     "sq_exponential" for random effects with a squared exponential
#'     covariance structure, (5) \code{"matern"}, for random effects
#'     with a Matern covariance structure, (5) \code{"individual"},
#'     for random effects that are restricted to being the same at all
#'     traps (equivalent to having an independent random effect on
#'     \code{lambda0} for each individual), or (6)
#'     \code{"lc_exponential"} for a linear combination of exponential
#'     covariance functions.
#' @param start A named list of parameter start values.
#' @param trace Logical. If \code{TRUE}, parameter values for each
#'     step of the optimisation algorithm are printed.
#' @param test Logical. If \code{TRUE}, the likelihood is calculated
#'     at parameter start values. If \code{FALSE}, a model is fitted.
#' 
#' @export
fit.sscr <- function(capt, traps, mask, resp = "binom", resp.pars = NULL, detfn = "hn",
                     cov.structure = "none", start = NULL, trace = FALSE, test = FALSE){
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
    trap.dists <- as.matrix(dist(traps))
    ## Number of traps.
    n.traps <- nrow(traps)
    ## Packaging the data up into a list.
    survey.data <- list(capt = capt,
                        mask.dists = mask.dists,
                        mask.area = mask.area,
                        n.mask = n.mask,
                        trap.dists = trap.dists,
                        n.traps = n.traps,
                        trace = trace)
    model.opts <- list(resp = resp, resp.pars = resp.pars, detfn = detfn,
                       cov.structure = cov.structure, start = start)
    ## Optimisation object constructor function.
    if (cov.structure == "none"){
        make.obj <- make.obj.none
    } else {
        make.obj <- make.obj.cov
    }
    ## Making optimisation object.
    opt.obj <- make.obj(survey.data, model.opts)
    ## Fitting model or testing likelihood.
    if (test){
        fit <- opt.obj$fn(opt.obj$par)
    } else {
        raw.fit <- nlminb(opt.obj$par, opt.obj$fn, opt.obj$gr)
        if (cov.structure == "none"){
            fit <- summary(sdreport(opt.obj), "report")[, 1]
        } else {
            fit <- opt.obj$organise(raw.fit)
        }
    }
    fit
}

