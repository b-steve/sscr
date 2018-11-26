#' Fitting an SCR model with second-order spatial dependence
#'
#' Fits an SSCR model. Estimation is by maximum likelihood. The
#' second-order spatial dependence is modelled via trap-level random
#' effects for each detected individual. The likelihood function is
#' calculated by integrating over these random effects using the
#' Laplace approximation.
#'
#' @param capt A capture history object. It should be a matrix, where
#'     the jth element of the ith row should provide a detection
#'     record of the ith individual at the jth detector.
#' @param traps A matrix with two columns, providing the Cartesian
#'     coordinates of the detector locations.
#' @param mask A mask object for integration over the survey area.
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
#' @param detfn.scale A character string, either \code{"er"} or
#'     \code{"prob"}. This indicates whether the detection function
#'     should provide the encounter rate (expected number of
#'     detections) or the probability of detection.
#' @param cov.structure Covariance structure of the random
#'     effects. The current options are (1) \code{"none"} for no
#'     random effects (regular SCR), (2) \code{"independent"}, for
#'     independent random effects (equivalent to counts of detections
#'     being overdispersed), (3) \code{"exponential"}, for random
#'     effects with an exponential covariance structure, (4)
#'     \code{"sq_exponential"} for random effects with a squared
#'     exponential covariance structure, (5) \code{"matern"}, for
#'     random effects with a Matern covariance structure, (5)
#'     \code{"individual"}, for random effects that are restricted to
#'     being the same at all traps (equivalent to having an
#'     independent random effect on \code{lambda0} for each
#'     individual), or (6) \code{"lc_exponential"} for a linear
#'     combination of exponential covariance functions.
#' @param re.scale A character string, either \code{"er"} or
#'     \code{"prob"}. This indicates whether the Gaussian random
#'     effects affect the encounter rate (expected number of
#'     detections) or the probability of detection.
#' @param start A named list of parameter start values.
#' @param toa A matrix with the same dimensions as \code{capt} that
#'     provides time-of-arrival information for acoustic detections.
#' @param trace Logical. If \code{TRUE}, parameter values for each
#'     step of the optimisation algorithm are printed.
#' @param test Logical. If \code{TRUE}, the negative log-likelihood is
#'     calculated at parameter start values. If \code{FALSE}, a model
#'     is fitted. Alternatively, a character string. If \code{"nll"},
#'     then the negative log-likelihood is calculated. If \code{"gr"},
#'     then the partial derivatives of the negative log-likelihood
#'     function with respect to the parameters is also calculated. If
#'     \code{"hess"} then the Hessian if also calculated.
#' @param test.conditional.n Logical. If \code{TRUE}, tests are
#'     carried out conditioning in the number of detections.
#' @param hess Logical. If \code{TRUE}, a Hessian is computed and a
#'     variance-covariance matrix is returned.
#' @param optim.fun A character string representing the R function to
#'     maximise the likelihood. This can be \code{"bobyqa"} (from
#'     package \code{minqa}), \code{"nlminb"}, or \code{"nlm"}.
#' @param manual.sep Logical. If \code{TRUE}, integration over
#'     individuals' spatially structured random effects is enforced in
#'     R, rather than automatically detected by TMB.
#' 
#' @export
fit.sscr <- function(capt, traps, mask, resp = "binom", resp.pars = NULL, detfn = "hn",
                     detfn.scale = "er", cov.structure = "none", re.scale = "er",
                     start = NULL, toa = NULL, trace = FALSE, test = FALSE,
                     test.conditional.n = TRUE, hess = FALSE, optim.fun = "nlminb",
                     manual.sep = TRUE){
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
                        toa = toa,
                        trace = trace)
    model.opts <- list(resp = resp, resp.pars = resp.pars, detfn = detfn,
                       detfn.scale = detfn.scale, cov.structure = cov.structure,
                       re.scale = re.scale, start = start, conditional.n = TRUE, Rhess = FALSE,
                       manual.sep = manual.sep)
    ## Optimisation object constructor function.
    if (cov.structure == "none"){
        any.cov <- FALSE
    } else {
        any.cov <- TRUE
    }
    opt.obj <- make.obj(survey.data, model.opts, any.cov)
    ## Fitting model or testing likelihood.
    if (is.logical(test)){
        if (test){
            test <- "nll"
        }
    }
    if (test %in% c("nll", "gr", "hess")){
        model.opts.test <-  list(resp = resp, resp.pars = resp.pars, detfn = detfn,
                                 detfn.scale = detfn.scale, cov.structure = cov.structure,
                                 re.scale = re.scale, start = start,
                                 conditional.n = test.conditional.n, Rhess = !test.conditional.n,
                                 manual.sep = manual.sep)
        opt.obj.test <- make.obj(survey.data, model.opts.test, any.cov)
        ## Setting up output list.
        fit <- list()
        ## Computing negative log-likelihood.
        fit$nll <- opt.obj.test$fn(opt.obj.test$par)
        if (test %in% c("gr", "hess")){
            ## Creating numerical gradient function for non-exact methods.
            if (is.null(opt.obj.test$gr)){
                opt.obj.test$gr <- function(x){
                    message("Determining test gradients numerically...")
                    grad(opt.obj.test$fn, x)
                }
            }
            ## Calculating 
            fit$gr <- opt.obj.test$gr(opt.obj.test$par)
            if (test %in% "hess"){
                fit$vcov <- opt.obj.test$vcov(opt.obj.test$par)
               
            }
        }
    } else {
        if (optim.fun == "nlminb"){
            raw.fit <- nlminb(opt.obj$par, opt.obj$fn, opt.obj$gr)
            fit <- list(ests = opt.obj$organise(raw.fit$par, raw.fit$objective),
                        grads = opt.obj$gr(raw.fit$par))
        } else if (optim.fun == "bobyqa"){
            raw.fit <- bobyqa(par = opt.obj$par, fn = opt.obj$fn)
            fit <- list(ests = opt.obj$organise(raw.fit$par, raw.fit$fval),
                        grads = opt.obj$gr(raw.fit$par),
                        bobyqa.code = raw.fit$ierr)
        } else if (optim.fun == "nlm"){
            f <- function(par){
                out <- opt.obj$fn(par)
                attr(out, "gradient") <- opt.obj$gr(par)
                out
            }
            raw.fit <- nlm(f, opt.obj$par)
            fit <- list(ests = opt.obj$organise(raw.fit$estimate, raw.fit$minimum),
                        grads = raw.fit$gradient)
        }
        if (hess){
            if (trace){
                cat("Computing Hessian...\n")
            } 
            model.opts.hess <-  list(resp = resp, resp.pars = resp.pars, detfn = detfn,
                                     detfn.scale = detfn.scale, cov.structure = cov.structure,
                                     re.scale = re.scale, start = fit$ests,
                                     conditional.n = FALSE, Rhess = TRUE, manual.sep = manual.sep)
            opt.obj.hess <- make.obj(survey.data, model.opts.hess, any.cov)
            fit$vcov <- opt.obj.hess$vcov(opt.obj.hess$par)
            fit$se <- sqrt(diag(fit$vcov))
        }
        
    }
    fit
}

