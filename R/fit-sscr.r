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
#'     elements. Either \code{"binom"} for a binomial distribution,
#'     \code{"pois"} for a Poisson distribution, or \code{"cmp"} for a
#'     Conway-Maxwell-Poisson distribution.
#' @param resp.pars A named vector of known, fixed parameters for the
#'     response distribution. If \code{resp} is \code{"binom"}, then
#'     this must have a single element named \code{"size"} giving the
#'     fixed number of trials; if this argument is not provided, then
#'     the default is 1.
#' @param detfn Detection function, given by a character string. Use
#'     \code{"hn"}, for halfnormal, \code{"hhn"} for hazard
#'     halfnormal, and \code{"hr"} for hazard rate.
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
#' @param re.multiplier A character string specifying how the expected
#'     encounter rate is calculated from the random effect. This is
#'     either (1) \code{"er"} for multiplication by the baseline
#'     encounter rate, or (2) \code{"prob"} for multiplication by the
#'     baseline encounter probability.
#' @param start A named vector of parameter start values.
#' @param fix A named vector of parameter values that are fixed in the
#'     model. Supercedes values provided in \code{start}.
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
#' 
#' @export
fit.sscr <- function(capt, traps, mask, resp = "binom",
                     resp.pars = NULL, detfn = "hhn",
                     cov.structure = "none", re.multiplier = "er",
                     start = NULL, fix = NULL, toa = NULL,
                     trace = FALSE, test = FALSE,
                     test.conditional.n = TRUE, hess = FALSE,
                     optim.fun = "nlminb"){
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
    ## Checking for identifiability problems.
    if (cov.structure != "none" & re.multiplier == "er" & detfn == "hhn"){
        if (is.null(fix)){
            fix <- numeric(0)
        }
        fix["mu.u"] <- 0
    }
    if (cov.structure != "none" & re.multiplier == "prob" & detfn == "hn"){
        if (is.null(fix)){
            fix <- numeric(0)
        }
        fix["g0"] <- 1
    }
    ## Sorting out fixed parameter values.
    if (is.null(start)){
        start <- list()
    }
    fix.names <- names(fix)
    ## Moving everything in fix to start.
    if (!is.null(fix)){
        for (i in 1:length(fix)){
            start[fix.names[i]] <- fix[i]
        }
    }
    start <- c(start, recursive = TRUE)
    ## Hard-coding manual separability to avoid having to maintain two
    ## templates.
    manual.sep <- TRUE
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
                       cov.structure = cov.structure, re.multiplier = re.multiplier,
                       start = start, fix.names = fix.names, conditional.n = TRUE,
                       Rhess = FALSE, manual.sep = manual.sep)
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
                                 cov.structure = cov.structure,
                                 re.multiplier = re.multiplier,
                                 start = start, fix.names = fix.names,
                                 conditional.n = test.conditional.n,
                                 Rhess = !test.conditional.n,
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
        ## Fitting the model with the appropriate optimiser.
        if (optim.fun == "nlminb"){
            raw.fit <- nlminb(opt.obj$par, opt.obj$fn, opt.obj$gr)
            ests.save <- opt.obj$organise(raw.fit$par, raw.fit$objective)
            fit <- list(ests = ests.save[[1]],
                        ests.link = ests.save[[2]],
                        grads = opt.obj$gr(raw.fit$par),
                        ll = -raw.fit$objective)
        } else if (optim.fun == "bobyqa"){
            raw.fit <- bobyqa(par = opt.obj$par, fn = opt.obj$fn)
            ests.save <- opt.obj$organise(raw.fit$par, raw.fit$fval)
            fit <- list(ests = ests.save[[1]],
                        ests.link = ests.save[[2]],
                        grads = opt.obj$gr(raw.fit$par),
                        bobyqa.code = raw.fit$ierr,
                        ll = -raw.fit$fval)
        } else if (optim.fun == "nlm"){
            f <- function(par){
                out <- opt.obj$fn(par)
                attr(out, "gradient") <- opt.obj$gr(par)
                out
            }
            raw.fit <- nlm(f, opt.obj$par)
            ests.save <- opt.obj$organise(raw.fit$estimate, raw.fit$minimum)
            fit <- list(ests = ests.save[[1]], ests.link = ests.save[[2]],
                        grads = raw.fit$gradient, ll = -raw.fit$minimum)
        }
        if (hess){
            if (trace){
                cat("Computing Hessian...\n")
            }
            model.opts.hess <- list(resp = resp, resp.pars = resp.pars, detfn = detfn,
                                    cov.structure = cov.structure,
                                    re.multiplier = re.multiplier,
                                    start = fit$ests, fix.names = fix.names,
                                    conditional.n = FALSE, Rhess = TRUE,
                                    manual.sep = manual.sep)
            opt.obj.hess <- make.obj(survey.data, model.opts.hess, any.cov)
            vcov.save <- opt.obj.hess$vcov(opt.obj.hess$par)
            fit$vcov <- vcov.save[[1]]
            fit$vcov.link <- vcov.save[[2]]
            fit$se <- sqrt(diag(fit$vcov))
            fit$se.link <- sqrt(diag(fit$vcov.link))
        }
    }
    fit$args <- list(traps = traps, mask = mask, resp = resp, resp.pars = resp.pars,
                     detfn = detfn, cov.structure = cov.structure,
                     re.multiplier = re.multiplier, det.names = opt.obj$det.names,
                     cov.names = opt.obj$cov.names)
    fit
}

