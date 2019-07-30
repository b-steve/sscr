#' Simulating SCR data with second-order spatial dependence
#'
#' Simulates SSCR data.
#'
#' @param fit A fitted model object returned by \link{fit.sscr}(). If
#'     provided, simulation is based on this fitted model, and all
#'     other arguments are ignored.
#' @inheritParams fit.sscr
#' @param D Animal density (individuals per hectare).
#' @param det.pars List of detection function parameters.
#' @param cov.pars List of covariance parameters.
#' @param toa.pars List of time-of-arrival parameters.
#' @param trunc.acs Logical. If \code{TRUE}, activity centres are
#'     truncated to the region covered by the mask.
#' @param keep.acs Logical. If \code{TRUE}, activity centre locations
#'     are returned along with the captures histories.
#' @param mask.spacing Optional. The distance between adjacent mask
#'     points. Providing a value will save processing time if activity
#'     centres are being truncated to the mask region.
#'
#' @export
sim.sscr <- function(fit = NULL, traps, mask, D, resp = NULL,
                     resp.pars = NULL, detfn = "hhn",
                     cov.structure = "none", re.multiplier = "er",
                     det.pars = NULL, cov.pars = NULL,
                     toa.pars = NULL, trunc.acs = TRUE, keep.acs = FALSE,
                     mask.spacing = NULL){
    ## Overwriting arguments from a model object, if provided.
    if (!missing(fit)){
        traps <- fit$args$traps
        mask <- fit$args$mask
        D <- fit$ests["D"]
        resp <- fit$args$resp
        resp.pars <- fit$args$resp.pars
        detfn <- fit$args$detfn
        cov.structure <- fit$args$cov.structure
        re.multiplier <- fit$args$re.multiplier
        det.names <- fit$args$det.names
        det.pars <- as.vector(fit$ests[det.names], mode = "list")
        if (cov.structure != "none"){
            cov.names <- fit$args$cov.names
            cov.pars <- as.vector(fit$ests[cov.names], mode = "list")
        } else {
            cov.pars <- NULL
        }
        if (any(names(fit$ests) == "sigma.toa")){
            toa.pars <- fit$ests["sigma.toa"]
        } else {
            toa.pars <- NULL
        }
    }
    if (is.null(resp.pars)){
        resp.pars <- 1
    }
    sim.toa <- !is.null(toa.pars)
    if (!is.null(resp)){
        if (sim.toa & (resp != "binom")){
            if (resp.pars$size != 1){
                stop("Time of arrival can only be simulated if the response distribution is Bernoulli.")
            }
        }
    }
    ## Sorting out identifiability.
    if (!(cov.structure %in% c("none", "OU")) & re.multiplier == "er" & detfn == "hhn"){
        if (is.null(cov.pars$mu.u)){
            cov.pars$mu.u <- 0
        } else {
            if (cov.pars$mu.u != 0){
                warning("The mu.u and lambda0 parameters are not identifiable for models with a hazard halfnormal detection function and an encounter-rate random effect multiplier. Setting nonzero mu.u is not recommended.")
            }
        }
        if (is.null(det.pars$lambda0)){
            stop("A value for lambda0 must be supplied.")
        }
    }
    if (!(cov.structure %in% c("none", "OU")) & re.multiplier == "prob" & detfn == "hn"){
        if (is.null(det.pars$g0)){
            cov.pars$g0 <- 1
        } else {
            if (det.pars$g0 != 1){
                warning("The mu.u and g0 parameters are not identifiable for models with a halfnormal detection function and a probability random effect multiplier. Setting g0 to anything but 1 is not recommended.")
            }
        }
        if (is.null(cov.pars$mu.u)){
            stop("A value of mu.u must be supplied.")
        }
    }
    ## Calculating mask spacing if not provided.
    if (is.null(mask.spacing)){
        mask.dists <- crossdist(mask[, 1], mask[, 2],
                                mask[, 1], mask[, 2])
        mask.spacing <- min(mask.dists[mask.dists > 0])
    }
    ## Finding extent of mask object.
    range.x <- range(mask[, 1]) + c(-mask.spacing, mask.spacing)/2
    range.y <- range(mask[, 2]) + c(-mask.spacing, mask.spacing)/2
    total.area <- diff(range.x)*diff(range.y)
    ## Extracting number of traps.
    n.traps <- nrow(traps)
    ## Simulating number of individuals.
    n <- rpois(1, D*total.area/10000)
    ## Simulating activity centre locations.
    acs.x <- runif(n, range.x[1], range.x[2])
    acs.y <- runif(n, range.y[1], range.y[2])
    acs <- cbind(acs.x, acs.y)
    ## Distances from activity centres to traps.
    ac.dists <- crossdist(acs[, 1], acs[, 2],
                          traps[, 1], traps[, 2])
    ## Truncating activity centres beyond the mask.
    if (trunc.acs){
        ## Function to check for mask-point proximity.
        mask.proximity <- function(ac){
            any(abs(ac[1] - mask[, 1]) < mask.spacing/2 &
                abs(ac[2] - mask[, 2]) < mask.spacing/2)
        }
        keep.ac <- apply(acs, 1, mask.proximity)
        ## Keeping only activity centres within a mask point's pixel.
        acs.trunc <- acs[!keep.ac, ]
        acs <- acs[keep.ac, ]
        ac.dists <- ac.dists[keep.ac, ]
        n <- nrow(acs)
    } else {
        acs.trunc <- NULL
    }
    if (cov.structure == "OU"){
        ## Extracting movement parameters.
        tau <- cov.pars$tau
        sigma <- cov.pars$sigma
        n.steps <- cov.pars$n.steps
        epsilon <- cov.pars$epsilon
        capt <- matrix(0, nrow = n, ncol = n.traps)
        for (i in 1:n){
            locs <- sim.ou(acs[i, ], tau, sigma, n.steps)
            capt[i, ] <- count.dets(locs, traps, epsilon)
        }
    } else {
        if (is.null(resp)){
            stop("A response distribution must be specified.")
        }
        ## Sorting out detection function.
        calc.detfn <- detfn.closure(detfn, det.pars)
        ## Distances between traps.
        trap.dists <- as.matrix(dist(traps))
        ## Calculating "baseline" encounter rates and probabilities.
        base.prob <- calc.detfn(ac.dists, er = FALSE)
        base.er <- calc.detfn(ac.dists, er = TRUE)
        ## Constructing covariance matrix.
        if (cov.structure == "none"){
            mu.u <- 0
            cov <- matrix(0, nrow = n.traps, ncol = n.traps)
        } else if (cov.structure == "independent"){
            ## Extracting parameters.
            mu.u <- cov.pars$mu.u
            sigma.u <- cov.pars$sigma.u
            ## Specifying covariance.
            cov <- matrix(0, nrow = n.traps, ncol = n.traps)
            diag(cov) <- sigma.u
        } else if (cov.structure == "exponential"){
            ## Extracting parameters.
            mu.u <- cov.pars$mu.u
            sigma.u <- cov.pars$sigma.u
            rho <- cov.pars$rho
            ## Specifying covariance.
            cov <- sigma.u^2*exp(-trap.dists/rho)
        } else if (cov.structure == "matern"){
            stop("Matern covariance not yet implemented.")
        } else if (cov.structure == "individual"){
            mu.u <- cov.pars$mu.u
            sigma.u <- cov.pars$sigma.u
            cov <- matrix(sigma.u^2, nrow = n.traps, ncol = n.traps)
        } else if (cov.structure == "lc_exponential"){
            stop("Linear combination of exponentials not yet implemented.")
        } else if (cov.structure == "sq_exponential"){
            ## Extracting parameters.
            mu.u <- cov.pars$mu.u
            sigma.u <- cov.pars$sigma.u
            rho <- cov.pars$rho
            ## Specifying covariance.
            cov <- sigma.u^2*exp(-(trap.dists^2)/(rho^2))
        }
        ## Simulating random effects.
        u.mat <- rmvnorm(n, rep(mu.u, n.traps), cov)
        ## Getting full encounter rates and probabilities.
        if (re.multiplier == "er"){
            full.er <- base.er*exp(u.mat)
        } else if (re.multiplier == "prob"){
            full.er <- base.prob*exp(u.mat)
        }
        full.prob <- 1 - exp(-full.er)
        ## Generating capture histories.
        if (resp == "pois"){
            capt <- matrix(rpois(n*n.traps, full.er), nrow = n)
        } else if (resp == "binom"){
            capt <- matrix(rbinom(n*n.traps, resp.pars, full.prob), nrow = n)
        }
        ## Generating times of arrival.
        if (sim.toa){
            n.dets <- sum(capt > 0)
            toa <- matrix(0, nrow = n, ncol = n.traps)
            toa[capt > 0] <- ac.dists[capt > 0]/toa.pars$sound.speed +
                rnorm(n.dets, 0, toa.pars$sigma.toa/1000)
        }
    }
    ## Removing undetected individuals.
    if (sim.toa){
        acs.missed <- acs[!apply(capt, 1, function(x) sum(x) > 0), ]
        acs <- acs[apply(capt, 1, function(x) sum(x) > 0), ]
        toa <- toa[apply(capt, 1, function(x) sum(x) > 0), ]
        capt <- capt[apply(capt, 1, function(x) sum(x) > 0), ]
        if (keep.acs){
            out <- list(capt = capt, toa = toa, acs = acs,
                        acs.trunc = acs.trunc, acs.missed = acs.missed)
        } else {
            out <- list(capt = capt, toa = toa)
        }
    } else {
        acs.missed <- acs[!apply(capt, 1, function(x) sum(x) > 0), ]
        acs <- acs[apply(capt, 1, function(x) sum(x) > 0), ]
        capt <- capt[apply(capt, 1, function(x) sum(x) > 0), ]
        if (keep.acs){
            out <- list(capt = capt, acs = acs, acs.trunc = acs.trunc,
                        acs.missed = acs.missed)
        } else {
            out <- capt
        }
    }
    out
}

#' Simulating an Ornstein-Uhlenbeck process
#'
#' Simulates animal movement from an Ornstein-Uhlenbeck process.
#'
#' @param mu Point of attraction (conceptually the home-range centre).
#' @param tau Temporal autocorrelation parameter.
#' @param sigma Spatial autocorrelation parameter.
#' @param start Coordinates of starting location. Default is to
#'     generate a location from the long-term stationary distribution.
#' @param n.steps Total number of time steps.
#' 
#' @author This is a modified version of a function provided by Theo
#'     Michelot.
#' 
#' @export
sim.ou <- function(mu, tau, sigma, n.steps, start = NULL){
    if (is.null(start)){
        start <- rmvnorm(1, mu, sigma^2*diag(2))
    }
    ## Changing Theo's parameterisation.
    b <- -1/tau
    v <- sigma^2
    out <- matrix(0, nrow = n.steps, ncol = 2)
    out[1, ] <- start
    for (i in 2:n.steps){
        out[i, ] <- rmvnorm(1, mu + exp(b)*(out[i - 1, ] - mu),
                            v*(1 - exp(2*b))*diag(2))
                            
    }
    out
}

## Function to count the number of detections by each trap.
count.dets <- function(locs, traps, epsilon){
    loc.dists <- crossdist(locs[, 1], locs[, 2],
                           traps[, 1], traps[, 2])
    apply(loc.dists, 2, function(x, epsilon) sum(x <= epsilon), epsilon = epsilon)
}

## Closure for detection function.
detfn.closure <- function(detfn, pars){
    if (detfn == "hn"){
        g0 <- pars$g0
        sigma <- pars$sigma
        out <- function(d, er = FALSE){
            out <- g0*exp(-d^2/(2*sigma^2))
            if (er){
                out <- -log(1 - out)
            }
            out
        }
    } else if (detfn == "hr"){
        g0 <- pars$g0
        sigma <- pars$sigma
        z <- pars$z
        out <- function(d, er = FALSE){
            out <- g0*(1 - exp(-((d/sigma)^-z)))
            if (er){
                out <- -log(1 - out)
            }
        }
    } else if (detfn == "hhn"){
        lambda0 <- pars$lambda0
        sigma <- pars$sigma
        out <- function(d, er = FALSE){
            out <- lambda0*exp(-d^2/(2*sigma^2))
            if (!er){
                out <- 1 - exp(-out)
            }
            out
        } 
    } else if (detfn == "hhr"){
        lambda0 <- pars$lambda0
        sigma <- pars$sigma
        z <- pars$z
        out <- function(d, er = FALSE){
            out <- -lambda0*(1 - exp(-((d/sigma)^-z)))
            if (!er){
                out <- 1 - exp(-out)
            }
            out
        }
    }
    out
}
