#' Simulating SCR data with second-order spatial dependence
#'
#' Simulates SSCR data.
#'
#' @inheritParams fit.sscr
#' @param D Animal density (individuals per hectare).
#' @param det.pars List of detection function parameters.
#' @param cov.pars List of covariance parameters.
#'
#' @export
sim.sscr <- function(traps, mask, resp, cov.structure, D, det.pars, cov.pars){
    ## Finding extent of mask object.
    range.x <- range(mask[, 1])
    range.y <- range(mask[, 2])
    total.area <- diff(range.x)*diff(range.y)
    ## Simulating number of individuals.
    n <- rpois(1, D*total.area/10000)
    ## Simulating activity centre locations.
    locs.x <- runif(n, range.x[1], range.x[2])
    locs.y <- runif(n, range.y[1], range.y[2])
    locs <- cbind(locs.x, locs.y)
    ## Distances from activity centres to traps.
    ac.dists <- crossdist(locs[, 1], locs[, 2],
                          traps[, 1], traps[, 2])
    ## Distances between traps.
    trap.dists <- as.matrix(dist(traps))
    ## Extracting parameters.
    lambda0 <- det.pars$lambda0
    sigma <- det.pars$sigma
    ## Calculating "baseline" encounter rates.
    base.ers <- lambda0*exp(-ac.dists^2/(2*sigma^2))
    ## Extracting number of traps.
    n.traps <- nrow(traps)
    ## Constructing covariance matrix.
    if (cov.structure == "none"){
        cov <- matrix(0, nrow = n.traps, ncol = n.traps)
    } else if (cov.structure == "independent"){
        ## Extracting parameters.
        sigma.u <- cov.pars$sigma.u
        ## Specifying covariance.
        cov <- matrix(0, nrow = n.traps, ncol = n.traps)
        diag(cov) <- sigma.u
    } else if (cov.structure == "exponential"){
        ## Extracting parameters.
        sigma.u <- cov.pars$sigma.u
        rho <- cov.pars$rho
        ## Specifying covariance.
        cov <- sigma.u*exp(-trap.dists/rho)
    } else if (cov.structure == "matern"){
        stop("Matern covariance not yet implemented.")
    } else if (cov.structure == "constant"){
        sigma.u <- cov.pars$sigma.u
        cov <- matrix(sigma.u, nrow = n.traps, ncol = n.traps)
        stop("Full covariance not yet implemented.")
    }
    ## Simulating random effects.
    u.mat <- rmvnorm(n, rep(0, n.traps), cov)
    ## Generating capture histories.
    full.ers <- exp(log(base.ers) + u.mat)
    capt <- matrix(rpois(n*n.traps, full.ers), nrow = n)
    ## Removing undetected individuals.
    capt <- capt[apply(capt, 1, function(x) sum(x) > 0), ]
    if (resp == "binom"){
        capt[capt > 0] <- 1
    } else if (resp != "pois"){
        stop("Response type not recognised.")
    }
    capt
}
