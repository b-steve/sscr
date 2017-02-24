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
sim.sscr <- function(traps, mask, resp, resp.pars, cov.structure, D, det.pars = NULL, cov.pars = NULL){
    ## Finding extent of mask object.
    range.x <- range(mask[, 1])
    range.y <- range(mask[, 2])
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
        ## Distances between traps.
        trap.dists <- as.matrix(dist(traps))
        ## Extracting parameters.
        lambda0 <- det.pars$lambda0
        sigma <- det.pars$sigma
        ## Calculating "baseline" encounter rates.
        base.ers <- lambda0*exp(-ac.dists^2/(2*sigma^2))
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
        } else if (cov.structure == "individual"){
            sigma.u <- cov.pars$sigma.u
            cov <- matrix(sigma.u, nrow = n.traps, ncol = n.traps)
            stop("Full covariance not yet implemented.")
        } else if (cov.structure == "lc_exponential"){
            stop("Linear combination of exponentials not yet implemented.")
        } else if (cov.structure == "sq_exponential"){
            ## Extracting parameters.
            sigma.u <- cov.pars$sigma.u
            rho <- cov.pars$rho
            ## Specifying covariance.
            cov <- sigma.u*exp(-(trap.dists^2)/(rho^2))
        }
        ## Simulating random effects.
        u.mat <- rmvnorm(n, rep(0, n.traps), cov)
        ## Generating capture histories.
        full.ers <- exp(log(base.ers) + u.mat)
        if (resp == "pois"){
            capt <- matrix(rpois(n*n.traps, full.ers), nrow = n)
        } else if (resp == "binom"){
            capt <- matrix(rbinom(n*n.traps, resp.pars, 1 - exp(-full.ers)), nrow = n)
        }
    }
    ## Removing undetected individuals.
    capt <- capt[apply(capt, 1, function(x) sum(x) > 0), ]
    capt
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

## ac <- c(150, 250)
## locs <- sim.ou(ac, 80, 75, 1000)
## traps <- test.data$traps
## mask <- test.data$mask
## epsilon <- 0.5*sqrt(10000*attr(mask, "area"))
## plot(mask, type = "n")
## points(locs, type = "o")
## points(ac[1], ac[2], col = "red", pch = 16)
## points(traps, col = "blue", pch = 16)

## count.dets(locs, traps, 0.7*epsilon)