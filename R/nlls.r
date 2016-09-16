## Negative log-likelihood for independent random effects.
ind.nll <- function(pars, survey.data){
    ## Extracting parameters.
    log.lambda0 <- pars[1]
    log.sigma <- pars[2]
    log.sigma.u <- pars[3]
    ## Extracting data.
    capt <- survey.data$capt
    mask.dists <- survey.data$mask.dists
    mask.area <- survey.data$mask.area
    n.mask <- survey.data$n.mask
    n.traps <- survey.data$n.traps
    n <- nrow(capt)
    trace <- survey.data$trace
    ## Calculating mask detection probabilities.
    det.probs <- numeric(n.mask)
    for (i in 1:n.mask){
        detprob.obj <- MakeADFun(data = list(dists = mask.dists[i, ],
                                             n_traps = n.traps,
                                             lambda0 = exp(log.lambda0),
                                             sigma = exp(log.sigma),
                                             sigma_u = exp(log.sigma.u)),
                                 parameters = list(u = rep(0, n.traps)),
                                 random = "u", DLL = "ind_detprob",
                                 silent = TRUE)
        det.probs[i] <- exp(-detprob.obj$fn())
    }
    ## Calculating negative log-likelihood.
    nll.obj <- MakeADFun(data = list(capt = capt,
                                     mask_dists = mask.dists,
                                     n = n,
                                     n_traps = n.traps,
                                     n_mask = n.mask,
                                     mask_area = mask.area,
                                     det_probs = det.probs,
                                     log_lambda0 = log.lambda0,
                                     log_sigma = log.sigma,
                                     log_sigma_u = log.sigma.u),
                         parameters = list(u = matrix(0, nrow = nrow(capt), ncol = n.traps)),
                         random = "u", DLL = "ind_nll", silent = TRUE)
    if (trace){
        cat("lambda0: ", exp(log.lambda0), ", sigma:", exp(log.sigma), ", sigma.u:", exp(log.sigma.u),
            ", nll: ", as.numeric(nll.obj$fn()), "\n", sep = "")
    }
    as.numeric(nll.obj$fn())
}

