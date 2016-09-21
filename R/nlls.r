## Negative log-likelihood for independent random effects.
cov.nll <- function(pars, survey.data, model.opts){
    ## Extracting model options.
    resp.id <- model.opts$resp.id
    cov.id <- model.opts$cov.id
    det.indices <- model.opts$det.indices
    cov.indices <- model.opts$cov.indices
    ## Extracting and unlinking parameters.
    det.pars <- exp(pars[det.indices])
    cov.pars <- exp(pars[cov.indices])
    ## Extracting data.
    capt <- survey.data$capt
    mask.dists <- survey.data$mask.dists
    mask.area <- survey.data$mask.area
    n.mask <- survey.data$n.mask
    trap.dists <- survey.data$trap.dists
    n.traps <- survey.data$n.traps
    n <- nrow(capt)
    trace <- survey.data$trace
    ## Calculating mask detection probabilities.
    det.probs <- numeric(n.mask)
    for (i in 1:n.mask){
        detprob.obj <- MakeADFun(data = list(mask_dists = mask.dists[i, ],
                                             trap_dists = trap.dists,
                                             n_traps = n.traps,
                                             cov_id = cov.id,
                                             det_pars = det.pars,
                                             cov_pars = cov.pars),
                                 parameters = list(u = rep(0, n.traps)),
                                 random = "u", DLL = "cov_detprob",
                                 silent = TRUE)
        det.probs[i] <- exp(-detprob.obj$fn())
    }
    ## Calculating negative log-likelihood.
    nll.obj <- MakeADFun(data = list(capt = capt,
                                     mask_dists = mask.dists,
                                     trap_dists = trap.dists,
                                     n = n,
                                     n_traps = n.traps,
                                     n_mask = n.mask,
                                     mask_area = mask.area,
                                     resp_id = resp.id,
                                     cov_id = cov.id,
                                     det_probs = det.probs,
                                     det_pars = det.pars,
                                     cov_pars = cov.pars),
                         parameters = list(u = matrix(0, nrow = nrow(capt), ncol = n.traps)),
                         random = "u", DLL = "cov_nll", silent = TRUE)
    if (trace){
        cat("Detection parameters: ", paste(det.pars, collapse = ", "),
            "; Covariance parameters: ", paste(cov.pars, collapse = ", "),
            "; nll: ", as.numeric(nll.obj$fn()), "\n", sep = "")
    }
    as.numeric(nll.obj$fn())
}

