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
        cat("Detection parameters: ", paste(format(round(det.pars, 2), nsmall = 2), collapse = ", "),
            "; Covariance parameters: ", paste(format(round(cov.pars, 2), nsmall = 2), collapse = ", "),
            "; nll: ", format(round(as.numeric(nll.obj$fn()), 2), nsmall = 2), "\n", sep = "")
    }
    as.numeric(nll.obj$fn())
}

cov.organise <- function(fit, survey.data, model.opts){
    pars <- fit$pars
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
    esa <- sum(det.probs)*mask.area
    D <- nrow(capt)/esa
    pars <- c(exp(pars), D = D, esa = esa)
    names(pars)[det.indices] <- paste("det_par_", 1:length(det.indices), sep = "")
    names(pars)[cov.indices] <- paste("cov_par_", 1:length(det.indices), sep = "")
    pars
}
