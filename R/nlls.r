## Negative log-likelihood for independent random effects.
cov.nll <- function(pars, survey.data, model.opts){
    ## Extracting model options.
    resp.id <- model.opts$resp.id
    resp.pars <- model.opts$resp.pars
    detfn.id <- model.opts$detfn.id
    cov.id <- model.opts$cov.id
    det.indices <- model.opts$det.indices
    cov.indices <- model.opts$cov.indices
    detfn.scale.id <- model.opts$detfn.scale.id
    ## Extracting and unlinking parameters.
    link.ids <- model.opts$link.ids
    par.link <- link.closure(link.ids)
    par.unlink <- unlink.closure(link.ids)
    det.pars <- par.unlink(pars, det.indices)
    cov.pars <- par.unlink(pars, cov.indices)
    ## Extracting data.
    capt <- survey.data$capt
    mask.dists <- survey.data$mask.dists
    mask.area <- survey.data$mask.area
    n.mask <- survey.data$n.mask
    trap.dists <- survey.data$trap.dists
    n.traps <- survey.data$n.traps
    n <- nrow(capt)
    trace <- survey.data$trace
    ## Set up latent variables.
    if (cov.id == 3){
        u.detprob <- 0
        u.nll <- matrix(0, nrow = n, ncol = 1)
    } else {
        u.detprob <- numeric(n.traps)
        u.nll <- matrix(0, nrow = n, ncol = n.traps)
    }
    ## Calculating mask detection probabilities.
    det.probs <- numeric(n.mask)
    for (i in 1:n.mask){
        detprob.obj <- MakeADFun(data = list(mask_dists = mask.dists[i, ],
                                             trap_dists = trap.dists,
                                             n_traps = n.traps,
                                             detfn_id = detfn.id,
                                             detfn_scale_id = detfn.scale.id,
                                             resp_id = resp.id,
                                             resp_pars = resp.pars,
                                             cov_id = cov.id,
                                             det_pars = det.pars,
                                             cov_pars = cov.pars),
                                 parameters = list(u = u.detprob),
                                 random = "u", DLL = "cov_detprob",
                                 silent = TRUE)
        det.probs[i] <- exp(-detprob.obj$fn())
    }
    if (any(is.nan(det.probs))){
        warning("At least one detection probability is NaN.")
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
                                     resp_pars = resp.pars,
                                     detfn_id = detfn.id,
                                     detfn_scale_id = detfn.scale.id,
                                     cov_id = cov.id,
                                     det_probs = det.probs,
                                     det_pars = det.pars,
                                     cov_pars = cov.pars),
                         parameters = list(u = u.nll),
                         random = "u", DLL = "cov_nll", silent = TRUE)
    out <- as.numeric(nll.obj$fn())
    if (trace){
        cat("Detection parameters: ", paste(format(round(det.pars, 2), nsmall = 2), collapse = ", "),
            "; Covariance parameters: ", paste(format(round(cov.pars, 2), nsmall = 2), collapse = ", "),
            "; nll: ", format(round(out, 2), nsmall = 2), "\n", sep = "")
    }
    out
}

cov.organise <- function(fit, survey.data, model.opts){
    pars <- fit$par
    ## Extracting model options.
    resp.id <- model.opts$resp.id
    resp.pars <- model.opts$resp.pars
    detfn.id <- model.opts$detfn.id
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
    ## Set up latent variables.
    if (cov.id == 3){
        u.detprob <- 0
    } else {
        u.detprob <- numeric(n.traps)
    }
    for (i in 1:n.mask){
        detprob.obj <- MakeADFun(data = list(mask_dists = mask.dists[i, ],
                                             trap_dists = trap.dists,
                                             n_traps = n.traps,
                                             detfn_id = detfn.id,
                                             resp_id = resp.id,
                                             resp_pars = resp.pars,
                                             cov_id = cov.id,
                                             det_pars = det.pars,
                                             cov_pars = cov.pars),
                                 parameters = list(u = u.detprob),
                                 random = "u", DLL = "cov_detprob",
                                 silent = TRUE)
        det.probs[i] <- exp(-detprob.obj$fn())
    }
    esa <- sum(det.probs)*mask.area
    D <- nrow(capt)/esa
    pars <- c(exp(pars), D = D, esa = esa)
    names(pars)[det.indices] <- paste("det_par_", 1:length(det.indices), sep = "")
    names(pars)[cov.indices] <- paste("cov_par_", 1:length(cov.indices), sep = "")
    pars
}
