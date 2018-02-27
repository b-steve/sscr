## Negative log-likelihood for independent random effects.
cov.nll <- function(pars, survey.data, model.opts, only.detprobs = FALSE){
    ## Extracting model options.
    resp.id <- model.opts$resp.id
    resp.pars <- model.opts$resp.pars
    detfn.id <- model.opts$detfn.id
    cov.id <- model.opts$cov.id
    det.indices <- model.opts$det.indices
    cov.indices <- model.opts$cov.indices
    D.indices <- model.opts$D.indices
    detfn.scale.id <- model.opts$detfn.scale.id
    re.scale.id <- model.opts$re.scale.id
    toa.id <- model.opts$toa.id
    conditional.n <- as.numeric(model.opts$conditional.n)
    ## Extracting and unlinking parameters.
    link.ids <- model.opts$link.ids
    par.link <- link.closure(link.ids)
    par.unlink <- unlink.closure(link.ids)
    det.pars <- par.unlink(pars, det.indices)
    cov.pars <- par.unlink(pars, cov.indices)
    if (conditional.n){
        D <- 1
    } else {
        D <- par.unlink(pars, D.indices)
    }
    sigma.toa <- 0
    if (toa.id == 1){
        toa.indices <- model.opts$toa.indices
        sigma.toa <- par.unlink(pars, toa.indices)
    }
    ## Extracting data.
    capt <- survey.data$capt
    n.dets <- survey.data$n.dets
    mask.dists <- survey.data$mask.dists
    mask.area <- survey.data$mask.area
    n.mask <- survey.data$n.mask
    trap.dists <- survey.data$trap.dists
    n.traps <- survey.data$n.traps
    toa.ssq <- survey.data$toa.ssq
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
                                             re_scale_id = re.scale.id,
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
    if (only.detprobs){
        out <- det.probs
    } else {
        ## Calculating negative log-likelihood.
        nll.obj <- MakeADFun(data = list(capt = capt,
                                         n_dets = n.dets,
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
                                         re_scale_id = re.scale.id,
                                         det_probs = det.probs,
                                         toa_id = toa.id,
                                         toa_ssq = toa.ssq,
                                         det_pars = det.pars,
                                         cov_pars = cov.pars,
                                         sigma_toa = sigma.toa/1000,
                                         conditional_n = conditional.n,
                                         D = D),
                             parameters = list(u = u.nll),
                             random = "u", DLL = "cov_nll", silent = TRUE)
        out <- as.numeric(nll.obj$fn())
        if (trace){
            cat("Detection parameters: ", paste(format(round(det.pars, 2), nsmall = 2), collapse = ", "),
                "; Covariance parameters: ", paste(format(round(cov.pars, 2), nsmall = 2), collapse = ", "),
                "; D: ", format(round(D, 2), nsmall = 2), sep = "")
            if (toa.id == 1){
                cat("; TOA parameter: ", format(round(sigma.toa, 2), nsmall = 2), sep = "")
            }
            cat("; nll: ", format(round(out, 2), nsmall = 2), "\n", sep = "")
        }
    }
    out
}

cov.organise <- function(fit, survey.data, model.opts, det.probs.fun){
    pars <- fit$par
    det.probs <- det.probs.fun(pars)
    mask.area <- survey.data$mask.area
    capt <- survey.data$capt
    esa <- sum(det.probs)*mask.area
    D <- nrow(capt)/esa
    link.ids <- model.opts$link.ids
    par.unlink <- unlink.closure(link.ids)
    ll <- -fit$objective
    pars <- c(par.unlink(pars), D, esa, ll)
    names(pars) <- c(model.opts$par.names, "D", "esa", "LL")
    pars    
}
