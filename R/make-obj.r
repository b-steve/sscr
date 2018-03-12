## Functions to make optimisation objects.
make.obj <- function(survey.data, model.opts, any.cov){
    ## Extracting data.
    capt <- survey.data$capt
    n.dets <- apply(capt, 1, function(x) sum(x > 0))
    mask.dists <- survey.data$mask.dists
    mask.area <- survey.data$mask.area
    n.mask <- survey.data$n.mask
    n.traps <- survey.data$n.traps
    trap.dists <- survey.data$trap.dists
    n <- nrow(capt)
    ## Extracting response type.
    resp <- model.opts$resp
    resp.pars <- model.opts$resp.pars
    resp.id <- switch(model.opts$resp, binom = 0, pois = 1)
    if (is.null(resp.pars)){
        resp.pars <- 1
    }
    ## Extracting detection function and scale.
    detfn <- model.opts$detfn
    detfn.id <- switch(detfn, hn = 0, hr = 1)
    detfn.scale <- model.opts$detfn.scale
    detfn.scale.id <- switch(detfn.scale, er = 0, prob = 1)
    ## Sorting out map for fixing parameters.
    map <- list()
    ## Extracting toa and sorting out indicator.
    toa <- survey.data$toa
    if (is.null(toa)){
        toa.ssq <- as.matrix(0)
        toa.id <- 0
        map[["link_sigma_toa"]] <- factor(NA)
    } else {
        ## Hard-coding speed of sound here.
        sound.speed <- 330
        toa.ssq <- make_toa_ssq(toa, t(mask.dists), sound.speed)
        toa.id <- 1
    }
    ## Extracting trace.
    trace <- survey.data$trace
    ## Indicator for conditional likelihood.
    conditional.n <- model.opts$conditional.n
    Rhess <- model.opts$Rhess
    ## Start values for optimisation.
    start <- model.opts$start
    start.names <- names(start)
    ## Indices and start values for detection function parameters.
    ## Link 0 is log, 1 is qlogis.
    ## For detection functions on the hazard scale.
    if (detfn.scale.id == 0){
        if (detfn.id == 0){
            ## .. With a halfnormal detection function.
            det.indices <- 1:2
            det.start <- numeric(2)
            det.start[1] <- ifelse(any(start.names == "lambda0"), start["lambda0"],
                                   max(capt)/(2*resp.pars))
            det.start[2] <- ifelse(any(start.names == "sigma"), start["sigma"],
                                   max(apply(mask.dists, 1, min))/5)
            det.link.ids <- c(0, 0)
            par.names <- c("lambda0", "sigma")
        } else if (detfn.id == 1){
            ## .. With a hazard rate detection function.
            det.indices <- 1:3
            det.start <- numeric(3)
            det.start[1] <- ifelse(any(start.names == "lambda0") , start["lambda0"],
                                   max(capt)/(2*resp.pars))
            det.start[2] <- ifelse(any(start.names == "sigma"), start["sigma"],
                                   max(apply(mask.dists, 1, min))/5)
            det.start[3] <- ifelse(any(start.names == "z") , start["z"], 1)                   
            det.link.ids <- c(0, 0, 0)
            par.names <- c("lambda0", "sigma", "z")
        }
        ## For detection functions on the probability scale.
    } else if (detfn.scale.id == 1){
        if (detfn.id == 0){
            ## .. With a halfnormal detection function.
            det.indices <- 1:2
            det.start <- numeric(2)
            det.start[1] <- ifelse(any(start.names == "g0") , start["g0"], 0.5)
            det.start[2] <- ifelse(any(start.names == "sigma"), start["sigma"],
                                   max(apply(mask.dists, 1, min))/5)
            det.link.ids <- c(1, 0)
            par.names <- c("g0", "sigma")
        } else if (detfn.id == 1){
            ## .. With a hazard rate detection function.
            det.indices <- 1:3
            det.start <- numeric(3)
            det.start[1] <- ifelse(any(start.names == "g0") , start["g0"], 0.5)
            det.start[2] <- ifelse(any(start.names == "sigma"), start["sigma"],
                                   max(apply(mask.dists, 1, min))/5)
            det.start[3] <- ifelse(any(start.names == "z") , start["z"], 1)           
            det.link.ids <- c(1, 0, 0)
            par.names <- c("g0", "sigma", "z")
        }
    }
    pars.start <- det.start
    link.ids <- det.link.ids
    #if (any.cov){
        ## Stuff only for covariance structures.
        cov.structure <- model.opts$cov.structure
        cov.id <- switch(model.opts$cov.structure,
                         independent = 0,
                         exponential = 1,
                         matern = 2,
                         individual = 3,
                         lc_exponential = 4,
                         sq_exponential = 5,
                         none = 6)
        re.scale <- model.opts$re.scale
        re.scale.id <- switch(re.scale, er = 0, prob = 1)     
        ## Indices and start values for covariance parameters.
        cov.index.start <- max(det.indices) + 1
        if (cov.id == 0){
            ## Independent.
            cov.indices <- cov.index.start
            cov.start <- numeric(1)
            cov.start[1] <- ifelse(any(start.names == "sigma.u"),
                                   start["sigma.u"], sd(capt))
            cov.link.ids <- 0
            par.names <- c(par.names, "sigma.u")
        } else if (cov.id == 1){
            ## Exponential.
            cov.indices <- cov.index.start:(cov.index.start + 1)
            cov.start <- numeric(2)
            cov.start[1] <- ifelse(any(start.names == "sigma.u"),
                                   start["sigma.u"], sd(capt))
            cov.start[2] <- ifelse(any(start.names == "rho"),
                                   start["rho"], mean(trap.dists))
            cov.link.ids <- c(0, 0)
            par.names <- c(par.names, "sigma.u", "rho")
        } else if (cov.id == 2){
            ## Matern.
            cov.indices <- cov.index.start:(cov.index.start + 2)
            stop("Matern covariance function not yet implemented.")
        } else if (cov.id == 3){
            ## Individual.
            cov.indices <- cov.index.start
            cov.start <- numeric(1)
            cov.start[1] <- ifelse(any(start.names == "sigma.u"),
                                   start["sigma.u"], sd(capt))
            cov.link.ids <- 0
            par.names <- c(par.names, "sigma.u")
        } else if (cov.id == 4){
            ## Linear combination of exponentials.
            cov.indices <- cov.index.start:(cov.index.start + 3)
            stop("Linear combination of exponentials not yet implemented.")
        } else if (cov.id == 5){
            ## Squared exponential.
            cov.indices <- cov.index.start:(cov.index.start + 1)
            cov.start <- numeric(2)
            cov.start[1] <- ifelse(any(start.names == "sigma.u"),
                                   start["sigma.u"], sd(capt))
            cov.start[2] <- ifelse(any(start.names == "rho"),
                                   start["rho"], mean(trap.dists))
            cov.link.ids <- c(0, 0)
            par.names <- c(par.names, "sigma.u", "rho")
        } else if (cov.id == 6){
            ## No covariance.
            cov.indices <- -1
            cov.start <- NULL
            cov.link.ids <- NULL
            map[["link_cov_pars"]] <- factor(NA)
        }
        pars.start <- c(pars.start, cov.start)
        link.ids <- c(link.ids, cov.link.ids)
                                        #}
    ## Start value for TOA parameter.
    toa.indices <- -1
    if (toa.id == 1){
        toa.indices <- length(pars.start) + 1
        toa.start <- ifelse(any(start.names == "sigma.toa"), start["sigma.toa"], 3)
        pars.start <- c(pars.start, toa.start)
        link.ids <- c(link.ids, 0)
        par.names <- c(par.names, "sigma.toa")
    }
    ## Start value for density parameter.
    D.indices <- -1
    if (!conditional.n){
        D.indices <- length(pars.start) + 1
        D.start <- ifelse(any(start.names == "D"),
                          start["D"], 5)
        pars.start <- c(pars.start, D.start)
        link.ids <- c(link.ids, 0)
        par.names <- c(par.names, "D")
        if (Rhess){
            map[["link_D"]] <- factor(NA)
        }
    } else {
        map[["link_D"]] <- factor(NA)
    }
    model.opts$link.ids <- link.ids
    model.opts$par.names <- par.names
    ## Getting par.link and par.unlink.
    par.link <- link.closure(link.ids)
    par.unlink <- unlink.closure(link.ids)
    par.dlink <- dlink.closure(link.ids)
    ## Converting parameters to link scale.
    link.pars.start <- par.link(pars.start)
    ##if (any.cov){
    detprob.objs <- list()
    ## Making detprob AD objects.
    if (cov.id == 6){
        u.detprob <- 0
        u.nll <- matrix(0, nrow = 1, ncol = 1)
        random.comp <- NULL
        map[["u"]] <- factor(NA)
    } else if (cov.id == 3){
        u.detprob <- 0
        u.nll <- matrix(0, nrow = n, ncol = 1)
        random.comp <- "u"
    } else {
        u.detprob <- numeric(n.traps)
        u.nll <- matrix(0, nrow = n, ncol = n.traps)
        random.comp <- "u"
    }
    for (i in 1:n.mask){
        detprob.objs[[i]] <- MakeADFun(data = list(mask_dists = mask.dists[i, ],
                                                   trap_dists = trap.dists,
                                                   n_traps = n.traps,
                                                   detfn_id = detfn.id,
                                                   detfn_scale_id = detfn.scale.id,
                                                   resp_id = resp.id,
                                                   resp_pars = resp.pars,
                                                   cov_id = cov.id,
                                                   re_scale_id = re.scale.id,
                                                   link_det_ids = link.ids[det.indices],
                                                   link_cov_ids = if (cov.id == 6) 0 else link.ids[cov.indices]),
                                       parameters = list(link_det_pars = link.pars.start[det.indices],
                                                         link_cov_pars = if (cov.id == 6) 1 else link.pars.start[cov.indices],
                                                         link_sigma_toa = ifelse(toa.id, link.pars.start[toa.indices], 1),
                                                         link_D = ifelse(conditional.n | Rhess, 1, link.pars.start[D.indices]),
                                                         u = u.detprob),
                                       map = map, random = random.comp, DLL = "cov_detprob", silent = TRUE)
    }
    get.fn.gr <- function(fun = "nll"){
        function(link.pars){
            if (Rhess){
                link.pars.tmb <- link.pars[-D.indices]
            } else {
                link.pars.tmb <- link.pars
            }
            det.probs <- numeric(n.mask)
            neglog.det.probs <- numeric(n.mask)
            neglog.det.probs.grads <- matrix(0, nrow = length(link.pars.tmb), ncol = n.mask)
            for (i in 1:n.mask){
                neglog.det.probs[i] <- detprob.objs[[i]]$fn(link.pars.tmb)
                det.probs[i] <- exp(-neglog.det.probs[i])
                neglog.det.probs.grads[, i] <- detprob.objs[[i]]$gr(link.pars.tmb)
            }
            if (fun != "det.probs"){
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
                                                 conditional_n = as.numeric(conditional.n | Rhess),
                                                 link_det_ids = link.ids[det.indices],
                                                 link_cov_ids = if (cov.id == 6) 0 else link.ids[cov.indices]),
                                     parameters = list(link_det_pars = link.pars[det.indices],
                                                       link_cov_pars = if (cov.id == 6) 1 else link.pars[cov.indices],
                                                       link_sigma_toa = ifelse(toa.id, link.pars[toa.indices], 1),
                                                       link_D = ifelse(conditional.n | Rhess, 1, link.pars[D.indices]),
                                                       u = u.nll),
                                     map = map, random = random.comp, DLL = "cov_nll", silent = TRUE)
            }
            if (fun == "nll"){
                out <- nll.obj$fn(link.pars.tmb)
                ## Log-likelihood component due to number of
                ## detected individuals, if calculated in R.
                if (Rhess){
                    out <- out - dpois(n, exp(link.pars[D.indices])*mask.area*
                                          sum(det.probs), TRUE)
                } else if (conditional.n){
                    ## Making the condtional-n likelihood look
                    ## like the non-conditional-n likelihood for
                    ## consistency (the conditional-n likelihood
                    ## is only used so that maximisation happens
                    ## over one fewer parameter).
                    out <- out - dpois(n, n, TRUE)
                }
                if (trace){
                    cat("Detection parameters: ", paste(format(round(par.unlink(link.pars, det.indices), 2), nsmall = 2),
                                                        collapse = ", "),
                        "; Covariance parameters: "[cov.id != 6], paste(format(round(par.unlink(link.pars, cov.indices), 2), nsmall = 2),
                                                                        collapse = ", ")[cov.id != 6])
                    if (!conditional.n){
                        cat("; D: ", format(round(par.unlink(link.pars, D.indices), 2), nsmall = 2), sep = "")
                    }
                    if (toa.id == 1){
                        cat("; TOA parameter: ", format(round(par.unlink(link.pars, toa.indices), 2), nsmall = 2), sep = "")
                    }
                    cat("; NLL: ", format(round(out, 2), nsmall = 2), "\n", sep = "")
                }
            } else if (fun == "gr"){
                if (Rhess){
                    out <- numeric(length(link.pars))
                    out[-D.indices] <- nll.obj$gr(link.pars.tmb) + exp(link.pars[D.indices])*mask.area*apply(neglog.det.probs.grads, 1, function(x) sum(-exp(-neglog.det.probs)*x))
                    out[D.indices] <- mask.area*sum(det.probs)*exp(link.pars[D.indices]) - n
                } else if (!conditional.n){
                    out <- nll.obj$gr(link.pars.tmb) + exp(link.pars[D.indices])*mask.area*apply(neglog.det.probs.grads, 1, function(x) sum(-exp(-neglog.det.probs)*x))
                } else {
                    out <- nll.obj$gr(link.pars.tmb) + n*apply(neglog.det.probs.grads, 1, function(x) sum(-exp(-neglog.det.probs)*x))/sum(det.probs)    
                }
                if (trace){
                    cat("Partial derivatives: ", paste(out, collapse = " "), "\n")
                }
            } else if (fun == "det.probs"){
                out <- det.probs
            }
            out
        }
    }
    obj.fn <- get.fn.gr(fun = "nll")
    obj.gr <- get.fn.gr(fun = "gr")
    obj.det.probs <- 
        obj.vcov <- vcov.closure(survey.data, model.opts, obj.fn, par.dlink, gr = obj.gr)
    obj.organise <- organise.closure(survey.data, model.opts, cov.organise, get.fn.gr(fun = "det.probs"))
    obj <- list(par = link.pars.start, fn = obj.fn, gr = obj.gr,
                vcov = obj.vcov, organise = obj.organise)
    obj
}

make.obj.error <- function(survey.data, resp){
    stop("Random effects structure not recognised. Use \"none\" or \"independent.\"")
}

## Closure to provide negative log-likelihood function without passing data.
nll.closure <- function(survey.data, model.opts, nll.fun){
    function(pars){
        nll.fun(pars, survey.data, model.opts)
    }
}

## Closure to provide negative log-likelihood gradient function without passing data.
grad.closure <- function(survey.data, model.opts, grad.fun){
    function(pars){
        grad.fun(pars, survey.data, model.opts)
    }
}

## Closure to provide variance-covariance calculation, given parameter estimates.
vcov.closure <- function(survey.data, model.opts, nll, dlink.fun, gr = NULL){
    function(pars){
        hess.link <- optimHess(pars, nll, gr = gr)
        vcov.link <- solve(hess.link)
        n.pars <- length(pars)
        jacobian <- diag(n.pars)
        diag(jacobian) <- dlink.fun(pars)
        out <- jacobian %*% vcov.link %*% t(jacobian)
        dimnames(out) <- list(model.opts$par.names, model.opts$par.names)
        out
    }
}

## Closure to provide organisation function stuff after model fitting.
organise.closure <- function(survey.data, model.opts, organise.fun, det.probs.fun = NULL){
    function(pars, objective){
        organise.fun(pars, objective, survey.data, model.opts, det.probs.fun)
    }
}

cov.organise <- function(pars, objective, survey.data, model.opts, det.probs.fun){
    det.probs <- det.probs.fun(pars)
    mask.area <- survey.data$mask.area
    capt <- survey.data$capt
    esa <- sum(det.probs)*mask.area
    D <- nrow(capt)/esa
    link.ids <- model.opts$link.ids
    par.unlink <- unlink.closure(link.ids)
    ll <- -objective
    pars <- c(par.unlink(pars), D, esa, ll)
    names(pars) <- c(model.opts$par.names, "D", "esa", "LL")
    pars    
}


## Closure to provide linking function without passing link ids.
link.closure <- function(link.ids){
    function(pars, which = NULL){
        n.pars <- length(pars)
        if (is.null(which)){
            which <- 1:n.pars
        }
        out <- numeric(n.pars)
        for (i in (1:n.pars)[which]){
            out[i] <- links[[link.ids[i] + 1]](pars[i])
        }
        out[which]
    }
}

unlink.closure <- function(link.ids){
    function(link.pars, which = NULL){
        n.pars <- length(link.pars)
        if (is.null(which)){
            which <- 1:n.pars
        }
        out <- numeric(n.pars)
        for (i in (1:n.pars)[which]){
            out[i] <- unlinks[[link.ids[i] + 1]](link.pars[i])
        }
        out[which]
    }
}

dlink.closure <- function(link.ids){
    function(link.pars, which = NULL){
        n.pars <- length(link.pars)
        if (is.null(which)){
            which <- 1:n.pars
        }
        out <- numeric(n.pars)
        for (i in (1:n.pars)[which]){
            out[i] <- dlinks[[link.ids[i] + 1]](link.pars[i])
        }
        out[which]
    }
}

links <- list(log, qlogis)
unlinks <- list(exp, plogis)
dlinks <- list(exp, dlogis)
