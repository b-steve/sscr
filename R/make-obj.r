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
    resp.id <- switch(model.opts$resp, binom = 0, pois = 1, cmp = 2, nb = 3, nba = 4)
    if (resp.id == 1){
        resp.pars <- 1
    }
    if (is.null(resp.pars)){
        if (resp.id == 0 | resp.id == 2){
            resp.pars <- 1
        } else if (resp.id == 3){
            resp.pars  <- 10
        } else if (resp.id == 4){
            resp.pars <- 2
        }
    }
    ## Extracting detection function and scale.
    detfn <- model.opts$detfn
    detfn.id <- switch(detfn, hn = 0, hr = 1, hhn = 2, hhr = 3)
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
    ## Fixed parameters.
    fix.names <- model.opts$fix.names
    ## Indicator for manual separability.
    manual.sep <- model.opts$manual.sep
    ## Uhh I dunno this stuff only works under certain combinations of
    ## manual.sep, conditional.n, and Rhess for confusing reasons that
    ## I keep forgetting, but I momentarily understand right now. I
    ## don't think this error should ever get triggered by standard
    ## use of fit.sscr(). If it does, then uh oh something went wrong.
    if (manual.sep){
        if (!(conditional.n | Rhess)){
            stop("For manual separability, either 'conditional.n' or 'Rhess' must be 'TRUE'.")
        }
    }
    ## Indices and start values for detection function parameters.
    ## Link 0 is log, 1 is qlogis, 2 is identity, 3 is log(x - 1).
    ## For detection functions on the hazard scale.
    if (detfn.id == 0){
        ## ... For a halfnormal detection function.
        det.indices <- 1:2
        det.start <- numeric(2)
        det.start[1] <- ifelse(any(start.names == "g0") , start["g0"], 0.5)
        det.start[2] <- ifelse(any(start.names == "sigma"), start["sigma"],
                               max(apply(mask.dists, 1, min))/5)
        det.link.ids <- c(1, 0)
        det.names <- c("g0", "sigma")
    } else if (detfn.id == 1){
        ## ... For a hazard rate detection function.
        det.indices <- 1:3
        det.start <- numeric(3)
        det.start[1] <- ifelse(any(start.names == "g0") , start["g0"], 0.5)
        det.start[2] <- ifelse(any(start.names == "sigma"), start["sigma"],
                               max(apply(mask.dists, 1, min))/5)
        det.start[3] <- ifelse(any(start.names == "z") , start["z"], 1)           
        det.link.ids <- c(1, 0, 0)
        det.names <- c("g0", "sigma", "z")
    } else if (detfn.id == 2){
        ## ... For a hazard halfnormal detection function.
        det.indices <- 1:2
        det.start <- numeric(2)
        det.start[1] <- ifelse(any(start.names == "lambda0"), start["lambda0"],
                        ifelse(resp.id == 0, max(capt)/(2*resp.pars), max(capt)/2))
        det.start[2] <- ifelse(any(start.names == "sigma"), start["sigma"],
                               max(apply(mask.dists, 1, min))/5)
        det.link.ids <- c(0, 0)
        det.names <- c("lambda0", "sigma")
    } else if (detfn.id == 3){
        ## ... For a hazard hazard rate (lol) detection function.
        det.indices <- 1:3
        det.start <- numeric(3)
        det.start[1] <- ifelse(any(start.names == "lambda0") , start["lambda0"],
                        ifelse(resp.id == 0, max(capt)/(2*resp.pars), max(capt)))
        det.start[2] <- ifelse(any(start.names == "sigma"), start["sigma"],
                               max(apply(mask.dists, 1, min))/5)
        det.start[3] <- ifelse(any(start.names == "z") , start["z"], 1)                   
        det.link.ids <- c(0, 0, 0)
        det.names <- c("lambda0", "sigma", "z")
    }
    par.names <- det.names
    pars.start <- det.start
    link.ids <- det.link.ids
    ## Setting map for detection function parameters.
    det.map <- factor(seq_along(det.start))
    det.map[det.names %in% fix.names] <- NA
    map[["link_det_pars"]] <- det.map
    ## Stuff only for covariance structures.
    cov.id <- switch(model.opts$cov.structure,
                     independent = 0,
                     exponential = 1,
                     matern = 2,
                     individual = 3,
                     lc_exponential = 4,
                     sq_exponential = 5,
                     none = 6)
    re.multiplier <- model.opts$re.multiplier
    mult.id <- switch(model.opts$re.multiplier,
                      er = 0, prob = 1)
    ## Indices and start values for covariance parameters.
    cov.index.start <- max(det.indices) + 1
    if (cov.id == 0){
        ## Independent.
        cov.indices <- cov.index.start:(cov.index.start + 1)
        cov.start <- numeric(2)
        cov.start[1] <- ifelse(any(start.names == "mu.u"),
                               start["mu.u"], 0)
        cov.start[2] <- ifelse(any(start.names == "sigma.u"),
                               start["sigma.u"], sd(capt))
        cov.link.ids <- c(2, 0)
        cov.names <- c("mu.u", "sigma.u")
    } else if (cov.id == 1){
        ## Exponential.
        cov.indices <- cov.index.start:(cov.index.start + 2)
        cov.start <- numeric(3)
        cov.start[1] <- ifelse(any(start.names == "mu.u"),
                               start["mu.u"], 0)
        cov.start[2] <- ifelse(any(start.names == "sigma.u"),
                               start["sigma.u"], sd(capt))
        cov.start[3] <- ifelse(any(start.names == "rho"),
                               start["rho"], min(trap.dists[trap.dists > 0]))
        cov.link.ids <- c(2, 0, 0)
        cov.names <- c("mu.u", "sigma.u", "rho")
    } else if (cov.id == 2){
        ## Matern.
        cov.indices <- cov.index.start:(cov.index.start + 3)
        cov.start <- numeric(4)
        cov.start[1] <- ifelse(any(start.names == "mu.u"),
                               start["mu.u"], 0)
        cov.start[2] <- ifelse(any(start.names == "sigma.u"),
                               start["sigma.u"], sd(capt))
        cov.start[3] <- ifelse(any(start.names == "phi"),
                               start["phi"], min(trap.dists[trap.dists > 0]))
        cov.start[4] <- ifelse(any(start.names == "kappa"),
                               start["kappa"], 1)
        cov.link.ids <- c(2, 0, 0, 0)
        cov.names <- c("mu.u", "sigma.u", "phi", "kappa")
    } else if (cov.id == 3){
        ## Individual.
        cov.indices <- cov.index.start:(cov.index.start + 1)
        cov.start <- numeric(2)
        cov.start[1] <- ifelse(any(start.names == "mu.u"),
                               start["mu.u"], 0)
        cov.start[2] <- ifelse(any(start.names == "sigma.u"),
                               start["sigma.u"], sd(capt))
        cov.link.ids <- c(2, 0)
        cov.names <- c("mu.u", "sigma.u")
    } else if (cov.id == 4){
        ## Linear combination of exponentials.
        cov.indices <- cov.index.start:(cov.index.start + 4)
        stop("Linear combination of exponentials not yet implemented.")
    } else if (cov.id == 5){
        ## Squared exponential.
        cov.indices <- cov.index.start:(cov.index.start + 2)
        cov.start <- numeric(3)
        cov.start[1] <- ifelse(any(start.names == "mu.u"),
                               start["mu.u"], 0)
        cov.start[2] <- ifelse(any(start.names == "sigma.u"),
                               start["sigma.u"], sd(capt))
        cov.start[3] <- ifelse(any(start.names == "rho"),
                               start["rho"], min(trap.dists[trap.dists > 0]))
        cov.link.ids <- c(2, 0, 0)
        cov.names <- c("mu.u", "sigma.u", "rho")
    } else if (cov.id == 6){
        ## No covariance.
        cov.indices <- -1
        cov.start <- NULL
        cov.link.ids <- NULL
        cov.names <- NULL
        map[["link_cov_pars"]] <- factor(NA)
    }
    ## Setting map for covariance parameters.
    if (cov.id != 6){
        cov.map <- factor(seq_along(cov.start))
        cov.map[cov.names %in% fix.names] <- NA
        map[["link_cov_pars"]] <- cov.map
    }
    par.names <- c(par.names, cov.names)
    pars.start <- c(pars.start, cov.start)
    link.ids <- c(link.ids, cov.link.ids)
    ## Indices and start values for response distribution parameters.
    if (resp.id == 0 | resp.id == 1){
        resp.indices <- length(pars.start) + 1
        resp.names <- "N"
        resp.start <- resp.pars
        resp.link.ids <- 2
        fix.names <- c(fix.names, "N")
        resp.map <- factor(NA)
    } else if (resp.id == 2){
        resp.indices <- length(pars.start) + 1
        resp.names <- "nu"
        resp.start <- resp.pars
        resp.link.ids <- 0
        resp.map <- factor(seq_along(resp.start))
    } else if (resp.id == 3){
        resp.indices <- length(pars.start) + 1
        resp.names <- "size"
        resp.start <- resp.pars
        resp.link.ids <- 0
        resp.map <- factor(seq_along(resp.start))
    } else if (resp.id == 4){
        resp.indices <- length(pars.start) + 1
        resp.names <- "tau"
        resp.start <- resp.pars
        resp.link.ids <- 3
        resp.map <- factor(seq_along(resp.start))
    }
    resp.map[resp.names %in% fix.names] <- NA
    par.names <- c(par.names, resp.names)
    pars.start <- c(pars.start, resp.start)
    link.ids <- c(link.ids, resp.link.ids)
    map[["link_resp_pars"]] <- resp.map
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
    ## Getting par.link and par.unlink.
    fixed <- par.names %in% fix.names
    par.link <- link.closure(link.ids, fixed)
    par.unlink <- unlink.closure(link.ids, fixed)
    par.dlink <- dlink.closure(link.ids, fixed)
    ## Converting parameters to link scale.
    link.pars.start <- par.link(pars.start)
    ## Total number of parameters.
    n.pars <- length(pars.start)
    ## Making vectors of fixed and unfixed parameters.
    pars.start.fixed <- pars.start[fixed]
    pars.start.unfixed <- pars.start[!(fixed)]
    link.pars.start.fixed <- link.pars.start[fixed]
    link.pars.start.unfixed <- link.pars.start[!(fixed)]
    detprob.objs <- list()
    ## Keep stuff that is important for some reason.
    keep <- rep(TRUE, n.pars)
    if (any(toa.indices != -1) & (length(levels(map$link_sigma_toa)) == 0)){
        keep[toa.indices] <- FALSE
    }
    if (any(det.indices != -1) & (length(levels(map$link_det_pars)) == 0)){
        keep[det.indices] <- FALSE
    }
    if (any(cov.indices != -1) & (length(levels(map$link_cov_pars)) == 0)){
        keep[cov.indices] <- FALSE
    }
    if (any(resp.indices != -1) & (length(levels(map$link_resp_pars)) == 0)){
        keep[resp.indices] <- FALSE
    }
    model.opts$par.names <- par.names[keep]
    fixed.keep <- fixed[keep]
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
                                                   resp_id = resp.id,
                                                   cov_id = cov.id,
                                                   mult_id = mult.id,
                                                   link_det_ids = link.ids[det.indices],
                                                   link_cov_ids = if (cov.id == 6) 0 else link.ids[cov.indices],
                                                   link_resp_ids = link.ids[resp.indices],
                                                   log_offset = 1e-10),
                                       parameters = list(link_det_pars = link.pars.start[det.indices],
                                                         link_cov_pars = if (cov.id == 6) 1 else link.pars.start[cov.indices],
                                                         link_resp_pars = link.pars.start[resp.indices],
                                                         link_sigma_toa = ifelse(toa.id, link.pars.start[toa.indices], 1),
                                                         link_D = ifelse(conditional.n | Rhess, 1, link.pars.start[D.indices]),
                                                         u = u.detprob),
                                       map = map, random = random.comp, DLL = "cov_detprob", silent = TRUE)
    }
    get.fn.gr <- function(fun = "nll"){
        function(link.pars.unfixed){
            link.pars <- numeric(n.pars)
            link.pars[fixed] <- link.pars.start.fixed
            link.pars[!fixed] <- link.pars.unfixed
            keep.tmb <- keep
            if (Rhess){
                keep.tmb[D.indices] <- FALSE
            }
            link.pars.tmb <- link.pars[keep.tmb]
            link.pars <- link.pars[keep]
            det.probs <- numeric(n.mask)
            neglog.det.probs <- numeric(n.mask)
            neglog.det.probs.grads <- matrix(0, nrow = length(link.pars.tmb), ncol = n.mask)
            for (i in 1:n.mask){
                neglog.det.probs[i] <- detprob.objs[[i]]$fn(link.pars.tmb)
                det.probs[i] <- exp(-neglog.det.probs[i])
                neglog.det.probs.grads[, i] <- detprob.objs[[i]]$gr(link.pars.tmb)
            }
            if (fun != "det.probs"){
                if (manual.sep){
                    ind.objs <- list()
                    for (i in 1:n){
                        if (cov.id == 6){
                            u.nll.ind <- 0
                        } else {
                            u.nll.ind <- u.nll[i, ]
                        }
                        if (is.null(toa)){
                            toa.ssq.ind <- 0
                        } else {
                            toa.ssq.ind <- toa.ssq[i, ]
                        }
                        ind.objs[[i]] <- MakeADFun(data = list(capt = capt[i, ],
                                                               n_dets = n.dets[i],
                                                               mask_dists = mask.dists,
                                                               trap_dists = trap.dists,
                                                               n_traps = n.traps,
                                                               n_mask = n.mask,
                                                               mask_area = mask.area,
                                                               resp_id = resp.id,
                                                               detfn_id = detfn.id,
                                                               cov_id = cov.id,
                                                               mult_id = mult.id,
                                                               det_probs = det.probs,
                                                               toa_id = toa.id,
                                                               toa_ssq = toa.ssq.ind,
                                                               link_det_ids = link.ids[det.indices],
                                                               link_cov_ids = if (cov.id == 6) 0 else link.ids[cov.indices],
                                                               link_resp_ids = link.ids[resp.indices]),
                                                   parameters = list(link_det_pars = link.pars[det.indices],
                                                                     link_cov_pars = if (cov.id == 6) 1 else link.pars[cov.indices],
                                                                     link_resp_pars = link.pars.start[resp.indices],
                                                                     link_sigma_toa = ifelse(toa.id, link.pars[toa.indices], 1),
                                                                     link_D = ifelse(conditional.n | Rhess, 1, link.pars[D.indices]),
                                                                     u = u.nll.ind),
                                                   map = map, random = random.comp, DLL = "cov_nll_sep", silent = TRUE)
                    }
                    nll.obj <- list()
                    nll.obj$fn <- function(pars){
                        nll.contribs <- numeric(n)
                        for (i in 1:n){
                            nll.contribs[i] <- ind.objs[[i]]$fn(pars)
                        }
                        if (trace){
                            cat("Calculating NLL...\n")
                        }
                        sum(nll.contribs)
                    }
                    nll.obj$gr <- function(pars){
                        gr.contribs <- matrix(0, nrow = n, ncol = length(pars))
                        for (i in 1:n){
                            gr.contribs[i, ] <- ind.objs[[i]]$gr(pars)
                        }
                        if (trace){
                            cat("Calculating gradients...\n")
                        }
                        apply(gr.contribs, 2, sum)
                    }
                } else {
                    stop("Non-manual separability no longer supported.")
                }
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
                    if (resp.id == 2){
                        cat("; nu: ", format(round(par.unlink(link.pars, resp.indices), 2), nsmall = 2), sep = "")
                    } else if (resp.id == 3){
                        cat("; size: ", format(round(par.unlink(link.pars, resp.indices), 2), nsmall = 2), sep = "")
                    } else if (resp.id == 4){
                        cat("; tau: ", format(round(par.unlink(link.pars, resp.indices), 2), nsmall = 2), sep = "")
                    }
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
                    D.indices.updated <- sum(keep[1:D.indices])
                    out[-D.indices.updated] <- nll.obj$gr(link.pars.tmb) + exp(link.pars[D.indices.updated])*mask.area*apply(neglog.det.probs.grads, 1, function(x) sum(-exp(-neglog.det.probs)*x))
                    out[D.indices.updated] <- mask.area*sum(det.probs)*exp(link.pars[D.indices.updated]) - n
                } else if (!conditional.n){
                    out <- nll.obj$gr(link.pars.tmb) + exp(link.pars[D.indices.updated])*mask.area*apply(neglog.det.probs.grads, 1, function(x) sum(-exp(-neglog.det.probs)*x))
                } else {
                    out <- nll.obj$gr(link.pars.tmb) + n*apply(neglog.det.probs.grads, 1, function(x) sum(-exp(-neglog.det.probs)*x))/sum(det.probs)    
                }
                out <- out[!fixed.keep]
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
    ## Closure to provide organisation function stuff after model fitting.
    organise.closure <- function(survey.data, model.opts, organise.fun, det.probs.fun = NULL){
        function(pars, objective){
            all.pars <- numeric(n.pars)
            all.pars[fixed] <- link.pars.start.fixed
            all.pars[!fixed] <- pars
            all.pars <- all.pars[keep]
            organise.fun(pars, all.pars, fixed, objective, survey.data, model.opts, det.probs.fun)
        }
    }
    ## Closure to provide variance-covariance calculation, given parameter estimates.
    vcov.closure <- function(survey.data, model.opts, nll, dlink.fun, gr = NULL){
        function(pars){
            hess.link <- optimHess(pars, nll, gr = gr)
            vcov.link <- solve(hess.link)
            n.pars <- length(pars)
            jacobian <- diag(n.pars)
            diag(jacobian) <- dlink.fun(pars, unfixed = TRUE)
            vcov <- jacobian %*% vcov.link %*% t(jacobian)
            dimnames(vcov) <- list(model.opts$par.names[!fixed.keep], model.opts$par.names[!fixed.keep])
            dimnames(vcov.link) <- list(paste(model.opts$par.names[!fixed.keep], "link", sep = "."),
                                        paste(model.opts$par.names[!fixed.keep], "link", sep = "."))
            out <- list(vcov = vcov, vcov.link = vcov.link)
            out
        }
    }
    obj.det.probs <- 
        obj.vcov <- vcov.closure(survey.data, model.opts, obj.fn, par.dlink, gr = obj.gr)
    obj.organise <- organise.closure(survey.data, model.opts, cov.organise, get.fn.gr(fun = "det.probs"))
    obj <- list(par = link.pars.start.unfixed, fn = obj.fn, gr = obj.gr,
                vcov = obj.vcov, organise = obj.organise, det.names = det.names,
                cov.names = cov.names)
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

## Function to organise stuff.
cov.organise <- function(pars, all.pars, fixed, objective, survey.data, model.opts, det.probs.fun){
    det.probs <- det.probs.fun(pars)
    mask.area <- survey.data$mask.area
    capt <- survey.data$capt
    esa <- sum(det.probs)*mask.area
    D <- nrow(capt)/esa
    link.ids <- model.opts$link.ids
    par.unlink <- unlink.closure(link.ids, fixed)
    ests.link <- c(all.pars, log(D))
    names(ests.link) <- paste(c(model.opts$par.names, "D"), "link", sep = ".")
    est <- c(par.unlink(all.pars), D, esa)
    names(est) <- c(model.opts$par.names, "D", "esa")
    list(ests = est, ests.link = ests.link)    
}


## Closure to provide linking function without passing link IDs.
link.closure <- function(link.ids, fixed){
    function(pars, which = NULL, unfixed = FALSE){
        if (unfixed){
            link.ids <- link.ids[!fixed]
        }
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

## Closure to provide unlinking function without passing link IDs.
unlink.closure <- function(link.ids, fixed){
    function(link.pars, which = NULL, unfixed = FALSE){
        if (unfixed){
            link.ids <- link.ids[!fixed]
        }
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

## Closure to provide the first derivative of the linking function without passing link IDs.
dlink.closure <- function(link.ids, fixed){
    function(link.pars, which = NULL, unfixed = FALSE){
        if (unfixed){
            link.ids <- link.ids[!fixed]
        }
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

links <- list(log, qlogis, identity, function(x) log(x - 1))
unlinks <- list(exp, plogis, identity, function(x) exp(x) + 1)
dlinks <- list(exp, dlogis, function(x) rep(1, length(x)), exp)
