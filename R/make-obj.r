## Functions to make optimisation objects.
make.obj.none <- function(survey.data, model.opts){
    ## Extracting data.
    capt <- survey.data$capt
    mask.dists <- survey.data$mask.dists
    mask.area <- survey.data$mask.area
    n.mask <- survey.data$n.mask
    n.traps <- survey.data$n.traps
    n <- nrow(capt)
    resp <- model.opts$resp
    detfn <- model.opts$detfn
    detfn.id <- switch(detfn, hn = 0, hr = 1)
    ## Packaging data for TMB template.
    data <- list(capt = capt,
                 mask_dists = mask.dists,
                 n = n,
                 n_traps = n.traps,
                 n_mask = n.mask,
                 mask_area = mask.area,
                 resp_id = switch(resp, binom = 0, pois = 1, count = 1),
                 detfn_id = detfn.id)
    ## Start values for optimisation.
        ## Indices and start values for detection function parameters.
    if (detfn.id == 0){
        det.start <- c(max(capt)/2, max(apply(mask.dists, 1, min))/5)
    } else if (resp.id == 1){
        det.start <- c(max(capt)/2, max(apply(mask.dists, 1, min))/5, 1)
    }
    log.det.pars <- log(det.start)
    ## Making optimisation object with TMB.
    obj <- MakeADFun(data = data, parameters = list(log_det_pars = log.det.pars),
                     DLL = "simple_nll", silent = TRUE)
    obj
}

make.obj.cov <- function(survey.data, model.opts){
    ## Extracting data.
    capt <- survey.data$capt
    mask.dists <- survey.data$mask.dists
    mask.area <- survey.data$mask.area
    n.mask <- survey.data$n.mask
    n.traps <- survey.data$n.traps
    trap.dists <- survey.data$trap.dists
    n <- nrow(capt)
    ## Indicator for response type.
    resp.id <- switch(model.opts$resp, binom = 0, pois = 1)
    model.opts$resp.id <- resp.id
    ## Indicator for detection functoin.
    detfn.id <- switch(model.opts$detfn, hn = 0, hr = 1)
    model.opts$detfn.id <- detfn.id
    ## Indicator for covariance structure.
    cov.id <- switch(model.opts$cov.structure,
                     independent = 0,
                     exponential = 1,
                     matern = 2,
                     full = 3)
    model.opts$cov.id <- cov.id
    ## Indices and start values for detection function parameters.
    if (resp.id == 0){
        det.indices <- 1:2
        det.start <- c(max(capt)/2, max(apply(mask.dists, 1, min))/5)
    } else if (resp.id == 1){
        det.indices <- 1:3
        det.start <- c(max(capt)/2, max(apply(mask.dists, 1, min))/5, 1)
    }
    ## Indices and start values for covariance parameters.
    cov.index.start <- max(det.indices) + 1
    if (cov.id == 0){
        cov.indices <- cov.index.start
        cov.start <- sd(capt)
    } else if (cov.id == 1){
        cov.indices <- cov.index.start:(cov.index.start + 1)
        cov.start <- c(sd(capt), mean(trap.dists))
    } else if (cov.id == 2){
        cov.indices <- cov.index.start:(cov.index.start + 2)
        cov.start <- c(1, 1, 1)
    } else if (cov.id == 3){
        cov.indices <- cov.index.start
        cov.start <- sd(capt)
    }
    model.opts$det.indices <- det.indices
    model.opts$cov.indices <- cov.indices
    pars <- log(c(det.start, cov.start))
    ## Optimisation object.
    list(par = pars, fn = nll.closure(pars, survey.data, model.opts, cov.nll),
         organise = organise.closure(pars, survey.data, model.opts, cov.organise))
}

make.obj.error <- function(survey.data, resp){
    stop("Random effects structure not recognised. Use \"none\" or \"independent.\"")
}

## Closure to provide negative log-likelihood function without passing data.
nll.closure <- function(pars, survey.data, model.opts, nll.fun){
    function(pars){
        nll.fun(pars, survey.data, model.opts)
    }
}

organise.closure <- function(fit, survey.data, model.opts, organise.fun){
    function(fit){
        organise.fun(fit, survey.data, model.opts)
    }
}
