## Functions to make optimisation objects.
make.obj.none <- function(survey.data, model.opts){
    ## Extracting data.
    capt <- survey.data$capt
    mask.dists <- survey.data$mask.dists
    mask.area <- survey.data$mask.area
    n.mask <- survey.data$n.mask
    n.traps <- survey.data$n.traps
    trap.dists <- survey.data$trap.dists
    trace <- survey.data$trace
    n <- nrow(capt)
    resp <- model.opts$resp
    resp.pars <- model.opts$resp.pars
    if (is.null(resp.pars)){
        resp.pars <- 1
    }
    detfn <- model.opts$detfn
    detfn.id <- switch(detfn, hn = 0, hr = 1)
    detfn.scale <- model.opts$detfn.scale
    detfn.scale.id <- switch(detfn.scale, er = 0, prob = 1)
    ## Start values for optimisation.
    start <- model.opts$start
    start.names <- names(start)
    ## Indices and start values for detection function parameters.
    ## Link 0 is log, 1 is qlogis.
    ## For detection functions on the hazard scale.
    if (detfn.scale.id == 0){
        if (detfn.id == 0){
            ## .. With a halfnormal detection function.
            det.start <- numeric(2)
            det.start[1] <- ifelse(any(start.names == "lambda0") , start["lambda0"],
                                   max(capt)/(2*resp.pars))
            det.start[2] <- ifelse(any(start.names == "sigma") , start["sigma"],
                                   min(trap.dists[trap.dists > 0]))
            link.ids <- c(0, 0)
        } else if (detfn.id == 1){
            ## .. With a hazard rate detection function.
            det.start <- numeric(3)
            det.start[1] <- ifelse(any(start.names == "lambda0") , start["lambda0"],
                                   max(capt)/(2*resp.pars))
            det.start[2] <- ifelse(any(start.names == "sigma") , start["sigma"],
                                   min(trap.dists[trap.dists > 0]))
            det.start[3] <- ifelse(any(start.names == "z") , start["z"], 1)                   
            link.ids <- c(0, 0, 0)
        }
        ## For detection functions on the probability scale.
    } else if (detfn.scale.id == 1){
        if (detfn.id == 0){
            ## .. With a halfnormal detection function.
            det.start <- numeric(2)
            det.start[1] <- ifelse(any(start.names == "g0") , start["g0"], 0.5)
            det.start[2] <- ifelse(any(start.names == "sigma") , start["sigma"],
                                   min(trap.dists[trap.dists > 0]))
            link.ids <- c(1, 0)
        } else if (detfn.id == 1){
            ## .. With a hazard rate detection function.
            det.start <- numeric(3)
            det.start[1] <- ifelse(any(start.names == "g0") , start["g0"], 0.5)
            det.start[2] <- ifelse(any(start.names == "sigma") , start["sigma"],
                                   min(trap.dists[trap.dists > 0]))
            det.start[3] <- ifelse(any(start.names == "z") , start["z"], 1)           
            link.ids <- c(1, 0, 0)
        }
    }
    ## Getting par.link and par.unlink.
    par.link <- link.closure(link.ids)
    par.unlink <- unlink.closure(link.ids)
    link.det.pars <- par.link(det.start)
    ## Packaging data for TMB template.
    data <- list(capt = capt,
                 mask_dists = mask.dists,
                 n = n,
                 n_traps = n.traps,
                 n_mask = n.mask,
                 mask_area = mask.area,
                 resp_id = switch(resp, binom = 0, pois = 1, count = 1),
                 resp_pars = resp.pars,
                 detfn_id = detfn.id,
                 detfn_scale_id = detfn.scale.id,
                 link_ids = link.ids)
    ## Making optimisation object with TMB.
    obj <- MakeADFun(data = data, parameters = list(link_det_pars = link.det.pars),
                     DLL = "simple_nll", silent = TRUE)
    if (trace){
        obj$fn.notrace <- obj$fn
        obj$fn <- function(x, ...){
            out <- obj$fn.notrace(x, ...)
            cat("Detection parameters: ", paste(format(round(par.unlink(x), 2), nsmall = 2), collapse = ", "),
                "; nll: ", format(round(as.numeric(out), 2), nsmall = 2), "\n", sep = "")
            out
        }
    }
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
    resp.pars <- model.opts$resp.pars
    if (is.null(resp.pars)){
        resp.pars <- 1
    }
    model.opts$resp.pars <- resp.pars
    ## Indicator for detection functoin.
    detfn.id <- switch(model.opts$detfn, hn = 0, hr = 1)
    model.opts$detfn.id <- detfn.id
    ## Indicator for covariance structure.
    cov.id <- switch(model.opts$cov.structure,
                     independent = 0,
                     exponential = 1,
                     matern = 2,
                     individual = 3,
                     lc_exponential = 4,
                     sq_exponential = 5)
    model.opts$cov.id <- cov.id
    start <- model.opts$start
    start.names <- names(start)
    ## Indices and start values for detection function parameters.
    if (detfn.id == 0){
        det.indices <- 1:2
        det.start <- numeric(2)
        det.start[1] <- ifelse(any(start.names == "lambda0"),
                               start["lambda0"],
                               max(capt)/(2*resp.pars))
        det.start[2] <- ifelse(any(start.names == "sigma"),
                               start["sigma"],
                               max(apply(mask.dists, 1, min))/5)
    } else if (detfn.id == 1){
        det.indices <- 1:3
        det.start <- numeric(3)
        det.start[1] <- ifelse(any(start.names == "lambda0"),
                               start["lambda0"],
                               max(capt)/(2*resp.pars))
        det.start[2] <- ifelse(any(start.names == "sigma"),
                               start["sigma"],
                               max(apply(mask.dists, 1, min))/5)
        det.start[3] <- ifelse(any(start.names == "z"),
                               start["z"], 1)
    }
    ## Indices and start values for covariance parameters.
    cov.index.start <- max(det.indices) + 1
    if (cov.id == 0){
        ## Independent.
        cov.indices <- cov.index.start
        cov.start <- numeric(1)
        cov.start[1] <- ifelse(any(start.names == "sigma.u"),
                               start["sigma.u"], sd(capt))
    } else if (cov.id == 1){
        ## Exponential.
        cov.indices <- cov.index.start:(cov.index.start + 1)
        cov.start <- numeric(2)
        cov.start[1] <- ifelse(any(start.names == "sigma.u"),
                               start["sigma.u"], sd(capt))
        cov.start[2] <- ifelse(any(start.names == "rho"),
                               start["rho"], mean(trap.dists))
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

## Closure to provide linking function without passing link ids.
link.closure <- function(link.ids){
    function(pars){
        n.pars <- length(pars)
        out <- numeric(n.pars)
        for (i in 1:n.pars){
            out[i] <- links[[link.ids[i] + 1]](pars[i])
        }
        out
    }
}

unlink.closure <- function(link.ids){
    function(link.pars){
        n.pars <- length(link.pars)
        out <- numeric(n.pars)
        for (i in 1:n.pars){
            out[i] <- unlinks[[link.ids[i] + 1]](link.pars[i])
        }
        out
    }
}

links <- list(log, qlogis)
unlinks <- list(exp, plogis)

## par.link <- function(pars, link.ids){
##     n.pars <- length(pars)
##     out <- numeric(n.pars)
##     for (i in 1:n.pars){
##         out[i] <- links[[link.ids[i] + 1]](pars[i])
##     }
##     out
## }

## par.unlink <- function(link.pars, link.ids){
##     n.pars <- length(link.pars)
##     out <- numeric(n.pars)
##     for (i in 1:n.pars){
##         out[i] <- unlinks[[link.ids[i] + 1]](link.pars[i])
##     }
##     out
## }
