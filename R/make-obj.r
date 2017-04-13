## Functions to make optimisation objects.
make.obj <- function(survey.data, model.opts, any.cov){
    ## Extracting data.
    capt <- survey.data$capt
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
    ## Extracting trace.
    trace <- survey.data$trace
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
        } else if (detfn.id == 1){
            ## .. With a hazard rate detection function.
            det.indices <- 1:3
            det.start <- numeric(3)
            det.start[1] <- ifelse(any(start.names == "g0") , start["g0"], 0.5)
            det.start[2] <- ifelse(any(start.names == "sigma"), start["sigma"],
                                   max(apply(mask.dists, 1, min))/5)
            det.start[3] <- ifelse(any(start.names == "z") , start["z"], 1)           
            det.link.ids <- c(1, 0, 0)
        }
    }
    pars.start <- det.start
    link.ids <- det.link.ids
    if (any.cov){
        ## Stuff only for covariance structures.
        cov.structure <- model.opts$cov.structure
        cov.id <- switch(model.opts$cov.structure,
                         independent = 0,
                         exponential = 1,
                         matern = 2,
                         individual = 3,
                         lc_exponential = 4,
                         sq_exponential = 5)
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
        } else if (cov.id == 1){
            ## Exponential.
            cov.indices <- cov.index.start:(cov.index.start + 1)
            cov.start <- numeric(2)
            cov.start[1] <- ifelse(any(start.names == "sigma.u"),
                                   start["sigma.u"], sd(capt))
            cov.start[2] <- ifelse(any(start.names == "rho"),
                                   start["rho"], mean(trap.dists))
            cov.link.ids <- c(0, 0)
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
        }
        model.opts$resp.id <- resp.id
        model.opts$resp.pars <- resp.pars
        model.opts$detfn.id <- detfn.id
        model.opts$det.indices <- det.indices
        model.opts$cov.indices <- cov.indices
        pars.start <- c(pars.start, cov.start)
        link.ids <- c(link.ids, cov.link.ids)
        model.opts$link.ids <- link.ids
        model.opts$cov.id <- cov.id
        model.opts$detfn.scale.id <- detfn.scale.id
        model.opts$re.scale.id <- re.scale.id
    }
    ## Getting par.link and par.unlink.
    par.link <- link.closure(link.ids)
    par.unlink <- unlink.closure(link.ids)
    ## Converting parameters to link scale.
    link.pars.start <- par.link(pars.start)
    ## Creating required object.
    if (any.cov){
        obj <-  list(par = link.pars.start, fn = nll.closure(survey.data,
                                                             model.opts, cov.nll),
                     organise = organise.closure(survey.data, model.opts, cov.organise))
    } else {
        ## Packaging data for TMB template.
        data <- list(capt = capt,
                     mask_dists = mask.dists,
                     n = n,
                     n_traps = n.traps,
                     n_mask = n.mask,
                     mask_area = mask.area,
                     resp_id = resp.id,
                     resp_pars = resp.pars,
                     detfn_id = detfn.id,
                     detfn_scale_id = detfn.scale.id,
                     link_ids = link.ids)
        ## Making optimisation object with TMB.
        obj <- MakeADFun(data = data, parameters = list(link_det_pars = link.pars.start),
                         DLL = "simple_nll", silent = TRUE)
        ## Making function for trace.
        if (trace){
            obj$fn.notrace <- obj$fn
            obj$fn <- function(x, ...){
                out <- obj$fn.notrace(x, ...)
                cat("Detection parameters: ", paste(format(round(par.unlink(x), 2), nsmall = 2), collapse = ", "),
                    "; nll: ", format(round(as.numeric(out), 2), nsmall = 2), "\n", sep = "")
                out
            }
        }
    }
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

organise.closure <- function(survey.data, model.opts, organise.fun){
    function(fit){
        organise.fun(fit, survey.data, model.opts)
    }
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

links <- list(log, qlogis)
unlinks <- list(exp, plogis)
