## Functions to make optimisation objects.
make.obj.none <- function(survey.data, resp){
    ## Extracting data.
    capt <- survey.data$capt
    mask.dists <- survey.data$mask.dists
    mask.area <- survey.data$mask.area
    n.mask <- survey.data$n.mask
    n.traps <- survey.data$n.traps
    n <- nrow(capt)
    ## Packaging data for TMB template.
    data <- list(capt = capt,
                 mask_dists = mask.dists,
                 n = n,
                 n_traps = n.traps,
                 n_mask = n.mask,
                 mask_area = mask.area,
                 resp_id = switch(resp, binom = 0, pois = 1))
    ## Start values for optimisation.
    pars <- list(log_lambda0 = log(n/mask.area),
                 log_sigma = log(max(apply(mask.dists, 1, min))/5))
    ## Making optimisation object with TMB.
    MakeADFun(data = data, parameters = pars, DLL = "simple_nll", silent = TRUE)
}

make.obj.independent <- function(survey.data, resp){
    ## Extracting data.
    ## Extracting data.
    capt <- survey.data$capt
    mask.dists <- survey.data$mask.dists
    mask.area <- survey.data$mask.area
    n.mask <- survey.data$n.mask
    n.traps <- survey.data$n.traps
    n <- nrow(capt)
    ## Start values for optimisation.
    #pars <- c(log.lambda0 = log(n/mask.area),
    #          log.sigma = log(max(apply(mask.dists, 1, min))/5),
    #          log.sigma.u = log(1))
    pars <- c(log.lambda0 = log(2),
              log.sigma = log(100),
              log.sigma.u = log(1))
    ## Optimisation object.
    list(par = pars, fn = nll.closure(pars, survey.data, ind.nll))
}

make.obj.error <- function(survey.data, resp){
    stop("Random effects structure not recognised. Use \"none\" or \"independent.\"")
}

## Closure to provide negative log-likelihood function without passing data.
nll.closure <- function(pars, survey.data, nll.fun){
    function(pars){
        nll.fun(pars, survey.data)
    }
}
