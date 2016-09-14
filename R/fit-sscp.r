#' Fitting an SCR model with second-order spatial dependence
#'
#' Fits an SSCR model. Estimation is by maximum likelihood. The
#' second-order spatial dependence is modelled via trap-level random
#' effects for each detected individual. The likelihood function is
#' calculated by integrating over these random effects using the
#' Laplace approximation.
#'
#' @param capt Capture history object.
#' @param traps Traps object.
#' @param mask Mask object.
#' @param resp Response distribution for capture history elements.
#' @param re.structure Structure of 
#' 
#' @export
fit.sscr <- function(capt, traps, mask, resp = "binom", re.structure = "none"){
    ## Loading DLLs.
    dll.dir <- paste(system.file(package = "sscr"), "/tmb/bin/", sep = "")
    for (i in paste(dll.dir, list.files(dll.dir), sep = "")){
        dyn.load(i)
    }
    ## Calculating mask distances.
    mask.dists <- crossdist(mask[, 1], mask[, 2],
                            traps[, 1], traps[, 2])
    ## Extracting mask area.
    mask.area <- attr(mask, "area")
    ## Number of mask points and traps.
    n.mask <- nrow(mask)
    n.traps <- nrow(traps)
    ## Packaging the data up into a list.
    survey.data <- list(capt = capt,
                        mask.dists = mask.dists,
                        mask.area = mask.area,
                        n.mask = n.mask,
                        n.traps = n.traps)
    ## Optimisation object constructor function.
    make.obj <- switch(re.structure,
                       none = make.obj.none,
                       independent = make.obj.independent,
                       make.obj.error)
    ## Making optimisation object.
    opt.obj <- make.obj(survey.data, resp)
    ## Fitting model.
    fit <- nlminb(opt.obj$par, opt.obj$fn, opt.obj$gr)
    fit
}

