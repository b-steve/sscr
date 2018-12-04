context("Testing model fits")

test_that(
    "Simulation,",
    {
        compile.sscr()
        ## Poisson response.
        set.seed(4321)
        sim.data <- sim.sscr(traps = test.data$traps, mask = test.data$mask, D = 1.5,
                             resp = "pois", detfn = "hhn",
                             cov.structure = "none", det.pars = list(lambda0 = 3, sigma = 50))
        expect_equivalent(sim.data[21, ], c(0, 0, 0, 0, 1, 4, 0, 4, 1))
        ## ... With detection function on the probability scale.
        set.seed(2468)
        sim.data <- sim.sscr(traps = test.data$traps, mask = test.data$mask, D = 1.5,
                             resp = "pois", detfn = "hn",
                             cov.structure = "none", det.pars = list(g0 = 0.9, sigma = 50))
        expect_equivalent(sim.data[11, ], c(0, 0, 0, 0, 2, 1, 0, 0, 0))
        ## Binomial response.
        set.seed(1234)
        sim.data <- sim.sscr(traps = test.data$traps, mask = test.data$mask, D = 1.5,
                             resp = "binom", resp.pars = 5, detfn = "hhn",
                              cov.structure = "none", det.pars = list(lambda0 = 3, sigma = 50))
        expect_equivalent(sim.data[8, ], c(1, 0, 0, 5, 1, 0, 1, 0, 0))
        ## With squared-exponential covariance function.
        set.seed(8642)
        sim.data <- sim.sscr(traps = test.data$traps, mask = test.data$mask, D = 1.5,
                             resp = "pois", detfn = "hhn",
                             cov.structure = "sq_exponential",, re.scale = "er",
                             det.pars = list(lambda0 = 3, sigma = 50),
                             cov.pars = list(sigma.u = 1.5, rho = 100))
        expect_equivalent(sim.data[13, ], c(0, 0, 0, 2, 5, 0, 1, 16, 0))
        ## With exponential covariance function on the probability scale.
        set.seed(8642)
        sim.data <- sim.sscr(traps = test.data$traps, mask = test.data$mask, D = 1.5,
                             resp = "pois", detfn = "hhn",
                             cov.structure = "exponential", re.scale = "prob",
                             det.pars = list(lambda0 = 3, sigma = 50),
                             cov.pars = list(sigma.u = 1.5, rho = 100))
        expect_equivalent(sim.data[11, ], c(0, 0, 0, 2, 4, 1, 2, 5, 1))
        ## With counts from an OU process.
        set.seed(7531)
        sim.data <- sim.sscr(traps = test.data$traps, mask = test.data$mask, D = 0.5,
                             cov.structure = "OU",
                             cov.pars = list(tau = 30, sigma = 50, n.steps = 100, epsilon = 10))
        expect_equivalent(sim.data[1, ], c(0, 1, 0, 0, 4, 0, 0, 0, 0))
        ## With time-of-arrival information.
        set.seed(3579)
        sim.data <- sim.sscr(traps = test.data$traps, mask = test.data$mask, D = 1.5,
                             resp = "binom", detfn = "hn",
                             cov.structure = "sq_exponential", re.scale = "er",
                             det.pars = list(g0 = 0.9, sigma = 50),
                             cov.pars = list(sigma.u = 1.5, rho = 100),
                             toa.pars = list(sigma.toa = 4, sound.speed = 330))
        expect_equivalent(sim.data$toa[6, ], c(0, 0, 0, 0, 0.201980061428043,
                                               0.20259893914751, 0, 0.2295964076072,
                                               0.223745109102043))
    })
    

test_that(
    "Fits with no random effects.",
    {
        ## Poisson.
        fit.pois <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, resp = "pois",
                             detfn = "hhn", hess = TRUE)
        expect_equivalent(fit.pois$ests, c(2.87645050233774, 93.388169123886, 0.454913613852446, 
                                           30.7750737144151, -128.981698520115), tol = 1e-4)
        expect_equivalent(fit.pois$se, c(0.5405034, 20.6192798, 0.1894562), tol = 1e-4)
        ## Poisson with detection function on probability scale.
        fit.pois.pr <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, resp = "pois", detfn = "hn")
        expect_equivalent(fit.pois.pr$est, c(0.999999999795205, 117.513845050955, 0.431769352119301, 
                                             32.424719196215, -124.387177667554), tol = 1e-4)
        ## Bernoulli.
        fit.bern <- fit.sscr(test.data$bin.capt, test.data$traps, test.data$mask, detfn = "hhn")
        expect_equivalent(fit.bern$ests, c(1.10544298833934, 144.339009111062, 0.312312925517947, 
                                           44.8268350622443, -73.7787592050756), tol = 1e-4)
        ## Bernoulli with detection function on probability scale.
        fit.bern.pr <- fit.sscr(test.data$bin.capt, test.data$traps, test.data$mask, detfn = "hn")
        expect_equivalent(fit.bern.pr$ests, c(0.709787201866111, 157.488436993617, 0.308950967387457, 
                                              45.3146339640442, -74.0448322726064), tol = 1e-4)
        ## Binomial.
        fit.binom <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, resp = "binom",
                              resp.pars = 10, detfn = "hhn")
        expect_equivalent(fit.binom$ests, c(0.519014374421606, 125.87869591362, 0.247980358116291, 
                                            56.4560842896866, -136.734734751027), tol = 1e-4)
        ## Binomial with detection function on probability scale.
        fit.binom.pr <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, resp = "binom",
                                 resp.pars = 10, detfn = "hn")
        expect_equivalent(fit.binom.pr$ests, c(0.604976358068237, 137.677620551755, 0.20620406615787, 
                                               67.8939084997556, -135.548541895984), tol = 1e-4)
    })

test_that(
    "Tests with random effects.",
    {   
        ## Testing exact gradients.
        fit.grad <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                             mask = test.data$mask, resp = "pois", detfn = "hhn",
                             test.conditional.n = FALSE, cov.structure = "exponential",
                             test = "gr")
        expect_equivalent(fit.grad$gr, c(53.936454081332, 424.759615492023, 41.4534610100935, -9.18292928534691, 
                                         260.095722145605),
                          tolerance = 1e-4, scale = 1)
        ## ... without manual separability.
        fit.grad.man <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                 mask = test.data$mask, resp = "pois", detfn = "hhn",
                                 test.conditional.n = FALSE, cov.structure = "exponential",
                                 test = "gr", manual.sep = FALSE)
        expect_equivalent(fit.grad.man$gr, c(53.936454081332, 424.759615492023, 41.4534610100935,
                                             -9.18292928534691, 260.095722145605),
                          tolerance = 1e-4, scale = 1)
        ## Independent random effects.
        test.ind.hn <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois", detfn = "hhn",
                                cov.structure = "independent", test = "nll")
        expect_equivalent(test.ind.hn$nll, 127.0738, tolerance = 1e-4, scale = 1)
        ## With detection function on probability scale.
        test.ind.hn.detpr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                   mask = test.data$mask, resp = "pois", detfn = "hn",
                                   cov.structure = "independent", test = "nll")
        expect_equivalent(test.ind.hn.detpr$nll, 123.8944, tolerance = 1e-4, scale = 1)
        ## With random effects on the probability scale.
        test.ind.hn.repr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                     mask = test.data$mask, resp = "pois", detfn = "hhn",
                                     cov.structure = "independent", re.scale = "prob",
                                     test = "nll")
        expect_equivalent(test.ind.hn.repr, 124.5118, tolerance = 1e-4, scale = 1)
        ## With detection function and random effects also on probability scale.
        test.ind.hn.detpr.repr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                           mask = test.data$mask, resp = "pois", detfn = "hn",
                                           cov.structure = "independent",
                                           re.scale = "prob", test = "nll")
        expect_equivalent(test.ind.hn.detpr.repr$nll, 129.4880, tolerance = 1e-4, scale = 1)
        ## Testing manual start values.
        test.sv.hn <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "pois", detfn = "hhn",
                               cov.structure = "independent",
                               start = c(lambda0 = 2, sigma.u = 4), test = "nll")
        expect_equivalent(test.sv.hn, 149.6288, tolerance = 1e-4, scale = 1)
        ## ... with hazard rate detection function.
        test.ind.hr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois", detfn = "hhr",
                                cov.structure = "independent", test = "nll")
        expect_equivalent(test.ind.hr, 152.7691, tolerance = 1e-4, scale = 1)
        ## Exponential covariance function.
        test.exp <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                             mask = test.data$mask, resp = "pois", detfn = "hhn",
                             cov.structure = "exponential", test = "nll")
        expect_equivalent(test.exp$nll, 125.2871, tolerance = 1e-4, scale = 1)
        ## ... without manual separability.
        test.exp.man <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                 mask = test.data$mask, resp = "pois", detfn = "hhn",
                                 cov.structure = "exponential", test = "nll",
                                 manual.sep = FALSE)
        expect_equivalent(test.exp.man$nll, 125.2871, tolerance = 1e-4, scale = 1)
        ## With detection function on probability scale.
        test.exp.detpr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois", detfn = "hn",
                                cov.structure = "exponential",
                                test = "nll")
        expect_equivalent(test.exp.detpr$nll, 122.8905, tolerance = 1e-4, scale = 1)
        ## With random effects on probability scale.
        test.exp.repr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                  mask = test.data$mask, resp = "pois", detfn = "hhn",
                                  cov.structure = "exponential", re.scale = "prob",
                                  test = "nll")
        expect_equivalent(test.exp.repr$nll, 125.0680, tolerance = 1e-4, scale = 1)
        ## Squared exponential covariance function (not sure we can trust this one).
        test.sqexp <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "pois", detfn = "hhn",
                               cov.structure = "sq_exponential", test = "nll")
        expect_equivalent(test.sqexp$nll, 124.6839, tolerance = 1e-4, scale = 1)
        ## Individual random effect.
        test.indiv <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "pois", detfn = "hhn",
                               cov.structure = "individual", test = "nll")
        expect_equivalent(test.indiv$nll, 125.1148, tolerance = 1e-4, scale = 1)
        ## Bernoulli response.
        test.bern <- fit.sscr(capt = test.data$bin.capt, traps = test.data$traps,
                              mask = test.data$mask, detfn = "hhn",
                              cov.structure = "exponential", test = "nll")
        expect_equivalent(test.bern$nll, 77.4819, tolerance = 1e-4, scale = 1)
        ## With detection function on probability scale.
        test.bern.detpr <- fit.sscr(capt = test.data$bin.capt, traps = test.data$traps,
                                    mask = test.data$mask, detfn = "hn",
                                    cov.structure = "exponential", test = "nll")
        expect_equivalent(test.bern.detpr$nll, 76.1733, tolerance = 1e-4, scale = 1)
        ## With random effects on probability scale.
        test.bern.repr <- fit.sscr(capt = test.data$bin.capt, traps = test.data$traps,
                                   mask = test.data$mask, detfn = "hhn",
                                   cov.structure = "exponential", re.scale = "prob",
                                   test = "nll")
        expect_equivalent(test.bern.repr$nll, 77.7858, tolerance = 1e-4, scale = 1)
        ## Binomial response.
        test.binom <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "binom",
                               resp.pars = 10, detfn = "hhn",
                               cov.structure = "exponential", test = "nll")
        expect_equivalent(test.binom$nll, 123.5817, tolerance = 1e-4, scale = 1)
    })
