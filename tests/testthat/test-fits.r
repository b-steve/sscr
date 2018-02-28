context("Testing model fits")

test_that(
    "Fits with no random effects.",
    {
        compile.sscr()
        ## Poisson.
        fit.pois <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "pois", hess = TRUE)
        expect_equivalent(fit.pois$ests, c(2.87645050233774, 93.388169123886, 0.454913613852446, 
                                           30.7750737144151, -128.981698520115), tol = 1e-4)
        expect_equivalent(fit.pois$se, c(0.5405034, 20.6192798, 0.1894562), tol = 1e-4)
        ## Poisson with detection function on probability scale.
        fit.pois.pr <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "pois", detfn.scale = "prob")
        expect_equivalent(fit.pois.pr$est, c(0.999999999795205, 117.513845050955, 0.431769352119301, 
                                             32.424719196215, -124.387177667554), tol = 1e-4)
        ## Bernoulli.
        fit.bern <- fit.sscr(test.data$bin.capt, test.data$traps, test.data$mask)
        expect_equivalent(fit.bern$ests, c(1.10544298833934, 144.339009111062, 0.312312925517947, 
                                           44.8268350622443, -73.7787592050756), tol = 1e-4)
        ## Bernoulli with detection function on probability scale.
        fit.bern.pr <- fit.sscr(test.data$bin.capt, test.data$traps, test.data$mask, detfn.scale = "prob")
        expect_equivalent(fit.bern.pr$ests, c(0.709787201866111, 157.488436993617, 0.308950967387457, 
                                              45.3146339640442, -74.0448322726064), tol = 1e-4)
        ## Binomial.
        fit.binom <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "binom", 10)
        expect_equivalent(fit.binom$ests, c(0.519014374421606, 125.87869591362, 0.247980358116291, 
                                            56.4560842896866, -136.734734751027), tol = 1e-4)
        ## Binomial with detection function on probability scale.
        fit.binom.pr <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "binom", 10, detfn.scale = "prob")
        expect_equivalent(fit.binom.pr$ests, c(0.604976358068237, 137.677620551755, 0.20620406615787, 
                                               67.8939084997556, -135.548541895984), tol = 1e-4)
    })

test_that(
    "Tests with random effects.",
    {   
        ## Full fitted model with maximisation and variance-covariance matrix.
        ## full.fit <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
        ##                      mask = test.data$mask, resp = "pois", detfn = "hn",
        ##                      cov.structure = "exponential", new = TRUE, hess = TRUE,
        ##                      Rhess = TRUE, trace = TRUE)
        ## > full.fit
        ## $pars
        ##      lambda0        sigma      sigma.u          rho            D          esa 
        ##    1.4623049  108.9016905    0.8922004   65.0390218    0.4138970   33.8248443 
        ##           LL 
        ## -118.4125173 
        
        ## $vcov
        ##              [,1]       [,2]         [,3]         [,4]         [,5]
        ## [1,]  0.339215087  -5.151117 -0.059296876   -1.2579096  0.005124181
        ## [2,] -5.151116994 789.518438  1.266556730  141.7894338 -4.398424594
        ## [3,] -0.059296876   1.266557  0.042107043    3.3862370 -0.004564574
        ## [4,] -1.257909642 141.789434  3.386237045 6321.0621113  0.249759871
        ## [5,]  0.005124181  -4.398425 -0.004564574    0.2497599  0.038847777
        
        ## Testing exact gradients.
        fit.grad <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                              mask = test.data$mask, resp = "pois", detfn = "hn",
                             cov.structure = "exponential", test = "gr")
        expect_equivalent(fit.grad$gr, c(59.6122557885348, 419.319087122365, 29.1110101667589, -13.3634887263877, 
                                         255.168266193995), tol = 1e-4)
        ## Independent random effects.
        test.ind.hn <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois", detfn = "hn",
                                cov.structure = "independent", test = "nll", trace = TRUE)
        expect_that(abs(test.ind.hn$nll - 124.8294) < 1e-4, is_true())
        ## With detection function on probability scale.
        test.ind.hn.detpr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                   mask = test.data$mask, resp = "pois", detfn = "hn",
                                   detfn.scale = "prob", cov.structure = "independent",
                                   test = "nll")
        expect_that(abs(test.ind.hn.detpr$nll - 121.65) < 1e-4, is_true())
        ## With random effects on the probability scale.
        test.ind.hn.repr <- test.ind.hn <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                                    mask = test.data$mask, resp = "pois", detfn = "hn",
                                                    cov.structure = "independent", re.scale = "prob",
                                                    test = "nll")
        expect_that(abs(test.ind.hn.repr$nll - 122.2674) < 1e-4, is_true())
        ## With detection function and random effects also on probability scale.
        test.ind.hn.detpr.repr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                           mask = test.data$mask, resp = "pois", detfn = "hn",
                                           detfn.scale = "prob", cov.structure = "independent",
                                           re.scale = "prob", test = "nll")
        expect_that(abs(test.ind.hn.detpr.repr$nll - 127.2436) < 1e-4, is_true())
        ## Testing manual start values.
        test.sv.hn <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "pois", detfn = "hn",
                               cov.structure = "independent",
                               start = c(lambda0 = 2, sigma.u = 4), test = "nll")
        expect_that(abs(test.sv.hn$nll - 147.3844) < 1e-4, is_true())
        ## ... with hazard rate detection function.
        test.ind.hr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois", detfn = "hr",
                                cov.structure = "independent", test = "nll")
        expect_that(abs(test.ind.hr$nll - 150.5247) < 1e-4, is_true())
        ## Exponential covariance function.
        test.exp <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                             mask = test.data$mask, resp = "pois",
                             cov.structure = "exponential", test = "nll")
        expect_that(abs(test.exp$nll - 122.9533) < 1e-4, is_true())
        ## With detection function on probability scale.
        test.exp.detpr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois",
                                cov.structure = "exponential", detfn.scale = "prob",
                                test = "nll")
        expect_that(abs(test.exp.detpr$nll - 120.397) < 1e-4, is_true())
        ## With random effects on probability scale.
        test.exp.repr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                  mask = test.data$mask, resp = "pois",
                                  cov.structure = "exponential", re.scale = "prob",
                                  test = "nll")
        expect_that(abs(test.exp.repr$nll - 123.3313) < 1e-4, is_true())
        ## Squared exponential covariance function (not sure we can trust this one).
        test.sqexp <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "pois",
                               cov.structure = "sq_exponential", test = "nll")
        expect_that(abs(test.sqexp$nll - 122.4532) < 1e-4, is_true())
        ## Individual random effect.
        test.indiv <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "pois",
                               cov.structure = "individual", test = "nll")
        expect_that(abs(test.indiv$nll - 122.8704) < 1e-4, is_true())
        ## Bernoulli response.
        test.bern <- fit.sscr(capt = test.data$bin.capt, traps = test.data$traps,
                              mask = test.data$mask, cov.structure = "exponential",
                              test = "nll")
        expect_that(abs(test.bern$nll - 75.09304) < 1e-4, is_true())
        ## With detection function on probability scale.
        test.bern.detpr <- fit.sscr(capt = test.data$bin.capt, traps = test.data$traps,
                                 mask = test.data$mask,  cov.structure = "exponential",
                                 detfn.scale = "prob", test = "nll")
        expect_that(abs(test.bern.detpr$nll - 73.86672) < 1e-4, is_true())
        ## With random effects on probability scale.
        test.bern.repr <- fit.sscr(capt = test.data$bin.capt, traps = test.data$traps,
                                   mask = test.data$mask, cov.structure = "exponential",
                                   re.scale = "prob", test = "nll")
        expect_that(abs(test.bern.repr$nll - 75.40248) < 1e-4, is_true())
        ## Binomial response.
        test.binom <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "binom",
                               resp.pars = 10, cov.structure = "exponential",
                               test = "nll")
        expect_that(abs(test.binom$nll - 121.3924) < 1e-4, is_true())
    })
