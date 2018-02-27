context("Testing model fits")

test_that(
    "Fits with no random effects.",
    {
        compile.sscr()
        ## Poisson.
        fit.pois <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "pois")
        expect_that(max(abs(fit.pois - c(2.876451, 93.388203, 30.7750916, 0.4549133, -126.7372800))) < 1e-4, is_true())
        ## Poisson with detection function on probability scale.
        fit.pois.pr <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "pois", detfn.scale = "prob")
        expect_that(max(abs(fit.pois.pr - c(1.0000000, 117.5138450, 32.4247192, 0.4317694, -122.1427591))) < 1e-4, is_true())
        ## Bernoulli.
        fit.bern <- fit.sscr(test.data$bin.capt, test.data$traps, test.data$mask)
        expect_that(max(abs(fit.bern - c(1.105443, 144.338956, 44.8268028, 0.3123132, -71.5343406))) < 1e-4, is_true())
        ## Bernoulli with detection function on probability scale.
        fit.bern.pr <- fit.sscr(test.data$bin.capt, test.data$traps, test.data$mask, detfn.scale = "prob")
        expect_that(max(abs(fit.bern.pr - c(0.7097872, 157.4884412, 45.3146356, 0.3089510, -71.8004137))) < 1e-4, is_true())
        ## Binomial.
        fit.binom <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "binom", 10)
        expect_that(max(abs(fit.binom - c(0.5190142, 125.8786894, 56.4560761, 0.2479804, -134.4903162))) < 1e-4, is_true())
        ## Binomial with detection function on probability scale.
        fit.binom.pr <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "binom", 10, detfn.scale = "prob")
        expect_that(max(abs(fit.binom.pr - c(0.6049765, 137.6776331, 67.8939222, 0.2062040, -133.3041233))) < 1e-4, is_true())
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
        
        ## Testing exact vs not-exact gradients.
        fit.exact <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                              mask = test.data$mask, resp = "pois", detfn = "hn",
                              cov.structure = "exponential", test = "gr",
                              test.conditional.n = FALSE)
        fit.rexact <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "pois", detfn = "hn",
                               cov.structure = "exponential", test = "gr",
                               test.conditional.n = FALSE, Rhess = TRUE)
        expect_that(max(abs(fit.exact$gr - c(59.61226, 419.31909, 29.11101, -13.36349, 255.16827))) < 1e-4, is_true())
        expect_that(max(abs(fit.rexact$gr - c(59.61226, 419.31909, 29.11101, -13.36349, 255.16827))) < 1e-4, is_true())
        ## ... and for conditional-n models.
        fit.exact.n <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois", detfn = "hn",
                                cov.structure = "exponential", test = "gr", 
                                test.conditional.n = TRUE)
        expect_that(max(abs(fit.exact.n$gr - c(4.058224, -2.900307, 5.428064, -1.282832))) < 1e-4, is_true())
        ## Independent random effects.
        test.ind.hn <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois", detfn = "hn",
                                cov.structure = "independent", test = TRUE)
        expect_that(abs(test.ind.hn$nll - 124.8294) < 1e-4, is_true())
        ## With detection function on probability scale.
        test.ind.hn.detpr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                   mask = test.data$mask, resp = "pois", detfn = "hn",
                                   detfn.scale = "prob", cov.structure = "independent",
                                   test = TRUE)
        expect_that(abs(test.ind.hn.detpr$nll - 121.65) < 1e-4, is_true())
        ## With random effects on the probability scale.
        test.ind.hn.repr <- test.ind.hn <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                                    mask = test.data$mask, resp = "pois", detfn = "hn",
                                                    cov.structure = "independent", re.scale = "prob",
                                                    test = TRUE)
        expect_that(abs(test.ind.hn.repr$nll - 122.2674) < 1e-4, is_true())
        ## With detection function and random effects also on probability scale.
        test.ind.hn.detpr.repr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                           mask = test.data$mask, resp = "pois", detfn = "hn",
                                           detfn.scale = "prob", cov.structure = "independent",
                                           re.scale = "prob", test = TRUE)
        expect_that(abs(test.ind.hn.detpr.repr$nll - 127.2436) < 1e-4, is_true())
        ## Testing manual start values.
        test.sv.hn <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "pois", detfn = "hn",
                               cov.structure = "independent",
                               start = c(lambda0 = 2, sigma.u = 4), test = TRUE)
        expect_that(abs(test.sv.hn$nll - 147.3844) < 1e-4, is_true())
        ## ... with hazard rate detection function.
        test.ind.hr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois", detfn = "hr",
                                cov.structure = "independent", test = TRUE)
        expect_that(abs(test.ind.hr$nll - 150.5247) < 1e-4, is_true())
        ## Exponential covariance function.
        test.exp <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                             mask = test.data$mask, resp = "pois",
                             cov.structure = "exponential", test = TRUE)
        expect_that(abs(test.exp$nll - 122.9533) < 1e-4, is_true())
        ## With detection function on probability scale.
        test.exp.detpr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois",
                                cov.structure = "exponential", detfn.scale = "prob",
                                test = TRUE)
        expect_that(abs(test.exp.detpr$nll - 120.397) < 1e-4, is_true())
        ## With random effects on probability scale.
        test.exp.repr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                  mask = test.data$mask, resp = "pois",
                                  cov.structure = "exponential", re.scale = "prob",
                                  test = TRUE)
        expect_that(abs(test.exp.repr$nll - 123.3313) < 1e-4, is_true())
        ## Squared exponential covariance function (not sure we can trust this one).
        test.sqexp <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "pois",
                               cov.structure = "sq_exponential", test = TRUE)
        expect_that(abs(test.sqexp$nll - 122.4532) < 1e-4, is_true())
        ## Individual random effect.
        test.indiv <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "pois",
                               cov.structure = "individual", test = TRUE)
        expect_that(abs(test.indiv$nll - 122.8704) < 1e-4, is_true())
        ## Bernoulli response.
        test.bern <- fit.sscr(capt = test.data$bin.capt, traps = test.data$traps,
                              mask = test.data$mask, cov.structure = "exponential",
                              test = TRUE)
        expect_that(abs(test.bern$nll - 75.09304) < 1e-4, is_true())
        ## With detection function on probability scale.
        test.bern.detpr <- fit.sscr(capt = test.data$bin.capt, traps = test.data$traps,
                                 mask = test.data$mask,  cov.structure = "exponential",
                                 detfn.scale = "prob", test = TRUE)
        expect_that(abs(test.bern.detpr$nll - 73.86672) < 1e-4, is_true())
        ## With random effects on probability scale.
        test.bern.repr <- fit.sscr(capt = test.data$bin.capt, traps = test.data$traps,
                                   mask = test.data$mask, cov.structure = "exponential",
                                   re.scale = "prob", test = TRUE)
        expect_that(abs(test.bern.repr$nll - 75.40248) < 1e-4, is_true())
        ## Binomial response.
        test.binom <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "binom",
                               resp.pars = 10, cov.structure = "exponential",
                               test = TRUE)
        expect_that(abs(test.binom$nll - 121.3924) < 1e-4, is_true())
    })
