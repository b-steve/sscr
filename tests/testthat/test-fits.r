context("Testing model fits")

test_that(
    "Fits with no random effects.",
    {
        compile.sscr()
        ## Poisson.
        fit.pois <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "pois")
        expect_that(max(abs(fit.pois - c(2.876451, 93.388203, 30.7750916, 0.4549133))) < 1e-4, is_true())
        ## Bernoulli.
        fit.bern <- fit.sscr(test.data$bin.capt, test.data$traps, test.data$mask)
        expect_that(max(abs(fit.bern - c(1.105443, 144.338956, 44.8268028, 0.3123132))) < 1e-4, is_true())
        ## Binomial.
        fit.binom <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "binom", 10)
        expect_that(max(abs(fit.binom - c(0.5190142, 125.8786894, 56.4560761, 0.2479804))) < 1e-4, is_true())
    })

test_that(
    "Tests with random effects.",
    {
        ## Independent random effects.
        test.ind.hn <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois", detfn = "hn",
                                cov.structure = "independent", test = TRUE)
        expect_that(abs(test.ind.hn - 124.8294) < 1e-4, is_true())
        ## ... with hazard rate detection function.
        test.ind.hr <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                                mask = test.data$mask, resp = "pois", detfn = "hr",
                                cov.structure = "independent", test = TRUE)
        expect_that(abs(test.ind.hr - 150.5247) < 1e-4, is_true())
        ## Exponential covariance function.
        test.exp <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                             mask = test.data$mask, resp = "pois",
                             cov.structure = "exponential", test = TRUE)
        expect_that(abs(test.exp - 122.9533) < 1e-4, is_true())
        ## Squared exponential covariance function.
        test.sqexp <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                             mask = test.data$mask, resp = "pois",
                             cov.structure = "sq_exponential", test = TRUE)
        expect_that(abs(test.sqexp - 123.7825) < 1e-4, is_true()) ## FAIL
        ## Individual random effect.
        test.indiv <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                             mask = test.data$mask, resp = "pois",
                             cov.structure = "individual", test = TRUE)
        expect_that(abs(test.indiv - 122.4532) < 1e-4, is_true()) ## FAIL
        ## ... with Bernoulli response.
        test.bern <- fit.sscr(capt = test.data$bin.capt, traps = test.data$traps,
                             mask = test.data$mask,
                             cov.structure = "exponential", test = TRUE)
        expect_that(abs(test.bern - 78.60182) < 1e-4, is_true()) ## FAIL
        ## ... with binomial response.
        test.binom <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                               mask = test.data$mask, resp = "binom",
                               resp.pars = 10, cov.structure = "exponential",
                               test = TRUE)
        expect_that(abs(test.binom - 121.3924) < 1e-4, is_true())
    })
