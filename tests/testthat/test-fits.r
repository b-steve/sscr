context("Testing model fits")

test_that(
    "Fits with no random effects.",
    {
        compile.sscr()
        fit.pois <- fit.sscr(test.data$capt, test.data$traps, test.data$mask, "pois")
        expect_that(max(abs(fit.pois - c(2.876451, 93.388203, 30.7750916, 0.4549133))) < 1e-5, is_true())
        fit.binom <- fit.sscr(test.data$bin.capt, test.data$traps, test.data$mask)
        expect_that(max(abs(fit.binom - c(1.105443, 144.338956, 44.8268028, 0.3123132))) < 1e-5, is_true())
    })

test_that(
    "Tests with random effects.",
    {
        ## Independent random effects.
        test.ind <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                             mask = test.data$mask, resp = "pois",
                             cov.structure = "independent", test = TRUE)
        expect_that(test.ind - 124.8294 < 1e-4, is_true())
        ## Exponential covariance function.
        test.exp <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
                             mask = test.data$mask, resp = "pois",
                             cov.structure = "exponential", test = TRUE)
        expect_that(test.exp - 122.9533 < 1e-4, is_true())
        ## Full dependence (individual-level random effect).
        ##test.full <- fit.sscr(capt = test.data$capt, traps = test.data$traps,
        ##                      mask = test.data$mask, resp = "pois",
        ##                      cov.structure = "full", test = TRUE, trace = TRUE)
        ##expect_that(test.full - 0, is_true())
    })
