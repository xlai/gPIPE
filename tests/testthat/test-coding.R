test_that("bisection prior reproduces pipe.design's beta.med (M&S V6 prior)", {
  # Scenario-A medians of Mander & Sweeting (2015) Sec. 3.4 on a 4 x 4 grid,
  # weak prior a + b = 1/16, in array order. `a` is what pipe.design 0.5.1
  # returns (validation/pipe_design_reference/ms2015_reference_spec.json).
  med <- matrix(c(0.04, 0.08, 0.12, 0.16, 0.1, 0.14, 0.18, 0.22,
                  0.16, 0.2, 0.24, 0.28, 0.22, 0.26, 0.3, 0.34), 4, 4)
  a <- c(0.028350830078125, 0.0290069580078125, 0.0294036865234375, 0.0297088623046875,
         0.0292205810546875, 0.029571533203125, 0.02984619140625, 0.0300750732421875,
         0.0297088623046875, 0.0299606323242188, 0.030181884765625, 0.030364990234375,
         0.0300750732421875, 0.0302734375, 0.03045654296875, 0.0306243896484375)
  pr <- beta_prior_from_median(med, 1 / 16)
  expect_equal(pr$a, a, tolerance = 1e-14)
  expect_equal(pr$a + pr$b, rep(1 / 16, 16))
})

test_that("bisection leaves a theta-median dose just on the safe side, as pipe.design does", {
  # V5: weak prior a + b = 1/36 with median 0.3 = theta. pipe.design gives
  # P(pi <= theta) - 0.5 = +9.3e-5, which decides the prior mode contour.
  pr <- beta_prior_from_median(0.3, 1 / 36)
  p <- beta_exceedance(0.3, pr$a, pr$b, 0, 0)
  expect_equal(0.5 - p, 9.27e-5, tolerance = 1e-3)
  expect_lt(p, 0.5)
})

test_that("exact prior hits the median and recycles the strength", {
  med <- c(0.05, 0.2, 0.5, 0.9)
  pr <- beta_prior_from_median(med, 1, method = "exact")
  expect_equal(stats::qbeta(0.5, pr$a, pr$b), med, tolerance = 1e-10)
  expect_equal(pr$a + pr$b, rep(1, 4))
  pr2 <- beta_prior_from_median(med, c(1, 2, 3, 4), method = "exact")
  expect_equal(pr2$a + pr2$b, c(1, 2, 3, 4))
  expect_error(beta_prior_from_median(med, c(1, 2)), "one entry per combination")
  expect_error(beta_prior_from_median(c(0, 0.5), 1), "strictly between")
  expect_error(beta_prior_from_median(0.5, 0), "positive")
})

test_that("per-combination priors give no NA posteriors (C1)", {
  # The old model looked priors up by dose name, so an unnamed vector of
  # priors gave NA. Here priors are indexed by combination.
  pr <- beta_prior_from_median(c(0.1, 0.2, 0.3, 0.4), 1)
  p <- beta_exceedance(0.3, pr$a, pr$b, n = c(3, 3, 0, 0), x = c(0, 1, 0, 0))
  expect_false(anyNA(p))
  expect_equal(p, 1 - stats::pbeta(0.3, pr$a + c(0, 1, 0, 0), pr$b + c(3, 2, 0, 0)))
  expect_error(beta_exceedance(0.3, 1, 1, n = 2, x = 3), "between 0 and `n`")
  expect_error(beta_exceedance(0.3, 1, 1, n = c(1, 2), x = 0), "one entry per combination")
})

test_that("gains follow the MTD and MED formulas", {
  p <- c(0.01, 0.3, 0.5, 0.8)
  w <- c(1, 2, 3, 4)
  eps <- 0.1
  expect_equal(pipe_gain(p, eps, w, "MTD"),
               w * (log((1 - eps) * p) - log(eps * (1 - p))))
  expect_equal(pipe_gain(p, eps, w, "MED"),
               w * (log(eps * p) - log((1 - eps) * (1 - p))))
  # eps = 0.5 is the original PIPE: d = w logit(p) in both settings
  expect_equal(pipe_gain(p, 0.5, w, "MTD"), w * stats::qlogis(p))
  expect_equal(pipe_gain(p, 0.5, w, "MED"), w * stats::qlogis(p))
  # eps < 0.5 makes the favourable call harder in both settings
  expect_true(all(pipe_gain(p, eps, 1, "MTD") > stats::qlogis(p)))
  expect_true(all(pipe_gain(p, eps, 1, "MED") < stats::qlogis(p)))
})

test_that("gains stay finite at p = 0 and p = 1", {
  d <- pipe_gain(c(0, 1), 0.5)
  expect_true(all(is.finite(d)))
  expect_equal(d, c(-1, 1) * stats::qlogis(1 - 1e-12), tolerance = 1e-6)
  expect_error(pipe_gain(0.5, 0), "eps")
  expect_error(pipe_gain(0.5, c(0.1, 0.2)), "eps")
})

test_that("weights are n + 1 for eps-tapering and 1 for the original PIPE", {
  n <- c(0, 2, 5)
  expect_equal(pipe_weights(n), c(1, 3, 6))
  expect_equal(pipe_weights(n, "1"), c(1, 1, 1))
})

test_that("taper by patients matches Cheung et al. (2022)", {
  # N = 125, eps = 0.05: 0.493 after the first subject (Table 1), eps at N
  expect_equal(taper_eps(1, 125, 0.05), 0.4928571, tolerance = 1e-6)
  expect_equal(taper_eps(125, 125, 0.05), 0.05)
  expect_equal(taper_eps(0, 125, 0.05), 125 / 126 * 0.5 + 1 / 126 * 0.05)
  expect_true(all(diff(taper_eps(0:125, 125, 0.05)) < 0))
  expect_error(taper_eps(126, 125, 0.05), "between 0 and `N`")
})
