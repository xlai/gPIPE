g22 <- DoseGrid(levels = c(A = 2, B = 2))
g23 <- DoseGrid(levels = c(A = 2, B = 3))

bin_design <- function(g = g22, N = 12) {
  cheung_design(g, theta = 0.3, eps = 0.5, target = "MTD", outcome = "binary",
                prior_a = 0.3, prior_b = 0.7, N = N)
}

test_that("calibration_scenarios() draws up-sets of every size and restores the RNG", {
  set.seed(5)
  before <- .Random.seed
  sc <- calibration_scenarios(g23, n = 200, seed = 1)
  expect_identical(.Random.seed, before)
  expect_length(sc, 200L)
  expect_true(all(vapply(sc, function(s) is_upset(g23, s), logical(1))))
  expect_true(all(vapply(sc, length, 1L) == 6L))
  sizes <- vapply(sc, sum, 1)
  expect_true(all(0:6 %in% sizes))
  expect_identical(calibration_scenarios(g23, n = 5, seed = 2), calibration_scenarios(g23, n = 5, seed = 2))
  # The 2 x 2 grid has six up-sets; all of them turn up.
  sc22 <- unique(calibration_scenarios(g22, n = 300, seed = 3))
  expect_length(sc22, 6L)
})

test_that("calibrate_eps() returns the largest eps in the bracket with FDR <= alpha", {
  des <- bin_design()
  cal <- calibrate_eps(des, n_scenarios = 4, nsim = 30, tol = 0.1, seed = 1)
  expect_s3_class(cal, "gPIPE::PipeCalibration")
  expect_equal(cal@design@eps, cal@eps)
  expect_true(cal@fdr <= cal@alpha)
  expect_true(cal@eps >= 0.01 && cal@eps <= 0.5)
  # Every evaluation is in the trace, in order, and the chosen one is among them.
  expect_false(is.unsorted(cal@trace$eps))
  expect_true(cal@eps %in% cal@trace$eps)
  expect_equal(cal@trace$fdr[cal@trace$eps == cal@eps], cal@fdr)
  # The next bracket point above the choice failed, or the choice is the upper end.
  above <- cal@trace[cal@trace$eps > cal@eps, ]
  expect_true(nrow(above) == 0 || all(above$fdr > cal@alpha))
  # Common random numbers: same seed, same answer; and the RNG is left alone.
  set.seed(9)
  before <- .Random.seed
  again <- calibrate_eps(des, n_scenarios = 4, nsim = 30, tol = 0.1, seed = 1)
  expect_identical(.Random.seed, before)
  expect_identical(again@trace, cal@trace)
  expect_output(print(cal), "FDR")
})

test_that("calibrate_eps() handles the two ends of the bracket", {
  des <- bin_design()
  easy <- calibrate_eps(des, n_scenarios = 3, nsim = 10, alpha = 0.999, seed = 1)
  expect_equal(easy@eps, 0.5)
  expect_equal(nrow(easy@trace), 1L)
  expect_warning(
    none <- calibrate_eps(des, n_scenarios = 3, nsim = 20, alpha = 1e-6, lower = 0.45, seed = 1),
    "no `eps` is feasible"
  )
  expect_true(is.na(none@eps))
  expect_equal(none@design@eps, 0.5)
  expect_output(print(none), "No eps")
})

test_that("calibrate_eps() works for a normal MED design and with given scenarios", {
  des <- cheung_design(g22, theta = 0, eps = 0.5, sigma0 = 1, N = 12)
  expect_error(calibrate_eps(des), "band")
  sc <- list(c(0, 0, 1, 1), c(0, 1, 1, 1))
  cal <- calibrate_eps(des, scenarios = sc, band = 1.5, nsim = 20, tol = 0.2, seed = 2)
  expect_equal(cal@n_scenarios, 2L)
  expect_true(is.na(cal@eps) || cal@fdr <= cal@alpha)
})

test_that("calibrate_eps() rejects bad input", {
  des <- bin_design()
  ms <- ms_design(g22, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 4, N = 12)
  expect_error(calibrate_eps(ms), "cheung")
  expect_error(calibrate_eps(des, scenarios = list(c(1, 0, 0, 0))), "up-set")
  expect_error(calibrate_eps(des, band = 0.8), "outside")
  expect_error(calibrate_eps(des, lower = 0.5, upper = 0.2), "lower")
  expect_error(calibrate_eps(g22), "PipeDesign")
})

test_that("calibrate_eps() has a default band of 0.1 for binary outcomes, settable by the user", {
  des <- bin_design()
  expect_equal(calibrate_eps(des, n_scenarios = 2, nsim = 5, alpha = 0.999, seed = 1)@band, 0.1)
  expect_equal(calibrate_eps(des, n_scenarios = 2, nsim = 5, alpha = 0.999, band = 0.15, seed = 1)@band, 0.15)
})

test_that("calibrate_eps() keeps the truths at theta +- band for a design with a shift", {
  des <- cheung_design(DoseGrid(levels = c(A = 3, B = 3)), theta = 0.3, eps = 0.5, target = "MTD",
                       outcome = "binary", prior_a = 0.3, prior_b = 0.7, N = 12, cohort = 3,
                       shift = 0.1)
  seen <- list()
  pipe_simulate_orig <- pipe_simulate
  local_mocked_bindings(pipe_simulate = function(design, truth, ...) {
    seen[[length(seen) + 1L]] <<- list(truth = truth)
    pipe_simulate_orig(design, truth, ...)
  })
  cal <- calibrate_eps(des, scenarios = list(rep(0L, 9), rep(1L, 9)), band = 0.1, nsim = 3,
                       tol = 0.4, seed = 1)
  expect_s3_class(cal, "gPIPE::PipeCalibration")
  expect_equal(seen[[1]]$truth, rep(0.2, 9))
  expect_equal(seen[[2]]$truth, rep(0.4, 9))
  expect_equal(cal@design@shift, 0.1)
})
