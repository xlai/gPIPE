g33 <- DoseGrid(levels = c(A = 3, B = 3))

test_that("PipeDesign validates its settings together", {
  des <- PipeDesign(grid = g33, theta = 0.3, N = 20, eps = 0.2)
  expect_true(S7::S7_inherits(des, PipeDesign))
  expect_error(PipeDesign(grid = g33, theta = 0.3, N = 20, eps = 2), "eps")
  expect_error(PipeDesign(grid = g33, theta = 0.3, N = 20, prior_a = c(1, 2)), "Beta priors")
  expect_error(PipeDesign(grid = g33, theta = 0.3, N = 20, rule = "ms", eps = 0.2), "eps = 0.5")
  expect_error(PipeDesign(grid = g33, theta = 0.3, N = 20, rule = "ms", outcome = "normal"),
               "binary")
  expect_error(PipeDesign(grid = g33, theta = 0.3, N = 20, target = "MXD"), "MTD")
  expect_error(des@eps <- 0, "eps")
})

test_that("convenience constructors set up the two designs", {
  ms <- ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 30)
  pr <- beta_prior_from_median(rep(0.3, 9), 1 / 9)
  expect_equal(ms@prior_a, pr$a)
  expect_equal(ms@rule, "ms")
  expect_equal(design_list(ms)$selection, "min_S")
  ch <- cheung_design(g33, theta = 0.5, eps = 0.05, N = 30)
  expect_equal(ch@target, "MED")
  expect_null(design_list(ch)$safety)
  expect_equal(design_list(ch)$start, 9L)
})

test_that("the constructors only add what is specific to each design", {
  # Everything else comes from the PipeDesign defaults, and `...` overrides them.
  ms <- ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 30)
  pr <- beta_prior_from_median(rep(0.3, 9), 1 / 9)
  expect_identical(ms, PipeDesign(grid = g33, theta = 0.3, target = "MTD", outcome = "binary",
                                  prior_a = pr$a, prior_b = pr$b, rule = "ms", N = 30))
  expect_equal(list(ms@cohort, ms@admis, ms@constraint, ms@safety), list(1, "closest", "neighbouring", 0.8))
  ms2 <- ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 30,
                   cohort = 3, safety = 0.9, selection = "inverse_S")
  expect_equal(list(ms2@cohort, ms2@safety, ms2@selection), list(3, 0.9, "inverse_S"))
  ch <- cheung_design(g33, theta = 0.5, eps = 0.05, N = 30)
  expect_identical(ch, PipeDesign(grid = g33, theta = 0.5, eps = 0.05, target = "MED",
                                  outcome = "normal", rule = "cheung", N = 30))
  ch2 <- cheung_design(g33, theta = 0.5, eps = 0.05, N = 30, target = "MTD", outcome = "binary",
                       selection = "prob", taper = FALSE)
  expect_equal(list(ch2@target, ch2@outcome, ch2@selection, ch2@taper),
               list("MTD", "binary", "prob", FALSE))
  expect_error(ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9), "N")
})

test_that("design_list() passes every property to the trial loop, and a safety rule only for M&S", {
  des <- PipeDesign(grid = g33, theta = 0.3, N = 20, eps = 0.2, taper = FALSE, cohort = 2,
                    weights = "1", selection = "prob", delta = 0.2, within_set = FALSE, start = 4)
  lst <- design_list(des)
  for (nm in c("target", "eps", "taper", "weights", "selection", "cohort", "delta",
               "within_set", "start", "rule", "admis", "constraint")) {
    expect_equal(lst[[nm]], S7::prop(des, nm), ignore_attr = TRUE, info = nm)
  }
  expect_identical(lst$N, 20L)
  expect_identical(lst$grid, g33)
  expect_equal(lst$model$family, "binary")
  expect_equal(lst$model$theta, 0.3)
  expect_null(lst$safety)  # the default 0.8 is ignored by the Cheung rule
  ms <- ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 20, safety = 0.9)
  expect_equal(design_list(ms)$safety, 0.9)
  expect_null(design_list(`@<-`(ms, "safety", NULL))$safety)
  nm <- PipeDesign(grid = g33, theta = 0.3, N = 20, outcome = "normal", sigma0 = 2)
  expect_equal(design_list(nm)$model$family, "normal")
})

test_that("M&S recommended doses are named rp2d, and estimate() reports them as the set", {
  ms <- ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 30)
  dat <- data.frame(A = c(1, 1, 2, 2), B = c(1, 1, 2, 2), y = c(0, 0, 0, 1))
  lst <- design_list(ms)
  x <- design_data(ms, dat)
  fin <- trial_final(lst, x$n, x$s1, x$s2, x$history)
  grid <- lst$grid
  dec <- ms_decide(grid, fin$p, x$n, S = lst$model$a + lst$model$b + x$n, history = x$history,
                   admis = lst$admis, constraint = lst$constraint, safety = lst$safety,
                   start = lst$start, final = TRUE)
  expect_true("rp2d" %in% names(dec))
  expect_false("rpII" %in% names(dec))
  expect_equal(estimate(ms, dat)@set, dec$rp2d)
})

test_that("trial data become sufficient statistics", {
  dat <- data.frame(A = c(1, 1, 2), B = c(1, 1, 3), y = c(0, 1, 2))
  x <- design_data(PipeDesign(grid = g33, theta = 0.3, N = 20, outcome = "normal"), dat)
  k <- combo_index(g33, rbind(c(1, 1), c(2, 3)))
  expect_equal(x$n[k], c(2, 1))
  expect_equal(x$s1[k], c(1, 2))
  expect_equal(x$s2[k], c(1, 4))
  expect_equal(x$history, k[c(1, 1, 2)])
  des <- PipeDesign(grid = g33, theta = 0.3, N = 20)
  expect_error(design_data(des, data.frame(A = 1, y = 0)), "columns A, B")
  expect_error(design_data(des, data.frame(A = 4, B = 1, y = 0)), "not in the grid")
  expect_error(design_data(des, data.frame(A = 1.5, B = 1, y = 0)), "whole numbers")
  expect_error(design_data(des, data.frame(A = 1, B = 1, y = 2)), "0 or 1")
  expect_error(design_data(des, data.frame(A = 1, B = 1, y = NA)), "missing")
})

test_that("an unnamed grid uses the default drug names A, B", {
  g <- DoseGrid(levels = c(3, 3))
  des <- ms_design(g, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 30)
  dat <- data.frame(A = c(1, 1), B = c(1, 1), y = c(0, 0))
  expect_true(S7::S7_inherits(estimate(des, dat), PipeEstimate))
  expect_true(next_dose(des, dat)$dose %in% seq_len(9))
  expect_error(estimate(des, data.frame(x = 1, y = 0)), "columns A, B and `y`")
})

test_that("estimate and next_dose agree with the trial loop", {
  des <- ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 30)
  dat <- data.frame(A = c(1, 1, 2, 2), B = c(1, 1, 2, 2), y = c(0, 0, 0, 1))
  x <- design_data(des, dat)
  d <- design_list(des)
  est <- estimate(des, dat)
  fin <- trial_final(d, x$n, x$s1, x$s2, x$history)
  expect_equal(est@gamma, as.integer(fin$gamma))
  expect_equal(est@q, fin$q)
  expect_output(print(est), "Recommended phase II doses")

  nd <- next_dose(des, dat)
  step <- trial_step(d, x$n, x$s1, x$s2, x$history)
  expect_equal(nd$candidates$prob, step$prob[step$prob > 0])
  expect_true(nd$dose %in% nd$candidates$dose)
  expect_equal(next_dose(des)$dose, 1L)
  full <- data.frame(A = rep(1, 30), B = rep(1, 30), y = 0)
  expect_error(next_dose(des, full), "already enrolled")
})

test_that("next_dose reports a safety stop", {
  des <- ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 30)
  dat <- data.frame(A = rep(1, 6), B = rep(1, 6), y = rep(1, 6))
  nd <- next_dose(des, dat)
  expect_true(nd$stopped)
  expect_true(is.na(nd$dose))
  expect_equal(estimate(des, dat)@set, integer())
})

test_that("simulate runs trials and keeps them only when asked", {
  des <- ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 12)
  truth <- matrix(seq(0.05, 0.6, length.out = 9), 3)
  a <- simulate(des, nsim = 4, seed = 1, truth = truth)
  b <- simulate(des, nsim = 4, seed = 1, truth = c(truth))
  expect_identical(a@summary, b@summary)
  expect_equal(dim(a@n), c(4L, 9L))
  expect_length(a@trials, 0)
  expect_length(simulate(des, nsim = 2, seed = 1, truth = truth, keep_trials = TRUE)@trials, 2)
  expect_output(print(a), "4 simulated PIPE trials")

  ch <- cheung_design(g33, theta = 0.5, eps = 0.1, N = 12)
  s <- simulate(ch, nsim = 3, seed = 2, truth = list(mean = matrix(1:9 / 9, 3), sd = 0.5))
  expect_equal(s@summary$n_enrolled, rep(12, 3))
})
