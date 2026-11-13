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

# ---- shift (#45) ------------------------------------------------------------

shift_design <- function(shift = 0, target = "MTD", ...) {
  if (target == "MTD") {
    cheung_design(g33, theta = 0.3, eps = 0.25, target = "MTD", outcome = "binary",
                  prior_a = 0.3, prior_b = 0.7, N = 30, cohort = 3, shift = shift, ...)
  } else {
    cheung_design(g33, theta = 0.5, eps = 0.1, target = "MED", outcome = "normal", sigma0 = 1,
                  N = 30, cohort = 3, shift = shift, ...)
  }
}
shift_truth <- list(MTD = c(0.1, 0.2, 0.3, 0.2, 0.3, 0.4, 0.3, 0.4, 0.5),
                    MED = list(mean = c(0.1, 0.3, 0.5, 0.3, 0.5, 0.7, 0.5, 0.7, 0.9), sd = 1))

test_that("shift = 0 is the default and changes nothing", {
  for (target in c("MTD", "MED")) {
    for (sel in c("equal", "prob")) {
      a <- simulate(shift_design(target = target, selection = sel), nsim = 10, seed = 5,
                    truth = shift_truth[[target]])
      b <- simulate(shift_design(0, target, selection = sel), nsim = 10, seed = 5,
                    truth = shift_truth[[target]])
      expect_identical(a@summary, b@summary)
      expect_identical(a@n, b@n)
      expect_identical(a@final_set, b@final_set)
    }
  }
  expect_equal(shift_design()@shift, 0)
})

test_that("shift = 0 reproduces results computed on dev before the shift existed", {
  # Computed on dev (commit 1dc7da0): seed 11, 20 trials, first two trials and first six picks.
  sim <- simulate(shift_design(), nsim = 20, seed = 11, truth = shift_truth$MTD)
  expect_identical(sim@n[1:2, ], matrix(c(3, 3, 0, 12, 3, 0, 0, 0, 0, 3, 9, 0, 0, 0, 6, 0, 9, 12), 2))
  expect_identical(sim@s1[1:2, ], matrix(c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 3, 5), 2))
  expect_identical(sim@summary$pick[1:6], c(8L, 5L, 2L, 3L, 2L, 4L))
})

test_that("estimate() reports P(pi > theta + shift) and next_dose() uses it", {
  dat <- data.frame(A = c(1, 1, 1, 2, 2, 3, 3), B = c(1, 1, 1, 2, 2, 3, 3),
                    y = c(0, 1, 0, 1, 1, 1, 0))
  des <- shift_design(0.1)
  est <- estimate(des, dat)
  n <- c(3, 0, 0, 0, 2, 0, 0, 0, 2)
  x <- c(1, 0, 0, 0, 2, 0, 0, 0, 1)
  expect_equal(est@p, stats::pbeta(0.4, 0.3 + x, 0.7 + n - x, lower.tail = FALSE))
  expect_false(isTRUE(all.equal(est@p, estimate(shift_design(0), dat)@p)))
  # Cell 1 has 1 event in 3 patients: posterior Beta(1.3, 2.7).
  expect_equal(est@p[1], stats::pbeta(0.4, 1.3, 2.7, lower.tail = FALSE))
  # next_dose() runs on the same threshold: its candidates come from the shifted p.
  set.seed(1)
  nd <- next_dose(des, dat)
  expect_equal(sum(nd$candidates$prob), 1)
  expect_output(print(est), "shift = 0.1")
  expect_output(print(est), "above theta \\+ shift")
  expect_no_match(capture.output(print(estimate(shift_design(0), dat))), "shift")
})

test_that("the single pick is centred on theta, not on theta + shift", {
  g1 <- DoseGrid(levels = c(A = 3))
  dat <- data.frame(A = rep(1:3, each = 20),
                    y = c(rep(1, 5), rep(0, 15), rep(1, 9), rep(0, 11), rep(1, 16), rep(0, 4)))
  mk <- function(shift) {
    cheung_design(g1, theta = 0.3, eps = 0.5, target = "MTD", outcome = "binary", prior_a = 1,
                  prior_b = 1, N = 100, shift = shift, within_set = FALSE)
  }
  # Posterior mass in [0.2, 0.4] is largest for dose 1 (5/20); in [0.4, 0.6] for dose 2 (9/20).
  mass <- function(lo, hi) {
    stats::pbeta(hi, 1 + c(5, 9), 1 + c(15, 11)) - stats::pbeta(lo, 1 + c(5, 9), 1 + c(15, 11))
  }
  expect_equal(which.max(mass(0.2, 0.4)), 1L)
  expect_equal(which.max(mass(0.4, 0.6)), 2L)
  expect_equal(estimate(mk(0.1), dat)@pick, 1L)
  # A target of 0.5 with the same delta picks dose 2, so the pick does follow theta.
  d2 <- cheung_design(g1, theta = 0.5, eps = 0.5, target = "MTD", outcome = "binary", prior_a = 1,
                      prior_b = 1, N = 100, within_set = FALSE)
  expect_equal(estimate(d2, dat)@pick, 2L)
})

test_that("summary() classifies the truth around theta, not theta + shift", {
  des <- shift_design(0.1)
  sim <- simulate(des, nsim = 8, seed = 2, truth = shift_truth$MTD)
  s <- summary(sim)
  # Truth exactly 0.3 (cells 3, 5, 7) is the target; theta + shift = 0.4 would pick cells 6 and 8.
  expect_equal(which(abs(sim@truth - 0.3) < 1e-9), c(3L, 5L, 7L))
  base <- summary(simulate(shift_design(0), nsim = 8, seed = 2, truth = shift_truth$MTD))
  expect_equal(as.character(s@by_combo$class), as.character(base@by_combo$class))
  expect_equal(s@by_combo$class[c(3, 5, 7)], base@by_combo$class[c(3, 5, 7)])
  expect_equal(as.character(s@by_combo$class)[c(1, 6, 9)], c("below", "above", "above"))
})

test_that("shift is validated", {
  expect_error(shift_design(c(0.1, 0.2)), "single finite number")
  expect_error(shift_design(NA_real_), "single finite number")
  expect_error(shift_design(Inf), "single finite number")
  expect_error(shift_design(0.7), "strictly between 0 and 1")
  expect_error(shift_design(-0.3), "strictly between 0 and 1")
  expect_equal(shift_design(-0.1)@shift, -0.1)
  expect_error(ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 30,
                         shift = 0.1), "M&S rule has no shift")
  expect_error(PipeDesign(grid = g33, theta = 0.3, N = 20, rule = "ms", eps = 0.5, shift = 0.1),
               "cheung")
  expect_s7_class(ms_design(g33, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 9, N = 30),
                  PipeDesign)
  # Normal outcomes have no bound on theta + shift.
  expect_equal(shift_design(2, "MED")@shift, 2)
  des <- shift_design(0.1)
  expect_error(des@shift <- 0.9, "strictly between")
})

update_none <- function(des) {
  PipeDesign(grid = des@grid, theta = des@theta, N = des@N, outcome = "binary", target = "MTD",
             rule = "cheung", eps = des@eps, cohort = des@cohort)
}

test_that("next_dose follows the start-up phase and then the usual rule", {
  des <- PipeDesign(grid = DoseGrid(levels = c(A = 5, B = 3)), theta = 0.3, N = 20,
                    outcome = "binary", target = "MTD", rule = "cheung", eps = 0.25,
                    startup = "diagonal")
  expect_equal(des@startup, "diagonal")
  expect_equal(next_dose(des)$dose, 1L)
  dat <- data.frame(A = c(1, 2), B = c(1, 2), y = c(0, 0))
  nd <- next_dose(des, dat)
  expect_equal(nd$dose, 13L)  # (3,3)
  expect_equal(nd$candidates$prob, 1)
  expect_false(nd$stopped)
  # the last start-up combination, then the usual rule
  dat5 <- data.frame(A = c(1, 2, 3, 4), B = c(1, 2, 3, 3), y = 0)
  expect_equal(next_dose(des, dat5)$dose, 15L)
  usual <- rbind(dat5, data.frame(A = 5, B = 3, y = 0))
  expect_gt(nrow(next_dose(des, usual)$candidates), 0)
  expect_equal(next_dose(des, usual)$candidates$prob, next_dose(update_none(des), usual)$candidates$prob)
  # a toxicity ends the phase
  tox <- data.frame(A = c(1, 2), B = c(1, 2), y = c(0, 1))
  expect_equal(next_dose(des, tox)$candidates$prob,
               next_dose(update_none(des), tox)$candidates$prob)
  # data off the path end the phase too
  off <- data.frame(A = c(1, 3), B = c(1, 1), y = c(0, 0))
  expect_equal(next_dose(des, off)$candidates$prob,
               next_dose(update_none(des), off)$candidates$prob)
  # cohorts of 2: the second patient of a cohort stays on its combination
  des2 <- PipeDesign(grid = DoseGrid(levels = c(A = 5, B = 3)), theta = 0.3, N = 20,
                     outcome = "binary", target = "MTD", rule = "cheung", eps = 0.25,
                     startup = "diagonal", cohort = 2)
  expect_equal(next_dose(des2, data.frame(A = 1, B = 1, y = 0))$dose, 1L)
  expect_equal(next_dose(des2, data.frame(A = c(1, 1), B = c(1, 1), y = 0))$dose, 7L)
  expect_equal(next_dose(des2, data.frame(A = c(1, 1, 2), B = c(1, 1, 2), y = 0))$dose, 7L)
})

test_that("PipeDesign checks the start-up option, and summary and plot work with it", {
  g <- DoseGrid(levels = c(A = 4, B = 4))
  mk <- function(...) PipeDesign(grid = g, theta = 0.3, N = 20, ...)
  expect_error(mk(startup = "diagonal", outcome = "binary", target = "MTD", rule = "ms"),
               "startup")
  expect_error(mk(startup = "diagonal", outcome = "normal", target = "MTD"), "binary")
  expect_error(mk(startup = "diagonal", outcome = "binary", target = "MED"), "MTD")
  expect_error(mk(startup = "x", outcome = "binary", target = "MTD"), "none")
  des <- mk(startup = "diagonal", outcome = "binary", target = "MTD", eps = 0.25, cohort = 2)
  truth <- c(outer(seq(0.05, 0.35, length.out = 4), seq(0, 0.3, length.out = 4), `+`))
  sim <- simulate(des, nsim = 10, seed = 3, truth = truth)
  expect_true(all(sim@summary$n_startup >= 1))
  expect_true(S7::S7_inherits(summary(sim), PipeSummary))
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(sim), "ggplot")
})
