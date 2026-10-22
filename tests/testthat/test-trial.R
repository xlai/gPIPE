beta_design <- function(rule = "cheung", ...) {
  g <- DoseGrid(levels = c(4, 4))
  pr <- beta_prior_from_median(rep(0.3, 16), 1)
  pipe_trial_design(g, outcome_beta(0.3, pr$a, pr$b), target = "MTD", N = 30, rule = rule, ...)
}

truth_4x4 <- c(outer(seq(0.05, 0.35, length.out = 4), seq(0, 0.3, length.out = 4), `+`))

test_that("eps-tapering trials always enrol N patients (C2)", {
  # The C2 case: prior mean at theta, so untried doses look toxic. The old loop
  # stopped after one patient in every trial.
  g <- DoseGrid(levels = c(4, 4))
  des <- pipe_trial_design(g, outcome_beta(0.3, 0.3, 0.7), "MTD", N = 30, eps = 0.25)
  sim <- pipe_simulate(des, truth_4x4, n_trials = 20, seed = 1)
  expect_true(all(sim$summary$n_enrolled == 30))
  expect_false(any(sim$summary$stopped))
  expect_equal(rowSums(sim$n), rep(30, 20))
})

test_that("the taper uses the patients enrolled so far and ends at eps (M6)", {
  used <- numeric()
  local_mocked_bindings(trial_estimate = function(design, p, n, eps) {
    used <<- c(used, eps)
    pipe_estimate(design$grid, pipe_gain(p, eps, pipe_weights(n, design$weights)), "MTD")
  })
  des <- beta_design(eps = 0.1, cohort = 3)
  set.seed(2)
  pipe_trial(des, truth_4x4)
  # after cohorts of 3 patients: n = 3, 6, ..., 27, then the final estimate
  expect_equal(used, c(taper_eps(seq(3, 27, by = 3), 30, 0.1), 0.1))
})

test_that("cohorts are cut at N and the trial is reproducible", {
  des <- beta_design(cohort = 4)
  a <- pipe_simulate(des, truth_4x4, n_trials = 3, seed = 7)
  b <- pipe_simulate(des, truth_4x4, n_trials = 3, seed = 7)
  expect_identical(a, b)
  t1 <- a$trials[[1]]
  expect_length(t1$doses, 8)
  expect_equal(sum(t1$n), 30)
  expect_equal(t1$doses[1], 1L)
  expect_equal(t1$final_set, pipe_final_set(des$grid, t1$gamma, "MTD"))
  expect_true(is.na(t1$pick) || t1$pick %in% t1$final_set)
})

test_that("Cheung doses come from S1 u S0 of the current estimate", {
  des <- beta_design(eps = 0.2, selection = "prob")
  local_mocked_bindings(draw_dose = function(x, prob) {
    expect_true(all(prob >= 0))
    x[which.max(prob)]
  })
  t1 <- pipe_trial(des, truth_4x4)
  expect_equal(t1$n_enrolled, 30)
})

test_that("M&S trials stop only through the safety rule, with no RPII", {
  des <- beta_design("ms", cohort = 2)
  sim <- pipe_simulate(des, rep(0.9, 16), n_trials = 10, seed = 3)
  expect_true(all(sim$summary$stopped))
  expect_true(all(sim$summary$n_enrolled < 30))
  expect_equal(sum(sim$final_set), 0)
  expect_true(all(is.na(sim$summary$pick)))
  safe <- pipe_simulate(des, rep(0.05, 16), n_trials = 5, seed = 4)
  expect_true(all(safe$summary$n_enrolled == 30))
  # escalation is by at most one level per drug
  for (t in safe$trials) {
    steps <- abs(diff(des$grid@combos[t$doses, , drop = FALSE]))
    expect_true(all(steps <= 1))
  }
})

test_that("normal outcomes: pooled t posterior, known sigma while df < 1", {
  n <- c(2, 1, 0, 0)
  s1 <- c(1.2, 0.9, 0, 0)
  s2 <- c(0.5^2 + 0.7^2, 0.81, 0, 0)
  p <- normal_exceedance(0.4, 1, n, s1, s2)
  sp <- sqrt(((0.5 - 0.6)^2 + (0.7 - 0.6)^2) / 1)
  expect_equal(p[1:2], stats::pt((c(0.6, 0.9) - 0.4) / (sp / sqrt(c(2, 1))), 1))
  expect_equal(p[3:4], c(0.5, 0.5))
  expect_equal(normal_exceedance(0.4, 2, c(1, 1, 0), c(1, 0, 0), c(1, 0, 0)),
               c(stats::pnorm(0.6 / 2), stats::pnorm(-0.4 / 2), 0.5))
  expect_error(outcome_normal(0.4, 0), "sigma0")
})

test_that("Cheung MED trial with normal outcomes", {
  g <- DoseGrid(levels = c(3, 3))
  truth <- list(mean = c(outer(c(0, 0.3, 0.6), c(0, 0.3, 0.6), `+`)), sd = 0.5)
  des <- pipe_trial_design(g, outcome_normal(0.5, 0.5), "MED", N = 24, eps = 0.1)
  expect_equal(des$start, 9L)
  sim <- pipe_simulate(des, truth, n_trials = 5, seed = 5)
  expect_true(all(sim$summary$n_enrolled == 24))
  expect_true(all(is.na(sim$summary$pick)))
  for (t in sim$trials) {
    expect_equal(t$final_set, pipe_final_set(g, t$gamma, "MED"))
  }
})

test_that("design checks", {
  g <- DoseGrid(levels = c(2, 2))
  m <- outcome_normal(0.5, 1)
  expect_error(pipe_trial_design(g, m, "MTD", N = 10, rule = "ms"), "binary")
  expect_error(pipe_trial_design(g, outcome_beta(0.3, 1, 1), "MED", N = 10, rule = "ms"), "MTD")
  expect_error(pipe_trial_design(g, m, "MED", N = 10, start = 5L), "start")
  expect_error(pipe_trial_design(g, outcome_beta(0.3, 1, 1), "MTD", N = 10, rule = "ms",
                                 eps = 0.2), "eps = 0.5")
  expect_error(pipe_trial_design(g, outcome_beta(0.3, 1, 1), "MTD", N = 10, rule = "ms",
                                 weights = "n+1"), "eps = 0.5")
  des <- pipe_trial_design(g, outcome_beta(0.3, 1, 1), "MTD", N = 10, rule = "ms")
  expect_equal(des$weights, "1")
  expect_equal(des$selection, "min_S")
})

test_that("default start: top dose for Cheung MED, combination 1 otherwise", {
  g <- DoseGrid(levels = c(3, 4))
  m <- outcome_normal(0.5, 1)
  b <- outcome_beta(0.3, 1, 1)
  expect_equal(pipe_trial_design(g, m, "MED", N = 10)$start, 12L)
  expect_equal(pipe_trial_design(g, b, "MTD", N = 10)$start, 1L)
  expect_equal(pipe_trial_design(g, b, "MTD", N = 10, rule = "ms")$start, 1L)
  # subset grid without the top corner: the maximal combination with the
  # largest index; without the bottom corner, the first kept combination
  keep <- matrix(TRUE, 3, 4)
  keep[3, 4] <- FALSE
  keep[1, 1] <- FALSE
  gs <- DoseGrid(levels = c(3, 4), keep = keep)
  top <- pipe_trial_design(gs, m, "MED", N = 10)$start
  expect_equal(unname(gs@combos[top, ]), c(2L, 4L))
  expect_true(top %in% maximal_elements(gs, rep(TRUE, n_combos(gs))))
  expect_equal(pipe_trial_design(gs, b, "MTD", N = 10)$start, 1L)
})
