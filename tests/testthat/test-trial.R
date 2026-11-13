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

test_that("M&S trials stop only through the safety rule, with no RP2D", {
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

test_that("an M&S design with the safety rule checks the grid when it is built", {
  g <- DoseGrid(levels = c(3, 3, 3))
  model <- outcome_beta(0.3, 0.3, 0.7)
  expect_silent(pipe_trial_design(g, model, "MTD", N = 10, rule = "ms"))
  local_mocked_bindings(upset_sweep_plan = function(grid, ...) stop("too wide", call. = FALSE))
  expect_error(pipe_trial_design(g, model, "MTD", N = 10, rule = "ms"), "too wide")
  expect_silent(pipe_trial_design(g, model, "MTD", N = 10, rule = "ms", safety = NULL))
  expect_silent(pipe_trial_design(g, model, "MTD", N = 10, rule = "cheung"))
})

startup_design <- function(grid, N = 12, ...) {
  pipe_trial_design(grid, outcome_beta(0.3, 1, 1), "MTD", N = N, startup = "diagonal", ...)
}

test_that("the diagonal start-up phase follows the diagonal and then the edge (5 x 3 grid)", {
  g <- DoseGrid(levels = c(5, 3))
  des <- startup_design(g, N = 8)
  # (1,1), (2,2), (3,3), (4,3), (5,3)
  expect_equal(des$startup_path, c(1L, 7L, 13L, 14L, 15L))
  expect_equal(unname(g@combos[des$startup_path, ]),
               cbind(1:5, c(1L, 2L, 3L, 3L, 3L)), ignore_attr = TRUE)
  t1 <- pipe_trial(des, rep(0, 15))
  expect_equal(t1$doses[1:5], c(1L, 7L, 13L, 14L, 15L))
  expect_equal(t1$n_startup, 5L)
  expect_equal(t1$n_enrolled, 8)
  # a trial that ends inside the phase counts the cohorts it treated
  expect_equal(pipe_trial(startup_design(g, N = 3), rep(0, 15))$n_startup, 3L)
})

test_that("the start-up phase ends at the first toxicity, and draws no random numbers", {
  g <- DoseGrid(levels = c(5, 3))
  des <- startup_design(g, N = 6)
  seen <- NULL
  local_mocked_bindings(draw_dose = function(x, prob) {
    seen <<- c(seen, length(x))
    x[1]
  })
  t1 <- pipe_trial(des, rep(1, 15))
  expect_equal(t1$n_startup, 1L)
  expect_equal(t1$doses[1], 1L)
  # the second cohort and all later ones come from the usual rule
  expect_length(seen, 5)
  # no toxicity: the five start-up doses draw nothing
  seen <- NULL
  pipe_trial(des, rep(0, 15))
  expect_length(seen, 1)
  set.seed(5)
  a <- .Random.seed
  startup_n <- pipe_trial(startup_design(g, N = 5), rep(0, 15))
  expect_identical(.Random.seed, a)
  expect_equal(startup_n$n_startup, 5L)
})

test_that("a toxicity after some cohorts ends the phase for good", {
  g <- DoseGrid(levels = c(5, 3))
  des <- startup_design(g, N = 10)
  # toxicity only at (3,3) = 13: the third cohort has events
  truth <- rep(0, 15)
  truth[13] <- 1
  local_mocked_bindings(draw_dose = function(x, prob) x[1])
  t1 <- pipe_trial(des, truth)
  expect_equal(t1$doses[1:3], c(1L, 7L, 13L))
  expect_equal(t1$n_startup, 3L)
  expect_equal(t1$n_enrolled, 10)
})

test_that("the phase ends where the diagonal step leaves a subset grid", {
  keep <- matrix(TRUE, 3, 3)
  keep[3, 3] <- FALSE
  keep[3, 2] <- FALSE
  g <- DoseGrid(levels = c(3, 3), keep = keep)
  des <- startup_design(g, N = 6)
  # (1,1), (2,2) are in; (3,3) is not, and (2,2) is not at the top
  expect_equal(unname(g@combos[des$startup_path, ]), cbind(1:2, 1:2), ignore_attr = TRUE)
  local_mocked_bindings(draw_dose = function(x, prob) x[1])
  t1 <- pipe_trial(des, rep(0, 7))
  expect_equal(t1$n_startup, 2L)
})

test_that("with three drugs one step raises all three", {
  g <- DoseGrid(levels = c(3, 2, 4))
  des <- startup_design(g, N = 10)
  # (1,1,1), (2,2,2), (3,2,3), (3,2,4)
  expect_equal(unname(g@combos[des$startup_path, ]),
               cbind(c(1L, 2L, 3L, 3L), c(1L, 2L, 2L, 2L), c(1L, 2L, 3L, 4L)),
               ignore_attr = TRUE)
  expect_equal(pipe_trial(des, rep(0, 24))$n_startup, 4L)
})

test_that("start-up design checks", {
  g <- DoseGrid(levels = c(3, 3))
  b <- outcome_beta(0.3, 1, 1)
  expect_error(pipe_trial_design(g, b, "MTD", N = 10, rule = "ms", startup = "diagonal"),
               "startup")
  expect_error(pipe_trial_design(g, outcome_normal(0.5, 1), "MTD", N = 10, startup = "diagonal"),
               "binary")
  expect_error(pipe_trial_design(g, b, "MED", N = 10, startup = "diagonal"), "MTD")
  expect_error(pipe_trial_design(g, b, "MTD", N = 10, startup = "up"), "none")
  expect_null(pipe_trial_design(g, b, "MTD", N = 10)$startup_path)
})

test_that("startup = 'none' reproduces the results of the version without the option", {
  g <- DoseGrid(levels = c(A = 4, B = 4))
  for (case in list(list("equal", c(171, 123, 36, 0, 126, 66, 0, 0, 60, 3, 0, 0, 15, 0, 0, 0),
                         c(1, 2, 2, 5, 5, 5, 5, 2, 5, 9), c(330, 90, 25)),
                    list("prob", c(363, 117, 27, 0, 60, 21, 0, 0, 12, 0, 0, 0, 0, 0, 0, 0),
                         c(1, 1, 1, 1, 1, 1, 1, 1, 1, 5), c(507, 61, 23)))) {
    des <- PipeDesign(grid = g, theta = 0.3, N = 30, eps = 0.25, outcome = "binary",
                      target = "MTD", rule = "cheung", selection = case[[1]], cohort = 3)
    s <- simulate(des, nsim = 20, seed = 11, truth = truth_4x4, keep_trials = TRUE)
    expect_equal(colSums(s@n), case[[2]])
    expect_equal(s@trials[[1]]$doses, case[[3]])
    expect_equal(c(sum(s@n[, 1:4]), sum(s@s1), sum(s@final_set)), case[[4]])
    expect_true(all(s@summary$n_startup == 0))
  }
})
