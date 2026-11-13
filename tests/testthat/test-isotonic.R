# Isotonic final pick (#54): the fit by repeated min-cuts against an
# enumeration of all up-sets, and the pick rule.

random_data <- function(K, seed) {
  set.seed(seed)
  n <- sample(0:6, K, replace = TRUE)
  n[sample(K, 2)] <- 0L
  n[sample(K, 2)] <- 4L
  s1 <- rbinom(K, n, sample(c(0.2, 0.5), K, replace = TRUE))
  list(n = n, s1 = s1)
}

test_that("the isotonic fit equals the enumeration reference on 2 and 3 drugs", {
  keep <- matrix(TRUE, 3, 3)
  keep[3, 1] <- keep[1, 3] <- FALSE
  specs <- list(list(levels = c(5, 3), keep = NULL), list(levels = c(3, 3), keep = NULL),
                list(levels = c(3, 3), keep = keep), list(levels = c(3, 2, 2), keep = NULL))
  for (spec in specs) {
    g <- DoseGrid(levels = spec$levels, keep = spec$keep)
    cells <- which(!is.na(c(g@lookup)))
    for (seed in 1:25) {
      d <- random_data(n_combos(g), seed)
      expect_equal(isotonic_fit(g, d$n, d$s1), ref_isotonic(spec$levels, d$n, d$s1, cells),
                   tolerance = 1e-10)
    }
  }
})

test_that("the fit is monotone, keeps NA for untried doses and equals monotone rates", {
  g <- DoseGrid(levels = c(5, 3))
  e <- g@edges
  for (seed in 30:40) {
    d <- random_data(15, seed)
    fit <- isotonic_fit(g, d$n, d$s1)
    expect_identical(is.na(fit), d$n == 0)
    # order between tried doses runs through untried ones too
    ok <- !is.na(fit)
    comp <- outer(seq_len(15), seq_len(15), function(i, j) {
      rowSums(g@combos[i, , drop = FALSE] <= g@combos[j, , drop = FALSE]) == 2L
    })
    for (i in which(ok)) for (j in which(ok)) if (comp[i, j]) expect_lte(fit[i], fit[j] + 1e-12)
  }
  n <- rep(10, 15)
  rate <- c(outer(c(.1, .2, .3, .4, .5), c(0, .1, .2), `+`))
  s1 <- rate * n
  expect_equal(isotonic_fit(g, n, s1), rate)
})

test_that("a one-drug grid is weighted PAVA", {
  g <- DoseGrid(levels = 4)
  n <- c(1, 2, 3, 1)
  s1 <- c(1, 0, 1, 1)
  expect_equal(isotonic_fit(g, n, s1), c(1, 1, 1, 1) / c(3, 3, 3, 1))
  expect_equal(isotonic_fit(g, n, s1), stats::isoreg(rep(s1 / n, n))$yf[cumsum(n)])
  # a decreasing run pools completely
  expect_equal(isotonic_fit(g, c(4, 4, 2, 2), c(3, 2, 1, 0)), rep(6 / 12, 4))
})

test_that("the isotonic pick is the candidate closest to theta", {
  g <- DoseGrid(levels = 4)
  n <- c(10, 10, 10, 10)
  s1 <- c(1, 3, 4, 7)
  expect_equal(isotonic_pick(g, 1:4, 0.3, n, s1), 2L)
  expect_equal(isotonic_pick(g, 1:4, 0.65, n, s1), 4L)
  # the target is theta, not a shifted value: 0.4 is closest to 0.45
  expect_equal(isotonic_pick(g, 1:4, 0.45, n, s1), 3L)
})

test_that("pick ties go to more patients, then to the lower index", {
  g <- DoseGrid(levels = 4)
  # estimates 0.1, 0.3, 0.4, 0.7; theta 0.35 is equally far from doses 2 and 3
  expect_equal(isotonic_pick(g, 1:4, 0.35, c(10, 10, 10, 10), c(1, 3, 4, 7)), 2L)
  expect_equal(isotonic_pick(g, 1:4, 0.35, c(10, 10, 20, 10), c(1, 3, 8, 7)), 3L)
  expect_equal(isotonic_pick(g, 1:4, 0.35, c(10, 20, 10, 10), c(1, 6, 4, 7)), 2L)
  # pooling creates equal estimates: doses 2 and 3 both at 0.5
  expect_equal(isotonic_pick(g, 1:4, 0.5, c(5, 10, 10, 5), c(0, 6, 4, 5)), 2L)
  expect_equal(isotonic_pick(g, 1:4, 0.5, c(5, 10, 12, 5), c(0, 6, 6, 5)), 3L)
})

test_that("the pick respects within_set and untried doses, and can be empty", {
  g <- DoseGrid(levels = 4)
  n <- c(10, 10, 10, 0)
  s1 <- c(1, 3, 4, 0)
  expect_equal(isotonic_pick(g, c(1, 4), 0.3, n, s1), 1L)
  expect_equal(isotonic_pick(g, c(1, 4), 0.3, n, s1, within_set = FALSE), 2L)
  expect_true(is.na(isotonic_pick(g, 4L, 0.3, n, s1)))
  expect_true(is.na(isotonic_pick(g, integer(), 0.3, n, s1)))
  expect_equal(isotonic_pick(g, integer(), 0.3, n, s1, within_set = FALSE), 2L)
  expect_true(is.na(isotonic_pick(g, 1:4, 0.3, integer(4), numeric(4))))
})

test_that("pick = 'isotonic' is only for binary MTD designs", {
  g <- DoseGrid(levels = c(3, 3))
  expect_error(pipe_trial_design(g, outcome_beta(0.3, 1, 1), "MED", N = 20, pick = "isotonic"),
               "binary outcome and MTD")
  expect_error(pipe_trial_design(g, outcome_normal(0.3, 1), "MED", N = 20, pick = "isotonic"),
               "binary outcome and MTD")
  expect_error(pipe_trial_design(g, outcome_normal(0.3, 1), "MTD", N = 20, pick = "isotonic"),
               "binary outcome and MTD")
  expect_error(pipe_trial_design(g, outcome_beta(0.3, 1, 1), "MTD", N = 20, pick = "other"))
  expect_error(PipeDesign(grid = g, target = "MED", theta = 0.3, N = 20, pick = "isotonic"),
               "binary outcome and MTD")
  expect_error(PipeDesign(grid = g, target = "MTD", theta = 0.3, outcome = "normal", N = 20,
                          pick = "isotonic"), "binary outcome and MTD")
  expect_no_error(PipeDesign(grid = g, target = "MTD", theta = 0.3, N = 20, pick = "isotonic"))
  expect_no_error(pipe_trial_design(g, outcome_beta(0.3, 1, 1), "MTD", N = 20, rule = "ms",
                                    pick = "isotonic"))
})

test_that("the default pick is unchanged from dev (values taken before #54)", {
  g <- DoseGrid(levels = c(4, 4))
  pr <- beta_prior_from_median(rep(0.3, 16), 1)
  truth <- c(outer(seq(0.05, 0.35, length.out = 4), seq(0, 0.3, length.out = 4), `+`))
  args <- list(g, outcome_beta(0.3, pr$a, pr$b), target = "MTD", N = 30, eps = 0.2, cohort = 3)
  des <- do.call(pipe_trial_design, args)
  des_i <- do.call(pipe_trial_design, c(args, pick = "interval"))
  a <- pipe_simulate(des, truth, 50, seed = 11, keep_trials = FALSE)
  b <- pipe_simulate(des_i, truth, 50, seed = 11, keep_trials = FALSE)
  expect_identical(a, b)
  expect_identical(a$summary$pick[1:20],
                   c(9L, 6L, 2L, 6L, 1L, 2L, 2L, 13L, 2L, 7L, 3L, 13L, 2L, 7L, 1L, 1L, 5L, 6L,
                     6L, 9L))
  expect_equal(sum(a$s1), 238)
  expect_equal(sum(a$final_set), 65)
  des_all <- do.call(pipe_trial_design, c(args, within_set = FALSE))
  c2 <- pipe_simulate(des_all, truth, 50, seed = 11, keep_trials = FALSE)
  expect_identical(c2$summary$pick[1:20],
                   c(2L, 1L, 6L, 1L, 5L, 3L, 5L, 6L, 3L, 2L, 1L, 3L, 6L, 2L, 2L, 5L, 9L, 2L, 1L,
                     2L))
})

test_that("the isotonic pick changes only the pick, not the trial", {
  g <- DoseGrid(levels = c(4, 4))
  pr <- beta_prior_from_median(rep(0.3, 16), 1)
  truth <- c(outer(seq(0.05, 0.35, length.out = 4), seq(0, 0.3, length.out = 4), `+`))
  args <- list(g, outcome_beta(0.3, pr$a, pr$b), target = "MTD", N = 30, eps = 0.2, cohort = 3,
               within_set = FALSE)
  a <- pipe_simulate(do.call(pipe_trial_design, args), truth, 30, seed = 3, keep_trials = FALSE)
  b <- pipe_simulate(do.call(pipe_trial_design, c(args, pick = "isotonic")), truth, 30, seed = 3,
                     keep_trials = FALSE)
  expect_identical(a$n, b$n)
  expect_identical(a$final_set, b$final_set)
  expect_true(all(b$n[cbind(seq_len(30), b$summary$pick)] > 0))
})
