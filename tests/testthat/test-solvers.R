# The package has one estimator (min-cut) and one routine for the M&S safety
# probabilities (a sweep). Both are checked here against a brute-force list of
# all up-sets, which only the tests build (#31).

# All up-sets of a grid, one 0/1 row each, from the cover relation.
oracle_upsets <- function(grid) {
  e <- grid@edges
  configs <- matrix(0L, nrow = 1, ncol = 0)
  for (k in seq_len(n_combos(grid))) {
    below <- e[e[, "upper"] == k, "lower"]
    forced <- if (length(below)) {
      rowSums(configs[, below, drop = FALSE]) > 0L
    } else {
      rep(FALSE, nrow(configs))
    }
    free <- configs[!forced, , drop = FALSE]
    configs <- rbind(cbind(configs, 1L), cbind(free, rep(0L, nrow(free))))
  }
  configs
}

# Best up-set by brute force, with the tie rule.
oracle_solve <- function(configs, d, largest) {
  score <- drop(configs %*% d)
  best <- which(score >= max(score) - 1e-10 * (1 + max(abs(score))))
  size <- rowSums(configs[best, , drop = FALSE])
  configs[best[if (largest) which.max(size) else which.min(size)], ]
}

# M&S q by brute force: contour probabilities prod_U p prod_notU (1 - p),
# normalised over the up-sets.
oracle_exceedance <- function(configs, p) {
  prob <- apply(configs, 1, function(u) prod(ifelse(u == 1, p, 1 - p)))
  drop(crossprod(configs, prob / sum(prob)))
}

solver_grids <- function() {
  keep2 <- matrix(TRUE, 4, 5)
  keep2[2, 3] <- FALSE
  keep2[4, 1] <- FALSE
  keep3 <- array(TRUE, c(3, 3, 2))
  keep3[2, 2, 1] <- FALSE
  list(DoseGrid(levels = 5), DoseGrid(levels = c(4, 4)), DoseGrid(levels = c(2, 5)),
       DoseGrid(levels = c(5, 3)), DoseGrid(levels = c(4, 5), keep = keep2),
       DoseGrid(levels = c(3, 3, 3)), DoseGrid(levels = c(3, 3, 2), keep = keep3),
       DoseGrid(levels = c(2, 2, 2, 2)))
}

test_that("the oracle lists every up-set once", {
  # binomial(I + J, I) in 2D; 980 monotone 3 x 3 x 3 arrays
  for (case in list(list(c(6, 6), 924), list(c(4, 4), 70), list(c(3, 5), 56),
                    list(c(1, 4), 5), list(c(3, 3, 3), 980))) {
    g <- DoseGrid(levels = case[[1]])
    C <- oracle_upsets(g)
    expect_equal(nrow(C), case[[2]])
    expect_true(all(apply(C, 1, is_upset, grid = g)))
    expect_equal(anyDuplicated(apply(C, 1, paste, collapse = "")), 0L)
  }
  # a grid with a combination left out: brute force over all 2^8 subsets
  keep <- matrix(TRUE, 3, 3)
  keep[2, 2] <- FALSE
  g <- DoseGrid(levels = c(3, 3), keep = keep)
  all_sets <- as.matrix(expand.grid(rep(list(0:1), 8)))
  expect_equal(nrow(oracle_upsets(g)), sum(apply(all_sets, 1, is_upset, grid = g)))
})

test_that("monotonic_matrices keeps its order and directions", {
  m <- monotonic_matrices(2, 2, "increasing")
  expect_length(m, 6)
  expect_equal(m[[1]], matrix(0, 2, 2))
  expect_equal(m[[2]], matrix(c(0, 0, 0, 1), 2, 2))
  expect_equal(m[[6]], matrix(1, 2, 2))
  d <- monotonic_matrices(2, 3, "decreasing")
  expect_length(d, 10)
  expect_true(all(vapply(d, function(x) all(diff(x) <= 0) && all(diff(t(x)) <= 0), TRUE)))
})

test_that("min-cut equals brute force on every kind of grid, both tie rules", {
  set.seed(18)
  for (g in solver_grids()) {
    C <- oracle_upsets(g)
    K <- n_combos(g)
    for (rep in 1:30) {
      # half the time integer gains, which give many exact ties
      d <- if (rep %% 2) stats::rnorm(K) else sample(-2:2, K, replace = TRUE)
      for (largest in c(TRUE, FALSE)) {
        expect_equal(solve_mincut(g, d, largest), oracle_solve(C, d, largest),
                     info = paste(paste(g@levels, collapse = "x"), rep, largest))
      }
    }
  }
})

test_that("exact ties go to the smallest favourable set", {
  g <- DoseGrid(levels = c(3, 4))
  expect_equal(pipe_estimate(g, rep(0, 12), "MTD"), rep(1L, 12))
  expect_equal(pipe_estimate(g, rep(0, 12), "MED"), rep(0L, 12))
})

test_that("pipe_estimate uses the same solver for any number of drugs", {
  set.seed(3)
  for (g in list(DoseGrid(levels = c(5, 5)), DoseGrid(levels = c(2, 2, 2)))) {
    d <- stats::rnorm(n_combos(g))
    expect_equal(pipe_estimate(g, d), solve_mincut(g, d, largest = TRUE))
    expect_equal(pipe_estimate(g, d, "MED"), solve_mincut(g, d, largest = FALSE))
    expect_true(is_upset(g, pipe_estimate(g, d)))
  }
  expect_error(pipe_estimate(DoseGrid(levels = c(5, 5)), stats::rnorm(8)),
               "one entry per combination")
  expect_error(pipe_estimate(DoseGrid(levels = c(2, 2)), c(1, NA, 0, 1)), "missing")
  # a grid far too large to list its up-sets
  g <- DoseGrid(levels = c(10, 10))
  expect_true(is_upset(g, pipe_estimate(g, stats::rnorm(100))))
})

test_that("upset_exceedance is the M&S contour average on every kind of grid", {
  set.seed(31)
  for (g in solver_grids()) {
    C <- oracle_upsets(g)
    K <- n_combos(g)
    for (rep in 1:10) {
      p <- stats::runif(K)
      expect_equal(upset_exceedance(g, p), oracle_exceedance(C, p), tolerance = 1e-12,
                   info = paste(paste(g@levels, collapse = "x"), rep))
    }
  }
})

test_that("upset_exceedance handles certain doses and checks its input", {
  g <- DoseGrid(levels = c(2, 3))
  q <- upset_exceedance(g, c(0, 1, 0, 1, 0, 1))
  expect_true(all(is.finite(q) & q >= 0 & q <= 1))
  # p = 0.5 everywhere: every up-set is equally likely
  C <- oracle_upsets(g)
  expect_equal(upset_exceedance(g, rep(0.5, 6)), colMeans(C))
  expect_error(upset_exceedance(g, rep(0.5, 5)), "one entry per combination")
  expect_error(upset_exceedance(g, c(0.5, NA, 0.5, 0.5, 0.5, 0.5)), "probabilities")
  expect_error(upset_exceedance(g, c(0.5, 1.2, 0.5, 0.5, 0.5, 0.5)), "probabilities")
})

test_that("the sweep needs few states and stops above its limit", {
  expect_equal(upset_sweep_plan(DoseGrid(levels = c(6, 6)))$peak, 16)
  expect_equal(upset_sweep_plan(DoseGrid(levels = c(4, 4, 4)))$peak, 372)
  # 4 x 4 x 4 has 232,848 up-sets; q is still a probability and rises with dose
  g <- DoseGrid(levels = c(4, 4, 4))
  set.seed(4)
  q <- upset_exceedance(g, stats::runif(64))
  expect_true(all(q >= 0 & q <= 1))
  expect_true(all(q[g@edges[, "lower"]] <= q[g@edges[, "upper"]] + 1e-12))
  expect_error(upset_exceedance(DoseGrid(levels = c(3, 3, 3)), rep(0.5, 27), max_states = 10),
               "max_states")
})

test_that("min-cut and the sweep reproduce the pipe.design snapshots", {
  for (setting in c("V5", "V6")) {
    states <- read_pipe_design_states(setting)
    g <- DoseGrid(levels = states[[1]]$spec$levels)
    for (s in states) {
      sp <- s$spec
      p <- beta_exceedance(sp$theta, sp$a, sp$b, sp$n, sp$r)
      expect_equal(pipe_estimate(g, pipe_gain(p, 0.5), "MTD"), s$expected$mode, info = s$id)
      q <- upset_exceedance(g, p)
      expect_equal(q, s$expected$q, tolerance = 1e-10, info = s$id)
      expect_equal(which(q < sp$safety), s$expected$safe, info = s$id)
    }
  }
})

test_that("subset grids agree with an oracle that does not use the cover edges", {
  # Up-sets straight from the definition: x in U and y >= x in every drug
  # imply y in U. All 2^K subsets are tried, so this does not rely on
  # grid@edges, which the solvers and oracle_upsets() both read.
  definition_upsets <- function(grid) {
    X <- grid@combos
    K <- nrow(X)
    above <- outer(seq_len(K), seq_len(K), Vectorize(function(i, j) all(X[j, ] >= X[i, ])))
    sets <- as.matrix(expand.grid(rep(list(0:1), K)))
    ok <- apply(sets, 1, function(u) all(u[col(above)[above]] >= u[row(above)[above]]))
    unname(sets[ok, , drop = FALSE])
  }
  keep2 <- matrix(TRUE, 3, 3)
  keep2[2, 2] <- FALSE
  keep3 <- array(TRUE, c(2, 2, 3))
  keep3[1, 2, 2] <- FALSE
  keep3[2, 1, 1] <- FALSE
  set.seed(24)
  for (g in list(DoseGrid(levels = c(3, 3), keep = keep2), DoseGrid(levels = c(2, 2, 3), keep = keep3))) {
    C <- definition_upsets(g)
    expect_equal(nrow(C), nrow(oracle_upsets(g)))
    K <- n_combos(g)
    for (rep in 1:20) {
      d <- if (rep %% 2) stats::rnorm(K) else sample(-2:2, K, replace = TRUE)
      p <- stats::runif(K)
      expect_equal(solve_mincut(g, d, TRUE), oracle_solve(C, d, TRUE))
      expect_equal(solve_mincut(g, d, FALSE), oracle_solve(C, d, FALSE))
      expect_equal(upset_exceedance(g, p), oracle_exceedance(C, p), tolerance = 1e-12)
    }
  }
})

test_that("the smallest grids work", {
  for (levels in list(1, c(1, 1), c(1, 1, 1), c(1, 4), c(3, 1, 1))) {
    g <- DoseGrid(levels = levels)
    K <- n_combos(g)
    C <- oracle_upsets(g)
    p <- seq(0.2, 0.7, length.out = K)
    expect_equal(upset_exceedance(g, p), oracle_exceedance(C, p), tolerance = 1e-12)
    d <- stats::qlogis(p)
    expect_equal(pipe_estimate(g, d), oracle_solve(C, d, TRUE))
    expect_equal(pipe_estimate(g, d, "MED"), oracle_solve(C, d, FALSE))
  }
  # one combination: q is p itself
  expect_equal(upset_exceedance(DoseGrid(levels = 1), 0.3), 0.3)
})

test_that("the sweep does not underflow when p is 0 or 1", {
  # One drug, 60 levels: surely toxic at the 30 lowest doses, surely safe at
  # the 30 highest. The 61 up-sets are the top s doses, s = 0..60. A sweep
  # that drops states far below the best one returns q = 0 or 1 here.
  g <- DoseGrid(levels = 60)
  p <- rep(c(1, 0), each = 30)
  C <- 1L * outer(0:60, 1:60, function(s, k) k > 60 - s)
  score <- drop(C %*% stats::qlogis(pmin(pmax(p, 1e-12), 1 - 1e-12)))
  w <- exp(score - max(score))
  expect_equal(upset_exceedance(g, p), drop(crossprod(C, w / sum(w))), tolerance = 1e-10)
  # the same contradiction on three drugs
  g3 <- DoseGrid(levels = c(3, 3, 3))
  C3 <- oracle_upsets(g3)
  for (p3 in list(as.numeric(rowSums(g3@combos) <= 5), rep(c(1, 0), length.out = 27))) {
    score <- drop(C3 %*% stats::qlogis(pmin(pmax(p3, 1e-12), 1 - 1e-12)))
    w <- exp(score - max(score))
    expect_equal(upset_exceedance(g3, p3), drop(crossprod(C3, w / sum(w))), tolerance = 1e-10)
  }
  expect_equal(unname(log_sum_exp_by(c(0, -2000, -2000, -Inf, -Inf), c(1, 2, 2, 3, 3))),
               c(0, -2000 + log(2), -Inf))
})

test_that("the tie tolerance follows the largest gain, not the total", {
  # 40 doses of one drug. Adding the lowest dose costs 1e-7, a real loss next
  # to gains of 100, so it stays out even though the gains sum to 3900.
  g <- DoseGrid(levels = 40)
  gamma <- solve_mincut(g, c(-1e-7, rep(100, 39)), largest = TRUE)
  expect_equal(gamma, c(0L, rep(1L, 39)))
  # weights n + 1 on a larger grid, against brute force
  g2 <- DoseGrid(levels = c(5, 5))
  C <- oracle_upsets(g2)
  set.seed(51)
  for (rep in 1:20) {
    n <- sample(0:30, 25, replace = TRUE)
    d <- pipe_gain(stats::runif(25), 0.3, n + 1, "MTD")
    expect_equal(solve_mincut(g2, d, TRUE), oracle_solve(C, d, TRUE))
  }
})

test_that("sweep plans are cached per grid and respect a lower state limit", {
  full <- DoseGrid(levels = c(3, 3, 2))
  keep <- array(TRUE, c(3, 3, 2))
  keep[2, 2, 1] <- FALSE
  part <- DoseGrid(levels = c(3, 3, 2), keep = keep)
  # same levels, different combinations: different plans
  expect_length(upset_sweep_plan(full)$steps, 18)
  expect_length(upset_sweep_plan(part)$steps, 17)
  expect_identical(upset_sweep_plan(full), upset_sweep_plan(full))
  # a cached plan above a new, lower limit is refused
  expect_error(upset_sweep_plan(full, max_states = 5), "max_states")
  expect_length(upset_sweep_plan(full)$steps, 18)
})
