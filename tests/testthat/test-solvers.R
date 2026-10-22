test_that("enumeration gives every up-set once", {
  # binomial(I + J, I) in 2D; 980 monotone 3 x 3 x 3 arrays
  for (case in list(list(c(6, 6), 924), list(c(4, 4), 70), list(c(3, 5), 56),
                    list(c(1, 4), 5), list(c(3, 3, 3), 980))) {
    g <- DoseGrid(levels = case[[1]])
    C <- grid_upsets(g)
    expect_equal(nrow(C), case[[2]])
    expect_true(all(apply(C, 1, is_upset, grid = g)))
    expect_equal(anyDuplicated(apply(C, 1, paste, collapse = "")), 0L)
  }
  expect_error(grid_upsets(DoseGrid(levels = c(4, 4)), max_configs = 10), "min-cut")
})

test_that("enumeration works on a grid with a combination left out", {
  keep <- matrix(TRUE, 3, 3)
  keep[2, 2] <- FALSE
  g <- DoseGrid(levels = c(3, 3), keep = keep)
  C <- grid_upsets(g)
  expect_true(all(apply(C, 1, is_upset, grid = g)))
  # brute force over all 2^8 subsets
  all_sets <- as.matrix(expand.grid(rep(list(0:1), 8)))
  expect_equal(nrow(C), sum(apply(all_sets, 1, is_upset, grid = g)))
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

test_that("min-cut equals enumeration on random grids, both tie rules", {
  set.seed(18)
  keep <- array(TRUE, c(3, 3, 2))
  keep[2, 2, 1] <- FALSE
  grids <- list(DoseGrid(levels = c(4, 4)), DoseGrid(levels = c(2, 5)),
                DoseGrid(levels = c(3, 3, 2)), DoseGrid(levels = c(3, 3, 2), keep = keep))
  for (g in grids) {
    C <- grid_upsets(g)
    K <- n_combos(g)
    for (rep in 1:30) {
      # half the time integer gains, which give many exact ties
      d <- if (rep %% 2) stats::rnorm(K) else sample(-2:2, K, replace = TRUE)
      for (largest in c(TRUE, FALSE)) {
        expect_equal(solve_mincut(g, d, largest), solve_enumeration(C, d, largest),
                     info = paste(paste(g@levels, collapse = "x"), rep, largest))
      }
    }
  }
})

test_that("exact ties go to the smallest favourable set", {
  for (method in c("enumerate", "mincut")) {
    g <- DoseGrid(levels = c(3, 4))
    expect_equal(pipe_estimate(g, rep(0, 12), "MTD", method), rep(1L, 12))
    expect_equal(pipe_estimate(g, rep(0, 12), "MED", method), rep(0L, 12))
  }
})

test_that("auto uses enumeration in 2D and min-cut otherwise", {
  g2 <- DoseGrid(levels = c(5, 5))
  g3 <- DoseGrid(levels = c(2, 2, 2))
  set.seed(3)
  d2 <- stats::rnorm(25)
  d3 <- stats::rnorm(8)
  expect_equal(pipe_estimate(g2, d2), pipe_estimate(g2, d2, method = "enumerate"))
  expect_equal(pipe_estimate(g3, d3), pipe_estimate(g3, d3, method = "mincut"))
  expect_error(pipe_estimate(g2, d3), "one entry per combination")
})

test_that("upset_exceedance is the M&S contour average", {
  g <- DoseGrid(levels = c(2, 3))
  C <- grid_upsets(g)
  p <- c(0.1, 0.4, 0.3, 0.6, 0.5, 0.9)
  # contour probability prod_U p prod_notU (1 - p), normalised
  prob <- apply(C, 1, function(u) prod(ifelse(u == 1, p, 1 - p)))
  prob <- prob / sum(prob)
  expect_equal(upset_exceedance(C, stats::qlogis(p)), drop(crossprod(C, prob)))
  expect_equal(upset_exceedance(C, pipe_gain(p, 0.5)), drop(crossprod(C, prob)))
})
