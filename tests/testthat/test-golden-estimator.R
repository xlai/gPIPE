# Golden tests (#17): properties of the PIPE estimator that any solver
# (enumeration or min-cut) must have.

# Brute-force optimum and tie rule: the largest optimal up-set for MTD, the
# smallest for MED.
brute_force_estimate <- function(levels, p, n, eps, weights, target) {
  w <- if (weights == "n+1") n + 1 else rep(1, length(p))
  ref_solve(levels, ref_gain(p, eps, w, target), largest = target == "MTD")
}

estimator_grids <- list(c(4, 4), c(3, 5), c(5, 2), c(2, 3, 3), c(3, 3, 3))

test_each_impl("estimate equals brute-force enumeration on random grids", "estimate", function(impl) {
  set.seed(17)
  for (levels in estimator_grids) {
    K <- prod(levels)
    for (rep in 1:20) {
      p <- stats::runif(K)
      n <- sample(0:6, K, replace = TRUE)
      eps <- stats::runif(1, 0.02, 0.5)
      for (target in c("MTD", "MED")) {
        for (weights in c("1", "n+1")) {
          expect_equal(impl$estimate(levels, p, n, eps, weights, target),
                       brute_force_estimate(levels, p, n, eps, weights, target),
                       info = sprintf("%s %s w=%s rep %d", paste(levels, collapse = "x"),
                                      target, weights, rep))
        }
      }
    }
  }
})

test_each_impl("exact ties go to the smallest favourable set", "estimate", function(impl) {
  for (levels in estimator_grids) {
    K <- prod(levels)
    # Every configuration scores the same: all toxic for MTD, nothing effective for MED.
    expect_equal(impl$estimate(levels, rep(0.5, K), rep(0, K), 0.5, "1", "MTD"), rep(1, K))
    expect_equal(impl$estimate(levels, rep(0.5, K), rep(0, K), 0.5, "1", "MED"), rep(0, K))
  }
  # Integer gains give many partial ties.
  set.seed(1717)
  for (levels in estimator_grids) {
    K <- prod(levels)
    for (rep in 1:20) {
      p <- stats::plogis(sample(-2:2, K, replace = TRUE))
      for (target in c("MTD", "MED")) {
        expect_equal(impl$estimate(levels, p, rep(0, K), 0.5, "1", target),
                     brute_force_estimate(levels, p, rep(0, K), 0.5, "1", target),
                     info = sprintf("%s %s rep %d", paste(levels, collapse = "x"), target, rep))
      }
    }
  }
})

test_each_impl("weights w = 1 and w = n + 1 give different estimates", "estimate", function(impl) {
  # Two doses of one drug: logit(0.6) = 0.41 at the lower dose with 9 patients,
  # logit(0.2) = -1.39 at the untried upper dose. With w = 1 the pooled gain
  # 0.41 - 1.39 < 0, so nothing exceeds; with w = n + 1 it is
  # 10 * 0.41 - 1.39 > 0, so both doses exceed.
  levels <- c(2, 1)
  p <- c(0.6, 0.2)
  n <- c(9, 0)
  for (target in c("MTD", "MED")) {
    expect_equal(impl$estimate(levels, p, n, 0.5, "1", target), c(0, 0), info = target)
    expect_equal(impl$estimate(levels, p, n, 0.5, "n+1", target), c(1, 1), info = target)
  }
})

test_each_impl("estimate does not depend on the order of the drugs", "estimate", function(impl) {
  set.seed(117)
  for (levels in estimator_grids) {
    K <- prod(levels)
    perm <- rev(seq_along(levels))
    for (rep in 1:10) {
      p <- stats::runif(K)
      n <- sample(0:6, K, replace = TRUE)
      eps <- stats::runif(1, 0.02, 0.5)
      for (target in c("MTD", "MED")) {
        swap <- function(x) c(aperm(array(x, levels), perm))
        got <- impl$estimate(levels[perm], swap(p), swap(n), eps, "n+1", target)
        expect_equal(got, swap(impl$estimate(levels, p, n, eps, "n+1", target)),
                     info = sprintf("%s %s rep %d", paste(levels, collapse = "x"), target, rep))
      }
    }
  }
})
