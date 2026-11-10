# Brute-force oracle for the solvers (#31): the package does not list up-sets,
# so the tests do, straight from the cover relation of a DoseGrid.

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
