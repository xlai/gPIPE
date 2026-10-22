# Reference implementation for the golden tests (#17).
#
# A deliberately slow, brute-force version of the PIPE rules, written only for
# the tests. It is checked against pipe.design 0.5.1 and Cheung et al. (2022)
# Table 1 in test-golden-reference.R, and the package is then checked against
# it (and against the same fixtures) through the adapters in helper-impl.R.
#
# Conventions (core design notes, #18):
# - a grid is given by `levels`, e.g. c(6, 6); combinations are ordered as in
#   R arrays, drug 1 fastest, so a 2D vector is c(m) with m[doseA, doseB];
# - gamma = 1 means pi_k > theta, and p_k = P(pi_k > theta | data);
# - the gamma = 1 set U is an up-set of the product order.

ref_combos <- function(levels) {
  combos <- as.matrix(expand.grid(lapply(levels, seq_len), KEEP.OUT.ATTRS = FALSE))
  unname(combos)
}

# Cover relation (unit steps) as a two-column matrix: lower, upper.
ref_covers <- function(levels) {
  combos <- ref_combos(levels)
  lookup <- array(seq_len(nrow(combos)), dim = levels)
  edges <- lapply(seq_along(levels), function(dim) {
    lower <- which(combos[, dim] < levels[dim])
    up <- combos[lower, , drop = FALSE]
    up[, dim] <- up[, dim] + 1L
    cbind(lower = lower, upper = lookup[up])
  })
  do.call(rbind, edges)
}

# All up-sets as a 0/1 matrix, one row per configuration. Lower covers come
# first in array order, so each cell only depends on cells already filled.
ref_upsets <- function(levels) {
  covers <- ref_covers(levels)
  configs <- matrix(0L, nrow = 1, ncol = 0)
  for (k in seq_len(prod(levels))) {
    below <- covers[covers[, "upper"] == k, "lower"]
    forced <- if (length(below)) apply(configs[, below, drop = FALSE], 1, max) == 1L else
      rep(FALSE, nrow(configs))
    free <- configs[!forced, , drop = FALSE]
    configs <- rbind(cbind(configs, 1L), cbind(free, rep(0L, nrow(free))))
  }
  configs
}

# Score difference d_k = gain(gamma_k = 1) - gain(gamma_k = 0).
ref_gain <- function(p, eps = 0.5, w = 1, target = c("MTD", "MED")) {
  target <- match.arg(target)
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  if (target == "MTD") {
    w * (log((1 - eps) * p) - log(eps * (1 - p)))
  } else {
    w * (log(eps * p) - log((1 - eps) * (1 - p)))
  }
}

# Best up-set for the score U -> sum_{k in U} d_k. Ties go to the smallest
# favourable set: the largest U for MTD, the smallest U for MED. The largest
# (smallest) optimal up-set is the union (intersection) of all optimal ones, so
# it is unique.
ref_solve <- function(levels, d, largest = TRUE, configs = ref_upsets(levels)) {
  score <- drop(configs %*% d)
  tol <- 1e-10 * (1 + max(abs(score)))
  best <- configs[score >= max(score) - tol, , drop = FALSE]
  size <- rowSums(best)
  pick <- if (largest) which(size == max(size)) else which(size == min(size))
  stopifnot(length(pick) == 1)
  best[pick, ]
}

ref_estimate <- function(levels, p, n = 0, eps = 0.5, weights = c("1", "n+1"),
                         target = c("MTD", "MED")) {
  weights <- match.arg(weights)
  target <- match.arg(target)
  w <- if (weights == "n+1") n + 1 else rep(1, length(p))
  ref_solve(levels, ref_gain(p, eps, w, target), largest = target == "MTD")
}

ref_minimal <- function(levels, inside) {
  e <- ref_covers(levels)
  setdiff(which(inside), e[inside[e[, "lower"]] & inside[e[, "upper"]], "upper"])
}

# Cheung et al. (2022) sampling grid, Algorithm 1 step 3: S1 = doses in U
# with a lower cover outside U, S0 = doses outside U with an upper cover in U,
# and the top dose when U is empty. Off-grid covers are ignored, so the bottom
# dose is in S1 whenever it is in U. Table 1 does not pin the off-grid case.
ref_cheung_candidates <- function(levels, gamma) {
  if (!any(gamma == 1)) return(prod(levels))
  e <- ref_covers(levels)
  s1 <- gamma == 1 & (seq_along(gamma) == 1 | seq_along(gamma) %in% e[gamma[e[, "lower"]] == 0, "upper"])
  s0 <- gamma == 0 & seq_along(gamma) %in% e[gamma[e[, "upper"]] == 1, "lower"]
  which(s1 | s0)
}

# Cheung's taper by patients: eps_n = (N - n)/(N + 1) * 0.5 + (n + 1)/(N + 1) * eps.
ref_taper <- function(n, N, eps) {
  (N - n) / (N + 1) * 0.5 + (n + 1) / (N + 1) * eps
}

# ---- Mander & Sweeting (2015) rules, as pipe.design 0.5.1 runs them ----

# Chebyshev distance <= 1 from combination k.
ref_box <- function(levels, k) {
  combos <- ref_combos(levels)
  which(apply(abs(sweep(combos, 2, combos[k, ])) <= 1, 1, all))
}

# Closest rule inside the admissible set A: doses outside A count as neither
# safe nor toxic.
ref_closest <- function(levels, gamma, A) {
  e <- ref_covers(levels)
  tox <- A & gamma == 1
  safe <- A & gamma == 0
  blocked_tox <- e[tox[e[, "upper"]] & tox[e[, "lower"]], "upper"]
  blocked_safe <- e[safe[e[, "lower"]] & safe[e[, "upper"]], "lower"]
  sort(c(setdiff(which(tox), blocked_tox), setdiff(which(safe), blocked_safe)))
}

# Adjacent rule on the full contour (2D only, as in pipe.design): off-grid
# counts as safe below and toxic above.
ref_adjacent <- function(levels, gamma) {
  m <- matrix(gamma, levels[1], levels[2])
  below <- rbind(0L, cbind(0L, m))
  above <- rbind(cbind(m, 1L), 1L)
  I <- levels[1]
  J <- levels[2]
  low_safe <- below[1:I, 2:(J + 1)] == 0 | below[2:(I + 1), 1:J] == 0 | below[1:I, 1:J] == 0
  up_tox <- above[2:(I + 1), 1:J] == 1 | above[1:I, 2:(J + 1)] == 1 | above[2:(I + 1), 2:(J + 1)] == 1
  which(c((m == 1 & low_safe) | (m == 0 & up_tox)))
}

# One M&S decision. `spec` holds levels, theta, a, b, n, r, history (one
# combination index per cohort, last = current dose), admis, final, and the
# safety threshold. Doses are combination indices.
ref_ms_state <- function(spec) {
  levels <- spec$levels
  K <- prod(levels)
  p <- 1 - stats::pbeta(spec$theta, spec$a + spec$r, spec$b + spec$n - spec$r)
  d <- ref_gain(p)
  configs <- ref_upsets(levels)
  score <- drop(configs %*% d)
  wts <- exp(score - max(score))
  wts <- wts / sum(wts)
  q <- drop(crossprod(configs, wts))
  gamma <- ref_solve(levels, d, largest = TRUE, configs = configs)

  safe <- q < spec$safety
  nbhd <- seq_len(K) %in% ref_box(levels, spec$history[length(spec$history)])
  A <- nbhd & safe
  out <- list(p = p, mode = gamma, q = q, safe = which(safe), neighbourhood = which(nbhd),
              admissible = which(A), stopped = !any(A))
  if (!out$stopped) {
    dominant <- if (spec$admis == "adjacent") ref_adjacent(levels, gamma) else integer()
    if (!length(intersect(dominant, which(A)))) dominant <- ref_closest(levels, gamma, A)
    cand <- intersect(dominant, which(A))
    S <- spec$a + spec$b + spec$n
    prob <- numeric(K)
    prob[cand] <- 1 / S[cand]
    out$dominant <- dominant
    out$candidates <- sort(cand)
    # Prior a + b is the same at every dose up to rounding, so ties in S are
    # ties in n; compare with a tolerance as pipe.design's 1/S does in effect.
    out$min_S <- sort(cand[S[cand] <= min(S[cand]) * (1 + 1e-9)])
    out$prob <- prob / sum(prob)
  }
  if (isTRUE(spec$final)) {
    out$rpII <- integer()
    if (!out$stopped) {
      closest <- ref_closest(levels, gamma, rep(TRUE, K))
      out$rpII <- sort(closest[gamma[closest] == 0 & spec$n[closest] > 0 & safe[closest]])
    }
  }
  out
}
