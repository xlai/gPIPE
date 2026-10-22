# Dose-finding rules (#18, step 4).
#
# Plain functions over combination vectors of a DoseGrid, with gamma = 1 for
# doses whose estimate exceeds theta (toxic for MTD, effective for MED). Two
# designs use them:
#
# - Cheung et al. (2022), epsilon-tapering: the sampling set S1 u S0;
# - Mander & Sweeting (2015), original PIPE: admissible set (constraint and
#   safety), closest or adjacent dominant doses, dose choice by sample size,
#   and the recommended phase II doses (RPII), as pipe.design 0.5.1 runs them.

# ---- Cheung (2022) ----------------------------------------------------------

# Sampling set of Cheung's Algorithm 1, step 3. S1: doses in U with a lower
# cover outside U, or with no lower cover at all (the bottom dose). S0: doses
# outside U with an upper cover in U. When U is empty, the top dose(s).
cheung_candidates <- function(grid, gamma) {
  e <- grid@edges
  K <- n_combos(grid)
  if (!any(gamma == 1)) {
    return(setdiff(seq_len(K), e[, "lower"]))
  }
  in_u <- gamma == 1
  bottom <- !seq_len(K) %in% e[, "upper"]
  s1 <- in_u & (bottom | seq_len(K) %in% e[!in_u[e[, "lower"]], "upper"])
  s0 <- !in_u & seq_len(K) %in% e[in_u[e[, "upper"]], "lower"]
  which(s1 | s0)
}

# ---- Final estimates --------------------------------------------------------

# Dose set reported at the end of an epsilon-tapering trial: for MTD the
# maximal safe doses, for MED the minimal effective doses.
pipe_final_set <- function(grid, gamma, target = c("MTD", "MED")) {
  target <- match.arg(target)
  if (target == "MTD") maximal_elements(grid, gamma == 0) else minimal_elements(grid, gamma == 1)
}

# Single dose from a set (Riviere et al. 2014): among the tried doses, the one
# with the largest posterior probability that pi lies in [theta - delta,
# theta + delta]. `within_set = FALSE` searches all tried doses instead.
# Returns NA when no dose qualifies.
riviere_pick <- function(set, theta, a, b, n, x, delta = 0.1, within_set = TRUE) {
  tried <- which(n > 0)
  pool <- if (within_set) intersect(set, tried) else tried
  if (!length(pool)) {
    return(NA_integer_)
  }
  a <- rep_len(a, length(n))[pool] + x[pool]
  b <- rep_len(b, length(n))[pool] + n[pool] - x[pool]
  prob <- stats::pbeta(theta + delta, a, b) - stats::pbeta(theta - delta, a, b)
  pool[which.max(prob)]
}

# ---- Mander & Sweeting (2015) -----------------------------------------------

# Doses allowed by the constraint, given the doses given so far (`history`,
# one combination index per cohort, last = current dose). Before the first
# cohort only `start` is allowed.
#   "neighbouring": Chebyshev distance <= 1 from the last dose;
#   "non-neighbouring": at most one level above some administered dose in
#   every drug (M&S Sec. 2.4.2).
ms_constraint <- function(grid, history, constraint = c("neighbouring", "non-neighbouring"),
                          start = 1L) {
  constraint <- match.arg(constraint)
  if (!length(history)) {
    return(as.integer(start))
  }
  if (constraint == "neighbouring") {
    return(grid_neighbourhood(grid, history[length(history)]))
  }
  given <- grid@combos[unique(history), , drop = FALSE]
  ok <- vapply(seq_len(n_combos(grid)), function(k) {
    any(colSums(t(given) + 1L >= grid@combos[k, ]) == ncol(given))
  }, logical(1))
  which(ok)
}

# Closest rule within the admissible set A (logical over combinations): doses
# outside A count as neither safe nor toxic. A toxic dose in A is dominant if
# no lower cover is toxic and in A; a safe dose in A if no upper cover is safe
# and in A.
ms_closest <- function(grid, gamma, A) {
  e <- grid@edges
  tox <- A & gamma == 1
  safe <- A & gamma == 0
  blocked_tox <- e[tox[e[, "upper"]] & tox[e[, "lower"]], "upper"]
  blocked_safe <- e[safe[e[, "lower"]] & safe[e[, "upper"]], "lower"]
  sort(c(setdiff(which(tox), blocked_tox), setdiff(which(safe), blocked_safe)))
}

# Adjacent rule on the full contour. A toxic dose is adjacent if a lower
# Chebyshev-1 neighbour is safe, and a safe dose if an upper one is toxic.
# Off-grid counts as safe below and toxic above; combinations left out of the
# grid count as neither.
ms_adjacent <- function(grid, gamma) {
  D <- length(grid@levels)
  shifts <- as.matrix(expand.grid(rep(list(0:1), D)))[-1, , drop = FALSE]
  low_safe <- up_tox <- logical(n_combos(grid))
  state_at <- function(pos) {
    off <- rowSums(pos < 1L | sweep(pos, 2, grid@levels, `>`)) > 0
    idx <- rep(NA_integer_, nrow(pos))
    idx[!off] <- grid@lookup[pos[!off, , drop = FALSE]]
    list(off = off, gamma = gamma[idx])
  }
  for (s in seq_len(nrow(shifts))) {
    below <- state_at(sweep(grid@combos, 2, shifts[s, ]))
    above <- state_at(sweep(grid@combos, 2, shifts[s, ], `+`))
    low_safe <- low_safe | below$off | (!is.na(below$gamma) & below$gamma == 0)
    up_tox <- up_tox | above$off | (!is.na(above$gamma) & above$gamma == 1)
  }
  which((gamma == 1 & low_safe) | (gamma == 0 & up_tox))
}

# Dose-choice weights over the candidates, from S = a + b + n:
#   "min_S": equal weight on the candidates with the smallest S, where S
#     values within a relative 1e-9 count as equal (a + b differ by rounding);
#   "inverse_S": weighted randomisation proportional to 1/S;
#   "equal": equal randomisation;
#   "inverse_sqrt_S": proportional to 1/sqrt(S), kept for back-compatibility.
ms_selection_weights <- function(candidates, S,
                                 selection = c("min_S", "inverse_S", "equal", "inverse_sqrt_S")) {
  selection <- match.arg(selection)
  prob <- numeric(length(S))
  s <- S[candidates]
  prob[candidates] <- switch(selection,
    min_S = as.numeric(s <= min(s) * (1 + 1e-9)),
    inverse_S = 1 / s,
    equal = 1,
    inverse_sqrt_S = 1 / sqrt(s)
  )
  prob / sum(prob)
}

# One M&S decision on fixed data. p = P(pi > theta | data), n patients per
# dose, S = a + b + n, history = doses given so far (one per cohort). With
# `safety = NULL` the safety rule is off and no q is computed.
#
# Returns p, the mode contour, q, the safe, constraint and admissible sets,
# and whether the trial stops (no admissible dose; pipe.design gives no RPII
# then). Unless stopped: the dominant and candidate doses, the min-S ties and
# the dose-choice weights for `selection`. With final = TRUE, the RPII:
# closest doses over the whole grid that are safe, tried and pass the safety
# rule.
ms_decide <- function(grid, p, n, S, history, admis = c("closest", "adjacent"),
                      constraint = "neighbouring", safety = 0.8, selection = "inverse_S",
                      start = 1L, final = FALSE) {
  admis <- match.arg(admis)
  K <- n_combos(grid)
  d <- pipe_gain(p, 0.5)
  if (is.null(safety)) {
    gamma <- pipe_estimate(grid, d, "MTD")
    q <- NULL
    safe <- rep(TRUE, K)
  } else {
    configs <- grid_upsets(grid)
    gamma <- solve_enumeration(configs, d, largest = TRUE)
    q <- upset_exceedance(configs, p)
    safe <- q < safety
  }
  allowed <- seq_len(K) %in% ms_constraint(grid, history, constraint, start)
  A <- allowed & safe
  out <- list(p = p, mode = gamma, q = q, safe = which(safe), neighbourhood = which(allowed),
              admissible = which(A), stopped = !any(A))
  if (!out$stopped) {
    dominant <- if (admis == "adjacent") ms_adjacent(grid, gamma) else integer()
    if (!length(intersect(dominant, which(A)))) {
      dominant <- ms_closest(grid, gamma, A)
    }
    cand <- sort(intersect(dominant, which(A)))
    out$dominant <- dominant
    out$candidates <- cand
    out$min_S <- which(ms_selection_weights(cand, S, "min_S") > 0)
    out$prob <- ms_selection_weights(cand, S, selection)
  }
  if (final) {
    out$rpII <- if (out$stopped) {
      integer()
    } else {
      closest <- ms_closest(grid, gamma, rep(TRUE, K))
      sort(closest[gamma[closest] == 0 & n[closest] > 0 & safe[closest]])
    }
  }
  out
}
