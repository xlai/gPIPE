# Outcome coding and PIPE scores (#18, step 2).
#
# Plain functions over combination vectors, in the array order of a DoseGrid
# (drug 1 fastest). They hold no state, so the trial loop can call them with
# integer and numeric vectors only.
#
# Coding: gamma_k = 1 means pi_k > theta, and p_k = P(pi_k > theta | data).
# The gamma = 1 set is an up-set of the partial order. `target` says which
# side is favourable: "MTD" favours gamma = 0 (safe), "MED" favours gamma = 1
# (effective). eps penalises declaring the favourable side, so eps < 0.5 is
# conservative in both settings.

# ---- Beta outcome model -----------------------------------------------------

# Beta(a, b) hyperparameters with a given median and strength a + b, one per
# combination. `median` and `strength` are vectors over the combinations
# (strength is recycled); a matrix of medians is read in array order, so
# priors are always indexed by combination, never by name.
#
# method = "bisection" is pipe.design 0.5.1's `beta.med`: bisection on the
# mean m in (0, 1) until |pbeta(median, K m, K (1 - m)) - 0.5| < 1e-4. Use it
# to reproduce pipe.design trials. method = "exact" solves the median to
# machine precision with uniroot.
#
# Mander & Sweeting (2015) call strength 1 per dose "strong" and a total
# strength of 1 over the grid (1 / K per dose) "weak".
beta_prior_from_median <- function(median, strength, method = c("bisection", "exact")) {
  method <- match.arg(method)
  median <- c(median)
  if (!is.numeric(median) || anyNA(median) || any(median <= 0 | median >= 1)) {
    stop("`median` must be strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(strength) || anyNA(strength) || any(strength <= 0)) {
    stop("`strength` must be positive.", call. = FALSE)
  }
  if (length(strength) != 1L && length(strength) != length(median)) {
    stop("`strength` must have length 1 or one entry per combination.", call. = FALSE)
  }
  strength <- rep_len(c(strength), length(median))
  solve <- if (method == "bisection") beta_mean_bisection else beta_mean_exact
  m <- mapply(solve, strength, median, USE.NAMES = FALSE)
  list(a = strength * m, b = strength * (1 - m))
}

# pipe.design 0.5.1, beta.med()'s inner betaprior1(), step for step, plus a
# cap on the halvings so an extreme median or strength cannot loop forever.
beta_mean_bisection <- function(K, x, p = 0.5, tol = 1e-4, max_iter = 200L) {
  m_lo <- 0
  m_hi <- 1
  for (i in seq_len(max_iter)) {
    m0 <- (m_lo + m_hi) / 2
    p0 <- stats::pbeta(x, K * m0, K * (1 - m0))
    if (p0 < p) m_hi <- m0 else m_lo <- m0
    if (abs(p0 - p) < tol) break
  }
  m0
}

beta_mean_exact <- function(K, x, p = 0.5) {
  f <- function(m) stats::pbeta(x, K * m, K * (1 - m)) - p
  stats::uniroot(f, c(.Machine$double.eps, 1 - .Machine$double.eps), tol = 1e-14)$root
}

# P(pi_k > theta | data) under independent Beta(a, b) priors with x events in
# n patients. All arguments are vectors over the combinations (a, b recycled).
beta_exceedance <- function(theta, a, b, n, x) {
  if (length(n) != length(x)) {
    stop("`n` and `x` must have one entry per combination.", call. = FALSE)
  }
  if (any(x < 0 | x > n)) {
    stop("`x` must lie between 0 and `n`.", call. = FALSE)
  }
  stats::pbeta(theta, a + x, b + n - x, lower.tail = FALSE)
}

# ---- Scores -----------------------------------------------------------------

# Per-dose weights: "n+1" for epsilon-tapering (Cheung et al. 2022), "1" for
# the original PIPE (Mander & Sweeting 2015).
pipe_weights <- function(n, weights = c("n+1", "1")) {
  weights <- match.arg(weights)
  if (weights == "n+1") n + 1 else rep(1, length(n))
}

# Score difference d_k = gain(gamma_k = 1) - gain(gamma_k = 0), where the
# gain is w log(eps P(state)) for the favourable state and
# w log((1 - eps) P(state)) for the other one. The score of an up-set U is
# sum(d[U]) plus a constant, so the solvers only need d.
#
#   MTD: d_k = w_k [log((1 - eps) p_k) - log(eps (1 - p_k))]
#   MED: d_k = w_k [log(eps p_k) - log((1 - eps) (1 - p_k))]
#
# With eps = 0.5 both reduce to w_k qlogis(p_k). p is clipped to
# [1e-12, 1 - 1e-12] so that certain doses keep a finite score.
pipe_gain <- function(p, eps, w = 1, target = c("MTD", "MED")) {
  target <- match.arg(target)
  check_eps(eps)
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  shift <- if (target == "MTD") log((1 - eps) / eps) else log(eps / (1 - eps))
  w * (log(p) - log1p(-p) + shift)
}

# Cheung's taper by patients: eps_n = (N - n)/(N + 1) * 0.5 + (n + 1)/(N + 1) * eps,
# with n the number of patients enrolled so far. It starts near 0.5 and
# reaches eps at n = N; the final estimate uses eps itself.
taper_eps <- function(n, N, eps) {
  check_eps(eps)
  if (any(n < 0 | n > N)) {
    stop("`n` must lie between 0 and `N`.", call. = FALSE)
  }
  (N - n) / (N + 1) * 0.5 + (n + 1) / (N + 1) * eps
}

check_eps <- function(eps) {
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0 || eps >= 1) {
    stop("`eps` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
}
