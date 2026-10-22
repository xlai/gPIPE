# Outcome models (#18, step 5).
#
# An outcome model is a plain list that the trial loop calls without S7
# dispatch. Data are kept as sufficient statistics per combination: n
# (patients), s1 (sum of outcomes) and s2 (sum of squared outcomes). Each
# model has
#   family:          "binary" or "normal";
#   exceedance(n, s1, s2): p = P(parameter_k > theta | data), per combination;
#   draw(k, size, truth): `size` outcomes at combination k under `truth`.
#
# `truth` is a vector of true probabilities for "binary", and a list with
# `mean` (per combination) and `sd` (scalar) for "normal".

# Independent Beta(a, b) priors with binary outcomes (Mander & Sweeting 2015).
# a and b are per combination (or scalars), e.g. from beta_prior_from_median().
outcome_beta <- function(theta, a, b) {
  list(
    family = "binary", theta = theta, a = a, b = b,
    exceedance = function(n, s1, s2) beta_exceedance(theta, a, b, n, s1),
    draw = function(k, size, truth) stats::rbinom(size, 1, truth[k])
  )
}

# Normal outcomes with a common unknown variance (Cheung et al. 2022 MED
# setting). Cheung's p_jk cannot be reproduced exactly (his sigma is not
# given), so gPIPE uses a pooled-variance t posterior under a flat prior:
#   P(mu_k > theta) is the t(df) CDF at (ybar_k - theta) / (s_p / sqrt(n_k)),
# with df = sum(n_k - 1) over tried doses. While df < 1 (or s_p = 0) the known
# `sigma0` is used with pnorm. Untried doses get p = 0.5.
outcome_normal <- function(theta, sigma0) {
  if (!is.numeric(sigma0) || length(sigma0) != 1L || !(sigma0 > 0)) {
    stop("`sigma0` must be a single positive number.", call. = FALSE)
  }
  list(
    family = "normal", theta = theta, sigma0 = sigma0,
    exceedance = function(n, s1, s2) normal_exceedance(theta, sigma0, n, s1, s2),
    draw = function(k, size, truth) stats::rnorm(size, truth$mean[k], truth$sd)
  )
}

normal_exceedance <- function(theta, sigma0, n, s1, s2) {
  p <- rep(0.5, length(n))
  tried <- n > 0
  if (!any(tried)) {
    return(p)
  }
  ybar <- s1[tried] / n[tried]
  df <- sum(n[tried] - 1)
  ss <- sum(pmax(s2[tried] - s1[tried]^2 / n[tried], 0))
  z <- if (df >= 1 && ss > 0) {
    sp <- sqrt(ss / df)
    stats::pt((ybar - theta) / (sp / sqrt(n[tried])), df)
  } else {
    stats::pnorm((ybar - theta) / (sigma0 / sqrt(n[tried])))
  }
  p[tried] <- z
  p
}
