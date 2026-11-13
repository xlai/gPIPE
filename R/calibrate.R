#' Calibration scenarios
#'
#' Draws monotone true configurations for [calibrate_eps()]. Each is a 0/1
#' vector over the combinations, 1 = above `theta`, and is an up-set of the
#' grid (the parameter never falls when a dose goes up). Every up-set can
#' occur, including all-0 and all-1; a configuration is built by adding
#' random combinations whose upper covers are already in it, so the draw is
#' not uniform over up-sets but covers all sizes.
#'
#' @param grid A [DoseGrid].
#' @param n Number of configurations.
#' @param seed Optional seed; the random-number state is restored afterwards.
#' @return A list of `n` integer vectors, each of length `n_combos(grid)`.
#' @examples
#' calibration_scenarios(DoseGrid(levels = c(A = 3, B = 3)), n = 3, seed = 1)
#' @export
calibration_scenarios <- function(grid, n = 20, seed = NULL) {
  K <- n_combos(grid)
  with_seed(seed, lapply(seq_len(n), function(i) {
    size <- sample.int(K + 1L, 1L) - 1L
    inside <- logical(K)
    ups <- lapply(seq_len(K), function(k) upper_covers(grid, k))
    for (step in seq_len(size)) {
      ok <- which(!inside & vapply(ups, function(u) all(inside[u]), logical(1)))
      inside[if (length(ok) == 1L) ok else sample(ok, 1L)] <- TRUE
    }
    as.integer(inside)
  }))
}

# Put the random-number state back when the calling function exits.
restore_rng <- function(env = parent.frame()) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  do.call(on.exit, list(bquote(
    if (is.null(.(old))) {
      if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", .(old), globalenv())
    }
  ), add = TRUE), envir = env)
}

# Run `expr` under `seed` (when given); the random-number state is restored.
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  restore_rng()
  set.seed(seed)
  expr
}

#' Calibrate epsilon
#'
#' Chooses the largest asymmetry `eps` for which the false discovery rate
#' stays at or below `alpha` over an ensemble of true scenarios (Cheung et al.
#' 2022). The false discovery rate is Cheung's eq. 7: the share of the final
#' dose set (minimal `gamma` = 1 combinations for MED, maximal `gamma` = 0
#' combinations for MTD) whose true `gamma` is wrong for the target, that is
#' 0 for MED and 1 for MTD. The true-positive rate (eq. 8) at that `eps` is
#' reported alongside. Only the true `gamma` of each scenario matters; there
#' is no separate notion of safe or toxic.
#'
#' Each scenario is a monotone configuration of true `gamma` (see
#' [calibration_scenarios()]). To generate data, combinations with `gamma` = 1
#' are placed at `theta + band` and the rest at `theta - band` for binary
#' outcomes. For normal outcomes the rest sit at `theta`, the least favourable
#' null, and the others at `theta + band`. The band only sets how far apart the
#' simulated truths are, so it also sets how hard the calibration is.
#' The true values do not move with the design's `shift`; they stay around
#' `theta`. With `shift` equal to `band`, the combinations above sit exactly at
#' the decision threshold `theta + shift` (the least favourable case, as in
#' Cheung's calibration for MED), and the others are two bands below it for
#' binary outcomes (one band for normal outcomes).
#'
#' The false discovery rate rises with `eps`, so `eps` is found by bisection
#' on `[lower, upper]`, not by a grid. Every evaluation uses the same seeds
#' for each scenario (common random numbers), which keeps the estimated rate
#' close to monotone in `eps` and makes the evaluations comparable. The search
#' stops when the bracket is narrower than `tol`. If the rate is within
#' `alpha` already at `upper`, `upper` is returned; if it exceeds `alpha` at
#' `lower`, there is no solution and `eps` is `NA` with a warning.
#'
#' Calibration is for the `"cheung"` rule; the Mander & Sweeting rule fixes
#' `eps = 0.5`.
#'
#' @param design A [PipeDesign] with `rule = "cheung"`. Its `eps` is replaced.
#' @param scenarios A list of 0/1 configurations (1 = above `theta`), each a
#'   vector over the combinations or an array over the full grid, or `NULL`
#'   to draw `n_scenarios` with [calibration_scenarios()].
#' @param n_scenarios Number of scenarios when `scenarios` is `NULL`.
#' @param nsim Trials per scenario and evaluation.
#' @param alpha Largest acceptable false discovery rate.
#' @param band Distance of the true values from `theta`, on the scale of the
#'   parameter. The default is 0.1 for binary outcomes, as in [summary()]
#'   for a simulation; it is required for normal outcomes, where it is the
#'   true effect.
#' @param sd True standard deviation for normal outcomes; defaults to the
#'   design's `sigma0`.
#' @param lower,upper Bracket for `eps`, within (0, 1).
#' @param tol Width of the final bracket.
#' @param seed Seed for the scenarios and the trials. The random-number state
#'   is restored afterwards.
#' @return A `PipeCalibration` with `@eps` (the selected value, or `NA`),
#'   `@fdr` and `@tpr` there, `@fdr_se` (Monte Carlo SE, treating trials as
#'   independent), `@design` (the input design with `eps` set), `@trace` (every
#'   `eps` tried, with its FDR, TPR and SE) and `@alpha`.
#' @examples
#' g <- DoseGrid(levels = c(A = 3, B = 3))
#' des <- cheung_design(g, theta = 0.3, eps = 0.5, target = "MTD", outcome = "binary",
#'                      prior_a = 0.3, prior_b = 0.7, N = 30)
#' cal <- calibrate_eps(des, n_scenarios = 3, nsim = 5, tol = 0.1, seed = 1)
#' cal
#' @export
calibrate_eps <- function(design, scenarios = NULL, n_scenarios = 20, nsim = 200, alpha = 0.05,
                          band = NULL, sd = NULL, lower = 0.01, upper = 0.5, tol = 0.01,
                          seed = 1) {
  if (!S7::S7_inherits(design, PipeDesign)) stop("`design` must be a PipeDesign.", call. = FALSE)
  restore_rng()
  if (design@rule != "cheung") {
    stop("Only the \"cheung\" rule has an `eps` to calibrate.", call. = FALSE)
  }
  check_eps(lower)
  check_eps(upper)
  if (!(lower < upper) || !(tol > 0) || !(alpha > 0 && alpha < 1)) {
    stop("Need 0 < lower < upper < 1, `tol` > 0 and 0 < `alpha` < 1.", call. = FALSE)
  }
  grid <- design@grid
  K <- n_combos(grid)
  binary <- design@outcome == "binary"
  if (is.null(band)) {
    if (!binary) stop("Give `band`, the true effect, for normal outcomes.", call. = FALSE)
    band <- 0.1
  }
  if (length(band) != 1L || !(band > 0)) stop("`band` must be a single positive number.", call. = FALSE)
  if (is.null(scenarios)) scenarios <- calibration_scenarios(grid, n_scenarios, seed)
  scenarios <- lapply(scenarios, function(s) {
    s <- if (length(dim(s)) > 1L) grid_from_array(grid, s) else c(s)
    if (length(s) != K || !all(s %in% 0:1) || !is_upset(grid, s)) {
      stop("Each scenario must be a 0/1 up-set over the combinations.", call. = FALSE)
    }
    s
  })
  truths <- lapply(scenarios, function(s) {
    if (binary) {
      design@theta + band * ifelse(s == 1, 1, -1)
    } else {
      list(mean = design@theta + band * s, sd = sd %||% design@sigma0)
    }
  })
  if (binary && any(unlist(truths) <= 0 | unlist(truths) >= 1)) {
    stop("`band` puts a true probability outside (0, 1).", call. = FALSE)
  }
  base <- design_list(design)
  fav_call <- if (design@target == "MED") 1 else 0

  evaluate <- function(eps) {
    d <- base
    d$eps <- eps
    per <- lapply(seq_along(scenarios), function(i) {
      sim <- pipe_simulate(d, truths[[i]], nsim, seed = seed + i, keep_trials = FALSE)
      declared <- sim$gamma == fav_call
      declared[sim$summary$stopped, ] <- FALSE
      fav <- scenarios[[i]] != 1
      if (design@target == "MED") fav <- !fav
      fdr <- rowSums(sim$final_set[, !fav, drop = FALSE]) / pmax(1, rowSums(sim$final_set))
      tpr <- if (any(fav)) rowSums(declared[, fav, drop = FALSE]) / sum(fav) else rep(NA_real_, nsim)
      list(fdr = fdr, tpr = tpr)
    })
    fdr <- unlist(lapply(per, `[[`, "fdr"))
    tpr <- unlist(lapply(per, `[[`, "tpr"))
    data.frame(eps = eps, fdr = mean(fdr), fdr_se = stats::sd(fdr) / sqrt(length(fdr)),
               tpr = mean(tpr, na.rm = TRUE))
  }

  trace <- new.env()
  trace$rows <- list()
  ev <- function(eps) {
    r <- evaluate(eps)
    trace$rows[[length(trace$rows) + 1L]] <- r
    r
  }
  best <- NULL
  hi <- ev(upper)
  if (hi$fdr <= alpha) {
    best <- hi
  } else {
    lo <- ev(lower)
    if (lo$fdr > alpha) {
      warning("The false discovery rate exceeds `alpha` even at `lower`; no `eps` is feasible.",
              call. = FALSE)
    } else {
      best <- lo
      a <- lower
      b <- upper
      while (b - a > tol) {
        mid <- (a + b) / 2
        r <- ev(mid)
        if (r$fdr <= alpha) {
          a <- mid
          best <- r
        } else {
          b <- mid
        }
      }
    }
  }
  trace <- do.call(rbind, trace$rows)
  trace <- trace[order(trace$eps), ]
  eps_hat <- if (is.null(best)) NA_real_ else best$eps
  PipeCalibration(
    design = if (is.na(eps_hat)) {
      design
    } else {
      S7::set_props(design, eps = eps_hat)
    },
    eps = eps_hat, alpha = alpha,
    fdr = if (is.null(best)) NA_real_ else best$fdr,
    fdr_se = if (is.null(best)) NA_real_ else best$fdr_se,
    tpr = if (is.null(best)) NA_real_ else best$tpr,
    trace = trace, n_scenarios = length(scenarios), nsim = nsim, band = band
  )
}

PipeCalibration := new_class(
  package = "gPIPE",
  properties = list(
    design = PipeDesign,
    eps = class_double,
    alpha = class_double,
    fdr = class_double,
    fdr_se = class_double,
    tpr = class_double,
    trace = class_data.frame,
    n_scenarios = class_integer,
    nsim = class_numeric,
    band = class_double
  )
)

method(print, PipeCalibration) <- function(x, digits = 3, ...) {
  cat(sprintf("Epsilon calibration (%s, theta = %g, band = %g): %d scenarios x %d trials per evaluation\n",
              x@design@target, x@design@theta, x@band, x@n_scenarios, as.integer(x@nsim)))
  if (is.na(x@eps)) {
    cat(sprintf("No eps keeps FDR <= %g.\n", x@alpha))
  } else {
    cat(sprintf("eps = %.3f: FDR = %.3f (SE %.3f, target <= %g), TPR = %.3f\n", x@eps, x@fdr,
                x@fdr_se, x@alpha, x@tpr))
  }
  cat("Evaluations:\n")
  print(x@trace, digits = digits, row.names = FALSE)
  invisible(x)
}
