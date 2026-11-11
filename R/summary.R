#' Summarise simulated PIPE trials
#'
#' Turns a [simulate()] result into the tables that papers report: where
#' patients were treated and which combinations were selected, and the
#' operating characteristics of the design. Each metric follows the definition
#' of the paper it comes from, so results can be compared with theirs.
#'
#' Combinations are first classed by the true scenario, using an indifference
#' band of half-width `band` around `theta`. For `target = "MTD"`,
#' *favourable* is a probability of at most `theta - band` (safe) and
#' *unfavourable* is at least `theta + band` (toxic). For `target = "MED"` the
#' sides swap: favourable (effective) is a mean of at least `theta + band`.
#' Combinations inside the band are *indifferent* and count for neither
#' side. With `band = 0`, "above `theta`" means a parameter above it by more
#' than a rounding error.
#'
#' The metrics (`@metrics`, with Monte Carlo standard errors) are:
#' * `fdr` (Cheung et al. 2022, eq. 7): the proportion of the final dose set
#'   (the MED or MTD set) that is unfavourable, 0 when the set is empty.
#'   Indifferent combinations count in the set's size but not as false.
#' * `tpr` (eq. 8): the proportion of favourable combinations that the final
#'   estimate declares favourable (`gamma` = 1 for MED, 0 for MTD). A
#'   stopped trial declares nothing, so it counts 0 in `fdr` and `tpr`; with
#'   many stops FDR is pulled down. FDR is 0 when no combination is
#'   unfavourable in the truth. Both are averaged over trials, and trials
#'   with no favourable combination are left out of `tpr`.
#' * `pcs` (Riviere et al. 2014, Wages): the proportion of trials whose single
#'   pick is a true MTD, a combination with `|truth - theta| <= mtd_tol`. Trials
#'   with no pick count as failures. Binary MTD designs only. It is `NA` when
#'   no combination is that close to `theta`; a wider `mtd_tol` gives a "within
#'   `mtd_tol` of the target" version.
#' * `pct_patients_favourable`, `pct_patients_unfavourable`: the share of
#'   patients treated at each kind of combination (an overdose rate for MTD),
#'   averaged over trials. This is not the pooled share in `@by_combo`, and
#'   the two differ when trials stop early.
#' * `p_set_unfavourable`: the proportion of trials whose final dose set
#'   contains an unfavourable combination.
#' * `mean_recommended` (Mander & Sweeting 2015): the mean number of
#'   combinations in the final set, counting 0 for a trial that stopped.
#'
#' `@bands` holds the Mander & Sweeting recommendation table: the final-set
#' combinations pooled over trials, as a percentage by true-parameter band,
#' with each trial that recommended nothing counted once as "None" in the
#' denominator. The default bands follow Table IV of that paper, by distance
#' from `theta` in percentage points, with the width set by `table_band`
#' (the paper's three bands, with the outer bands split by side). Only the recommendation percentages
#' are given, not their experimentation percentages. Give `breaks`
#' (increasing, left-closed cut points covering the true values, e.g.
#' `c(0, .15, .25, .35, .46, 1)` for their Table III) to use your own.
#'
#' @param object A `PipeSimulation` from [simulate()].
#' @param band Half-width of the indifference band. `NULL` uses
#'   `band_frac * theta` for MTD and `0` for MED.
#' @param band_frac Width of the default MTD band as a fraction of `theta`.
#'   The default 0.25 is the decided calibration band; 0.5 is the sensitivity
#'   analysis. Not used when `band` is given.
#' @param mtd_tol Tolerance for "true MTD" in `pcs` only: a combination is a
#'   true MTD when `|truth - theta| <= mtd_tol`. Separate from the bands; the
#'   default is a rounding tolerance, as in Riviere et al., whose true MTDs sit
#'   exactly at the target.
#' @param table_band Half-width, on the probability scale, of the "within"
#'   bands of the default recommendation table. The default 0.1 is the 10
#'   percentage points of Mander & Sweeting's Table IV; the row labels follow
#'   it. Use `table_band = band` to tie the table to the indifference band.
#' @param breaks Cut points for the recommendation table, or `NULL` for the
#'   bands set by `table_band`.
#' @param ... Not used.
#' @return A `PipeSummary` with `@by_combo` (one row per combination: truth,
#'   percentage of trials selecting it, mean patients, share of patients,
#'   percentage of trials that tried it), `@trials` (mean patients enrolled,
#'   percentage stopped early), `@metrics` and `@bands`.
#' @examples
#' g <- DoseGrid(levels = c(A = 2, B = 2))
#' des <- ms_design(g, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 4, N = 12)
#' sim <- simulate(des, nsim = 20, seed = 1, truth = c(0.1, 0.2, 0.3, 0.5))
#' summary(sim, band = 0)
#' @name summary.PipeSimulation
NULL

PipeSummary := new_class(
  package = "gPIPE",
  properties = list(
    design = PipeDesign,
    band = class_double,
    by_combo = class_data.frame,
    trials = class_data.frame,
    metrics = class_data.frame,
    bands = class_data.frame
  )
)

# Class each combination as favourable, unfavourable or indifferent (see
# `summary.PipeSimulation`). Favourable is a deliberate side, so a value
# exactly at the edge of the band is on that side, within `tol`.
classify_truth <- function(target, theta, truth, band, tol = 1e-9) {
  above <- truth > theta + tol
  if (band == 0) {
    # No indifference zone: a parameter above theta (beyond `tol`) is "above".
    fav <- if (target == "MTD") !above else above
    unfav <- !fav
  } else {
    low <- truth <= theta - band + tol
    high <- truth >= theta + band - tol
    fav <- if (target == "MTD") low else high
    unfav <- if (target == "MTD") high else low
  }
  factor(ifelse(fav, "favourable", ifelse(unfav, "unfavourable", "indifferent")),
         levels = c("favourable", "unfavourable", "indifferent"))
}

mc_row <- function(metric, x) {
  x <- x[!is.na(x)]
  data.frame(metric = metric, estimate = if (length(x)) mean(x) else NA_real_,
             mc_se = if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_)
}

# Mander & Sweeting's recommendation table. `counts` are the final-set
# memberships pooled over trials; `n_none` is the number of trials that
# recommended nothing, counted once as "None" in the denominator (Table IV).
# `width` is the half-width of the "within" bands, and sets their labels.
ms_band_table <- function(counts, truth, theta, breaks, n_none, width = 0.1) {
  if (is.null(breaks)) {
    # Rounded so that a value exactly `width` from theta stays in "within".
    d <- round(truth - theta, 9)
    w <- round(width, 9)
    pts <- format(100 * width)
    lab <- c(sprintf("more than %s below", pts), sprintf("within %s below", pts), "at theta",
             sprintf("within %s above", pts), sprintf("more than %s above", pts))
    idx <- ifelse(d == 0, 3L, ifelse(d < -w, 1L, ifelse(d < 0, 2L, ifelse(d > w, 5L, 4L))))
    band <- factor(lab[idx], levels = lab)
  } else {
    band <- cut(truth, breaks, right = FALSE, include.lowest = TRUE)
  }
  n_combos <- as.vector(table(band))
  total <- as.vector(tapply(counts, band, sum))
  total[is.na(total)] <- 0
  total <- c(total, n_none)
  data.frame(band = c(levels(band), "None"), n_combos = c(n_combos, NA_integer_),
             recommended = total, pct_recommended = 100 * total / sum(total))
}

check_breaks <- function(breaks, truth) {
  if (!is.numeric(breaks) || length(breaks) < 2L || anyNA(breaks) || is.unsorted(breaks, strictly = TRUE)) {
    stop("`breaks` must be at least two increasing numbers.", call. = FALSE)
  }
  if (min(breaks) > min(truth) || max(breaks) < max(truth)) {
    stop("`breaks` must cover the true values, from ", format(min(truth)), " to ", format(max(truth)),
         ".", call. = FALSE)
  }
}

method(summary, PipeSimulation) <- function(object, band = NULL, band_frac = 0.25, mtd_tol = 1e-8,
                                            table_band = 0.1, breaks = NULL, ...) {
  des <- object@design
  K <- n_combos(des@grid)
  nsim <- nrow(object@summary)
  truth <- object@truth
  if (length(truth) != K) {
    stop("The simulation has no `truth` for every combination.", call. = FALSE)
  }
  if (is.null(band)) band <- if (des@target == "MTD") band_frac * des@theta else 0
  if (!is.numeric(band) || length(band) != 1L || is.na(band) || band < 0) {
    stop("`band` must be a single number, 0 or more.", call. = FALSE)
  }
  if (!is.numeric(table_band) || length(table_band) != 1L || is.na(table_band) || table_band <= 0) {
    stop("`table_band` must be a single positive number.", call. = FALSE)
  }
  if (!is.null(breaks)) check_breaks(breaks, truth)
  cls <- classify_truth(des@target, des@theta, truth, band)
  fav <- cls == "favourable"
  unfav <- cls == "unfavourable"

  stopped <- object@summary$stopped
  fin <- object@final_set
  # FDR is over the final dose set (Cheung eq. 7); TPR over every combination
  # the final estimate declares favourable (eq. 8). A stopped trial declares nothing.
  # With no unfavourable combination in the truth nothing can be a false
  # discovery, so FDR is 0 (set explicitly; an empty set also gives 0).
  fdr <- if (any(unfav)) rowSums(fin[, unfav, drop = FALSE]) / pmax(1, rowSums(fin)) else numeric(nsim)
  declared <- if (des@target == "MED") object@gamma == 1 else object@gamma == 0
  declared[stopped, ] <- FALSE
  tpr <- if (any(fav)) rowSums(declared[, fav, drop = FALSE]) / sum(fav) else rep(NA_real_, nsim)

  pcs <- NULL
  true_mtd <- abs(truth - des@theta) <= mtd_tol
  if (des@outcome == "binary" && des@target == "MTD") {
    pick <- object@summary$pick
    pcs <- if (any(true_mtd)) as.numeric(!is.na(pick) & true_mtd[pick]) else rep(NA_real_, nsim)
  }
  n_tot <- rowSums(object@n)
  share <- function(sel) if (any(sel)) 100 * rowSums(object@n[, sel, drop = FALSE]) / n_tot else rep(NA_real_, nsim)
  set_unfav <- if (any(unfav)) {
    as.numeric(rowSums(fin[, unfav, drop = FALSE]) > 0)
  } else {
    rep(NA_real_, nsim)
  }
  metrics <- rbind(
    mc_row("fdr", fdr), mc_row("tpr", tpr),
    if (!is.null(pcs)) mc_row("pcs", pcs),
    mc_row("pct_patients_favourable", share(fav)),
    mc_row("pct_patients_unfavourable", share(unfav)),
    mc_row("p_set_unfavourable", set_unfav),
    mc_row("mean_recommended", rowSums(fin))
  )

  by_combo <- data.frame(
    combo = combo_labels(des@grid), truth = truth, class = cls,
    pct_selected = 100 * colMeans(fin), mean_n = colMeans(object@n),
    pct_patients = 100 * colSums(object@n) / sum(object@n),
    pct_trials_tried = 100 * colMeans(object@n > 0)
  )
  PipeSummary(
    design = des, band = band, by_combo = by_combo,
    trials = data.frame(nsim = nsim, mean_n = mean(object@summary$n_enrolled),
                        pct_stopped = 100 * mean(stopped)),
    metrics = metrics,
    bands = ms_band_table(colSums(fin), truth, des@theta, breaks, sum(rowSums(fin) == 0), table_band)
  )
}

method(print, PipeSummary) <- function(x, digits = 3, ...) {
  d <- x@design
  cat(sprintf("Summary of %d simulated PIPE trials (%s, theta = %g, band = %g)\n",
              x@trials$nsim, d@target, d@theta, x@band))
  cat(sprintf("Mean patients: %.1f; stopped early: %.1f%%\n\n", x@trials$mean_n,
              x@trials$pct_stopped))
  cat("By combination:\n")
  print(x@by_combo, digits = digits, row.names = FALSE)
  cat("\nOperating characteristics (with Monte Carlo SE):\n")
  print(x@metrics, digits = digits, row.names = FALSE)
  cat("\nRecommended combinations by true-parameter band (pooled over trials):\n")
  print(x@bands, digits = digits, row.names = FALSE)
  invisible(x)
}
