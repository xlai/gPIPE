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
#' side. With `band = 0`, "above `theta`" means a parameter strictly above it.
#'
#' The metrics (`@metrics`, with Monte Carlo standard errors) are:
#' * `fdr`, `tpr` (Cheung et al. 2022): among the combinations the final
#'   estimate declares favourable (`gamma` = 1 for MED, 0 for MTD), the
#'   proportion that are unfavourable (0 when none is declared), and the
#'   proportion of favourable combinations that are declared. Averaged over
#'   trials; trials with no favourable combination are left out of `tpr`.
#' * `pcs` (Riviere et al. 2014, Wages): the proportion of trials whose single
#'   pick is a true MTD, a combination with `|truth - theta| <= tol`. Trials
#'   with no pick count as failures. Binary MTD designs only.
#' * `pct_patients_favourable`, `pct_patients_unfavourable`: the share of
#'   patients treated at each kind of combination (an overdose rate for MTD).
#' * `p_set_unfavourable`: the proportion of trials whose final dose set
#'   contains an unfavourable combination.
#' * `mean_recommended` (Mander & Sweeting 2015): the mean number of
#'   combinations in the final set, counting 0 for a trial that stopped.
#'
#' `@bands` holds the Mander & Sweeting recommendation table: the final-set
#' combinations pooled over trials, as a percentage by true-parameter band.
#' The default bands are Table IV of that paper, by distance from `theta` in
#' percentage points; give `breaks` (left-closed cut points on the true
#' parameter, e.g. `c(0, .15, .25, .35, .46, 1)` for their Table III) to use
#' your own.
#'
#' @param object A `PipeSimulation` from [simulate()].
#' @param band Half-width of the indifference band. `NULL` uses `0.25 * theta`
#'   for MTD (the calibration band) and `0` for MED.
#' @param tol Tolerance for "true MTD" in `pcs`.
#' @param breaks Cut points for the recommendation table, or `NULL` for the
#'   Mander & Sweeting Table IV bands.
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
  above <- truth > theta
  if (band == 0) {
    # No indifference zone: a parameter strictly above theta is "above".
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

ms_band_table <- function(counts, truth, theta, breaks) {
  if (is.null(breaks)) {
    d <- 100 * (truth - theta)
    d[abs(d) < 1e-7] <- 0
    lab <- c("more than 10 below", "1-10 below", "at theta", "1-10 above", "more than 10 above")
    band <- factor(lab[ifelse(d == 0, 3L, ifelse(d < -10, 1L, ifelse(d < 0, 2L,
                                                                      ifelse(d > 10, 5L, 4L))))],
                   levels = lab)
  } else {
    band <- cut(truth, breaks, right = FALSE, include.lowest = TRUE)
  }
  n_combos <- as.vector(table(band))
  total <- tapply(counts, band, sum)
  total[is.na(total)] <- 0
  data.frame(band = levels(band), n_combos = n_combos, recommended = as.vector(total),
             pct_recommended = if (sum(total) > 0) 100 * as.vector(total) / sum(total) else NA_real_)
}

method(summary, PipeSimulation) <- function(object, band = NULL, tol = 1e-8, breaks = NULL, ...) {
  des <- object@design
  K <- n_combos(des@grid)
  nsim <- nrow(object@summary)
  truth <- object@truth
  if (length(truth) != K) {
    stop("The simulation has no `truth` for every combination.", call. = FALSE)
  }
  if (is.null(band)) band <- if (des@target == "MTD") 0.25 * des@theta else 0
  if (length(band) != 1L || band < 0) stop("`band` must be a single number, 0 or more.", call. = FALSE)
  cls <- classify_truth(des@target, des@theta, truth, band)
  fav <- cls == "favourable"
  unfav <- cls == "unfavourable"

  stopped <- object@summary$stopped
  fin <- object@final_set
  # Combinations declared favourable by the final estimate; none after a stop.
  declared <- if (des@target == "MED") object@gamma == 1 else object@gamma == 0
  declared[stopped, ] <- FALSE
  n_decl <- rowSums(declared)
  fdr <- rowSums(declared[, unfav, drop = FALSE]) / pmax(1, n_decl)
  tpr <- if (any(fav)) rowSums(declared[, fav, drop = FALSE]) / sum(fav) else rep(NA_real_, nsim)

  pcs <- NULL
  true_mtd <- abs(truth - des@theta) <= tol
  if (des@outcome == "binary" && des@target == "MTD") {
    pick <- object@summary$pick
    pcs <- if (any(true_mtd)) as.numeric(!is.na(pick) & true_mtd[pick]) else rep(NA_real_, nsim)
  }
  n_tot <- rowSums(object@n)
  share <- function(sel) if (any(sel)) 100 * rowSums(object@n[, sel, drop = FALSE]) / n_tot else rep(NA_real_, nsim)
  metrics <- rbind(
    mc_row("fdr", fdr), mc_row("tpr", tpr),
    if (!is.null(pcs)) mc_row("pcs", pcs),
    mc_row("pct_patients_favourable", share(fav)),
    mc_row("pct_patients_unfavourable", share(unfav)),
    mc_row("p_set_unfavourable", if (any(unfav)) as.numeric(rowSums(fin[, unfav, drop = FALSE]) > 0) else rep(NA_real_, nsim)),
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
    bands = ms_band_table(colSums(fin), truth, des@theta, breaks)
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
