# Component sensitivity study of PIPE with epsilon-tapering (MTD, binary
# toxicity): the scenario, the configuration table, the design built from a
# row, and the summary of one simulation. Used by calibrate.R, run.R and
# report.R. Run those from this folder.
#
# To add a configuration, add a row in mtd_configs(). Nothing else changes.

mtd_theta <- 0.30
mtd_band <- 0.1
mtd_nsim <- 10000L

# Riviere et al. (2014), scenario 1: 5 levels of drug A x 3 levels of drug B.
# Rows = drug A level, columns = drug B level. Written column by column, so
# this is also gPIPE's combination order (A varies fastest: A1.B1, A2.B1, ...).
mtd_truth <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
                      0.10, 0.15, 0.30, 0.45, 0.55,
                      0.15, 0.30, 0.45, 0.50, 0.60), nrow = 5, ncol = 3)
mtd_grid <- function() DoseGrid(levels = c(A = nrow(mtd_truth), B = ncol(mtd_truth)))

# Every monotone true configuration of the grid, as 0/1 vectors over the
# combinations (1 = above theta): for each level of drug B, the top c levels
# of drug A are above, with c not falling as B rises. 56 on the 5 x 3 grid,
# all-0 and all-1 included. Cheung et al. (2022) calibrate over all such
# configurations of their grid in the same way.
mtd_upsets <- function(nA = nrow(mtd_truth), nB = ncol(mtd_truth)) {
  cuts <- as.matrix(expand.grid(rep(list(0:nA), nB)))
  cuts <- cuts[apply(cuts, 1, function(x) all(diff(x) >= 0)), , drop = FALSE]
  lapply(seq_len(nrow(cuts)), function(i) {
    as.integer(c(vapply(cuts[i, ], function(k) rev(seq_len(nA) <= k), logical(nA))))
  })
}

# True MTDs: the combinations with truth = theta (Riviere's definition).
mtd_is_mtd <- function(truth = c(mtd_truth), theta = mtd_theta) abs(truth - theta) < 1e-8
# The true contour: toxicity at or above theta.
mtd_is_contour <- function(truth = c(mtd_truth), theta = mtd_theta) truth >= theta - 1e-8

# One row per configuration.
#   id, factor:   name of the row, and the factor it varies ("baseline",
#                 "eps", "reference", ...);
#   rule:         "cheung" or "ms";
#   eps:          epsilon (cheung only). NA means the value calibrated by
#                 calibrate.R, read from results/calibration_baseline.csv;
#   taper, weights, selection: as in PipeDesign (for "ms", selection is the
#                 M&S choice and taper and weights are not used);
#   prior_median, prior_strength, prior_method: flat Beta prior on every
#                 combination, with this median and strength a + b;
#   N, cohort, start (NA = the default, combination 1), delta, within_set;
#   admis, constraint, safety: M&S options ("ms" only; safety NA = none);
#   seed:         seed of this configuration's simulation (unique per row).
mtd_configs <- function() {
  K <- nrow(mtd_truth) * ncol(mtd_truth)
  row <- function(id, factor, rule = "cheung", eps = NA_real_, taper = TRUE, weights = "n+1",
                  selection = "equal", prior_median = mtd_theta, prior_strength = 1 / K,
                  prior_method = "exact", N = 60, cohort = 3, start = NA_real_, delta = 0.1,
                  within_set = TRUE, admis = NA_character_, constraint = NA_character_,
                  safety = NA_real_, seed) {
    data.frame(id = id, factor = factor, rule = rule, eps = eps, taper = taper, weights = weights,
               selection = selection, prior_median = prior_median, prior_strength = prior_strength,
               prior_method = prior_method, N = N, cohort = cohort, start = start, delta = delta,
               within_set = within_set, admis = admis, constraint = constraint, safety = safety,
               seed = seed, stringsAsFactors = FALSE)
  }
  rbind(
    row("baseline", "baseline", seed = 4501),
    row("eps_0.5", "eps", eps = 0.5, seed = 4502),
    # The M&S design as ms_design() has it by default, except the prior,
    # which is the same flat prior as the baseline (bisection is its default
    # way to find the median).
    row("ref_ms", "reference", rule = "ms", eps = NA_real_, weights = NA_character_,
        taper = NA, selection = "min_S", prior_method = "bisection", admis = "closest",
        constraint = "neighbouring", safety = 0.8, seed = 4503)
  )
}

# Epsilon of a configuration: its own value, or for NA the calibrated value
# stored in results/calibration_baseline.csv. NA for the M&S rule.
mtd_eps <- function(cfg, dir = "results") {
  if (cfg$rule != "cheung") return(NA_real_)
  if (!is.na(cfg$eps)) return(cfg$eps)
  f <- file.path(dir, "calibration_baseline.csv")
  if (!file.exists(f)) stop("Run calibrate.R first: ", f, " is missing.", call. = FALSE)
  cal <- utils::read.csv(f)
  eps <- unique(cal$eps[cal$selected])
  if (length(eps) != 1L) stop("The calibration file does not name one selected eps.", call. = FALSE)
  eps
}

# Turn a row of mtd_configs() into a PipeDesign. `eps` is the value from
# mtd_eps().
mtd_design <- function(cfg, eps = mtd_eps(cfg)) {
  grid <- mtd_grid()
  K <- n_combos(grid)
  start <- if (is.na(cfg$start)) NULL else cfg$start
  if (cfg$rule == "ms") {
    return(ms_design(grid, theta = mtd_theta, prior_median = cfg$prior_median,
                     prior_strength = cfg$prior_strength, prior_method = cfg$prior_method,
                     N = cfg$N, cohort = cfg$cohort, start = start, admis = cfg$admis,
                     constraint = cfg$constraint, safety = if (is.na(cfg$safety)) NULL else cfg$safety,
                     selection = cfg$selection, delta = cfg$delta, within_set = cfg$within_set))
  }
  pr <- gPIPE:::beta_prior_from_median(rep(cfg$prior_median, K), cfg$prior_strength,
                                       method = cfg$prior_method)
  cheung_design(grid, theta = mtd_theta, eps = eps, target = "MTD", outcome = "binary",
                taper = cfg$taper, weights = cfg$weights, selection = cfg$selection,
                prior_a = pr$a, prior_b = pr$b, N = cfg$N, cohort = cfg$cohort, start = start,
                delta = cfg$delta, within_set = cfg$within_set)
}

# ---- Mechanism metrics, one value per trial ---------------------------------
# `doses` is the combination given to each cohort, in order (a trial's
# `doses` from simulate(..., keep_trials = TRUE)); `levels` is the grid's
# matrix of levels, one row per combination (grid@combos).

# Number of the first cohort given a combination on the true contour (truth
# >= theta); the first cohort is number 1. NA when no cohort reaches it.
first_contour_cohort <- function(doses, contour) {
  hit <- which(contour[doses])
  if (length(hit)) hit[1] else NA_integer_
}

# Number of different combinations given.
n_distinct_doses <- function(doses) length(unique(doses))

# TRUE when some cohort after the first is a skip: its combination is more
# than one level above every earlier combination in some drug. Not a skip: it
# is within one level above, in every drug, of at least one earlier cohort's
# combination (a lower level in some drug is within one level). With
# `two_sided` the test is |difference| <= 1 in every drug instead.
has_skip <- function(doses, levels, two_sided = FALSE) {
  for (t in seq_along(doses)[-1]) {
    d <- levels[doses[t], ]
    before <- levels[doses[seq_len(t - 1)], , drop = FALSE]
    diff <- sweep(before, 2, d, function(b, x) x - b)  # new level - earlier level, per drug
    if (two_sided) diff <- abs(diff)
    near <- rowSums(diff > 1) == 0
    if (!any(near)) return(TRUE)
  }
  FALSE
}

# Mean and Monte Carlo SE of a vector of per-trial values (NA left out).
mean_se <- function(x) {
  x <- x[!is.na(x)]
  c(mean(x), if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_)
}

# Per-trial mechanism values from a simulation run with keep_trials = TRUE.
mtd_trial_metrics <- function(sim, truth = sim@truth, theta = mtd_theta, two_sided = FALSE) {
  levels <- sim@design@grid@combos
  contour <- mtd_is_contour(truth, theta)
  mtd <- mtd_is_mtd(truth, theta)
  doses <- lapply(sim@trials, `[[`, "doses")
  data.frame(
    first_contour_cohort = vapply(doses, first_contour_cohort, numeric(1), contour = contour),
    n_distinct = vapply(doses, n_distinct_doses, numeric(1)),
    skip = vapply(doses, has_skip, numeric(1), levels = levels, two_sided = two_sided),
    n_dlt = rowSums(sim@s1),
    set_has_mtd = as.numeric(rowSums(sim@final_set[, mtd, drop = FALSE]) > 0),
    n_mtd_tried = rowSums(sim@n[, mtd, drop = FALSE] > 0)
  )
}

# Units: fdr, tpr, pcs, p_set_false are proportions (0 to 1); names starting
# with pct_ are percentages (0 to 100); the others are counts or means.
mtd_metric_names <- function() {
  list(
    overall = c("fdr", "tpr", "pcs", "pct_patients_above", "pct_patients_in_band",
                "pct_patients_below", "p_set_false", "mean_recommended", "pct_stopped",
                "mean_patients"),
    by_combo = c("pct_selected", "pct_pick", "mean_n"),
    mechanism = c("first_contour_cohort", "pct_never_contour", "n_distinct", "pct_skip",
                  "n_dlt", "pct_set_has_mtd", "n_mtd_tried")
  )
}

# Summary of one simulation (run with keep_trials = TRUE) in long format:
# id, metric, combo (a combination label, or "" for a whole-trial metric),
# estimate and its Monte Carlo SE.
mtd_summarise <- function(id, sim, band = mtd_band) {
  s <- summary(sim, band = band)
  labels <- s@by_combo$combo
  nsim <- nrow(sim@summary)
  truth <- sim@truth
  long <- function(metric, est, se, combo = "") {
    data.frame(id = id, metric = metric, combo = combo, estimate = est, mc_se = se,
               stringsAsFactors = FALSE)
  }
  m <- s@metrics
  out <- list(long(m$metric, m$estimate, m$mc_se))

  stopped <- sim@summary$stopped
  enrolled <- sim@summary$n_enrolled
  out$trials <- rbind(
    long("pct_stopped", 100 * mean(stopped), 100 * sqrt(mean(stopped) * (1 - mean(stopped)) / nsim)),
    long("mean_patients", mean(enrolled), stats::sd(enrolled) / sqrt(nsim))
  )

  # Per combination. pct_selected: final-set membership. pct_pick: the single
  # pick (Riviere), so these add up to pcs over the true MTDs.
  p_sel <- colMeans(sim@final_set)
  pick <- sim@summary$pick
  p_pick <- tabulate(pick[!is.na(pick)], length(labels)) / nsim
  out$combo <- rbind(
    long("pct_selected", 100 * p_sel, 100 * sqrt(p_sel * (1 - p_sel) / nsim), labels),
    long("pct_pick", 100 * p_pick, 100 * sqrt(p_pick * (1 - p_pick) / nsim), labels),
    long("mean_n", colMeans(sim@n), apply(sim@n, 2, stats::sd) / sqrt(nsim), labels)
  )

  tm <- mtd_trial_metrics(sim, truth)
  never <- is.na(tm$first_contour_cohort)
  est <- function(x) mean_se(x)
  rows <- list(
    first_contour_cohort = est(tm$first_contour_cohort),
    pct_never_contour = 100 * est(as.numeric(never)),
    n_distinct = est(tm$n_distinct),
    pct_skip = 100 * est(tm$skip),
    n_dlt = est(tm$n_dlt),
    pct_set_has_mtd = 100 * est(tm$set_has_mtd),
    n_mtd_tried = est(tm$n_mtd_tried)
  )
  out$mech <- long(names(rows), vapply(rows, `[`, 0, 1), vapply(rows, `[`, 0, 2))
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
