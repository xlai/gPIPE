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
# True toxicity of Riviere scenarios 2 and 3 (same combination order), used
# as increasing priors that are too high (scenario 2) or too low (scenario 3)
# for scenario 1. Taken from riviere2014_scenarios_all14.csv.
mtd_truth_s2 <- c(0.15, 0.30, 0.45, 0.50, 0.60,
                  0.30, 0.45, 0.50, 0.60, 0.75,
                  0.45, 0.55, 0.60, 0.70, 0.80)
mtd_truth_s3 <- c(0.02, 0.07, 0.10, 0.15, 0.30,
                  0.07, 0.10, 0.15, 0.30, 0.45,
                  0.10, 0.15, 0.30, 0.45, 0.55)
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
#   prior_median, prior_strength, prior_method: Beta prior with this median
#                 and strength a + b;
#   prior_form:   "flat" (the same median at every combination), or an
#                 increasing prior whose medians are the truth of
#                 "scenario1" ("correct"), "scenario3" (too low) or
#                 "scenario2" (too high); then prior_median is not used.
#                 Beta(1,1) is flat with median 0.5, strength 2; Beta(0.5,0.5)
#                 is flat with median 0.5, strength 1 (exact method);
#   calibration:  file in results/ that gives epsilon when eps is NA;
#   shift:        shift of the decision threshold (theta + shift), cheung
#                 only; 0 is the usual rule;
#   startup:      "none" (the usual rule from the first cohort) or "diagonal"
#                 (diagonal start-up phase until the first toxicity), cheung
#                 only;
#   N, cohort, start (NA = the default, combination 1), delta, within_set;
#   admis, constraint, safety: M&S options ("ms" only; safety NA = none);
#   seed:         seed of this configuration's simulation (unique per row).
mtd_configs <- function() {
  K <- nrow(mtd_truth) * ncol(mtd_truth)
  row <- function(id, factor, rule = "cheung", eps = NA_real_, taper = TRUE, weights = "n+1",
                  selection = "equal", prior_median = mtd_theta, prior_strength = 1 / K,
                  prior_method = "exact", prior_form = "flat",
                  calibration = "calibration_baseline.csv", N = 60, cohort = 3, start = NA_real_, delta = 0.1,
                  within_set = TRUE, admis = NA_character_, constraint = NA_character_,
                  safety = NA_real_, shift = 0, startup = "none", seed) {
    data.frame(id = id, factor = factor, rule = rule, eps = eps, taper = taper, weights = weights,
               selection = selection, prior_median = prior_median, prior_strength = prior_strength,
               prior_method = prior_method, prior_form = prior_form,
               calibration = calibration, N = N, cohort = cohort, start = start, delta = delta,
               within_set = within_set, admis = admis, constraint = constraint, safety = safety,
               shift = shift, startup = startup, seed = seed,
               stringsAsFactors = FALSE)
  }
  rbind(
    row("baseline", "baseline", seed = 4501),
    row("eps_0.5", "eps", eps = 0.5, seed = 4502),
    # The M&S design as ms_design() has it by default, except the prior,
    # which is the same flat prior as the baseline (bisection is its default
    # way to find the median).
    row("ref_ms", "reference", rule = "ms", eps = NA_real_, weights = NA_character_,
        taper = NA, selection = "min_S", prior_method = "bisection", admis = "closest",
        constraint = "neighbouring", safety = 0.8, seed = 4503),
    # One factor at a time from the baseline. Epsilon stays at the baseline
    # calibration (NA) unless the row is about epsilon or the calibration.
    row("eps_0.4", "eps", eps = 0.4, seed = 4504),
    row("eps_0.3", "eps", eps = 0.3, seed = 4505),
    row("eps_0.15", "eps", eps = 0.15, seed = 4506),
    row("eps_0.05", "eps", eps = 0.05, seed = 4507),
    row("alpha_0.10", "calibration alpha", calibration = "calibration_alpha_0.10.csv", seed = 4508),
    row("alpha_0.20", "calibration alpha", calibration = "calibration_alpha_0.20.csv", seed = 4509),
    row("taper_off", "taper", taper = FALSE, seed = 4510),
    row("weights_1", "weights", weights = "1", seed = 4511),
    row("prior_strength_1", "prior strength", prior_strength = 1, seed = 4512),
    row("prior_strength_4", "prior strength", prior_strength = 4, seed = 4513),
    row("prior_beta_1_1", "prior form", prior_median = 0.5, prior_strength = 2, seed = 4514),
    row("prior_beta_0.5_0.5", "prior form", prior_median = 0.5, prior_strength = 1, seed = 4515),
    row("prior_correct", "prior form", prior_form = "scenario1", seed = 4516),
    row("prior_too_low", "prior form", prior_form = "scenario3", seed = 4517),
    row("prior_too_high", "prior form", prior_form = "scenario2", seed = 4518),
    row("prior_correct_strength_1", "prior form x strength", prior_form = "scenario1",
        prior_strength = 1, seed = 4519),
    row("selection_prob", "selection", selection = "prob", seed = 4520),
    row("cohort_1", "cohort", cohort = 1, seed = 4521),
    row("N_30", "N", N = 30, seed = 4522),
    row("N_120", "N", N = 120, seed = 4523),
    row("start_top", "start", start = 15, seed = 4524),
    row("single_pick_off", "single pick", within_set = FALSE, seed = 4525),
    # Selection of the next dose from the sampling set.
    row("selection_fewest", "selection", selection = "fewest", seed = 4526),
    row("selection_inverse_n", "selection", selection = "inverse_n", seed = 4527),
    row("selection_interval", "selection", selection = "interval", seed = 4528),
    # Band shift: the shift is the band b of its calibration (truths at
    # theta +- b, so the toxic ones sit on the threshold theta + b), and
    # epsilon is calibrated again for that design. summary() keeps band 0.1.
    row("shift_0.1", "band shift", shift = 0.1, calibration = "calibration_shift_0.1.csv", seed = 4529),
    row("shift_0.075", "band shift", shift = 0.075, calibration = "calibration_shift_0.075.csv",
        seed = 4530),
    row("shift_0.05", "band shift", shift = 0.05, calibration = "calibration_shift_0.05.csv", seed = 4531),
    row("shift_0.025", "band shift", shift = 0.025, calibration = "calibration_shift_0.025.csv",
        seed = 4532),
    # The shift of shift_0.1 with the baseline epsilon (not recalibrated).
    row("shift_0.1_eps_fixed", "band shift, fixed eps", shift = 0.1, seed = 4533),
    # Exploratory rows: each changes more than one thing from the baseline.
    row("eps_0.4_pick_all", "combination", eps = 0.4, within_set = FALSE, seed = 4534),
    row("eps_0.3_pick_all", "combination", eps = 0.3, within_set = FALSE, seed = 4535),
    row("shift_0.1_pick_all", "combination", shift = 0.1, calibration = "calibration_shift_0.1.csv",
        within_set = FALSE, seed = 4536),
    row("shift_0.1_pick_all_fewest", "combination", shift = 0.1, selection = "fewest",
        calibration = "calibration_shift_0.1.csv", within_set = FALSE, seed = 4537),
    row("shift_0.1_pick_all_interval", "combination", shift = 0.1, selection = "interval",
        calibration = "calibration_shift_0.1.csv", within_set = FALSE, seed = 4538),
    row("shift_0.05_pick_all", "combination", shift = 0.05, calibration = "calibration_shift_0.05.csv",
        within_set = FALSE, seed = 4539),
    # Start-up phase (#51). Each "_startup" row is the row of the same name
    # without its suffix, with startup = "diagonal" and its own seed. The
    # one row below has no start-up phase: it is the partner of
    # shift_0.1_eps_fixed_pick_all_startup.
    row("shift_0.1_eps_fixed_pick_all", "startup", shift = 0.1, within_set = FALSE, seed = 4540),
    row("baseline_startup", "startup", startup = "diagonal", seed = 4541),
    row("single_pick_off_startup", "startup", startup = "diagonal", within_set = FALSE, seed = 4542),
    row("eps_0.3_startup", "startup", eps = 0.3, startup = "diagonal", seed = 4543),
    row("eps_0.4_startup", "startup", eps = 0.4, startup = "diagonal", seed = 4544),
    row("eps_0.3_pick_all_startup", "startup", eps = 0.3, startup = "diagonal", within_set = FALSE,
        seed = 4545),
    row("eps_0.4_pick_all_startup", "startup", eps = 0.4, startup = "diagonal", within_set = FALSE,
        seed = 4546),
    row("selection_fewest_startup", "startup", selection = "fewest", startup = "diagonal", seed = 4547),
    row("selection_inverse_n_startup", "startup", selection = "inverse_n", startup = "diagonal",
        seed = 4548),
    row("shift_0.1_pick_all_fewest_startup", "startup", shift = 0.1, selection = "fewest",
        calibration = "calibration_shift_0.1.csv", startup = "diagonal", within_set = FALSE, seed = 4549),
    row("shift_0.1_eps_fixed_startup", "startup", shift = 0.1, startup = "diagonal", seed = 4550),
    row("shift_0.1_eps_fixed_pick_all_startup", "startup", shift = 0.1, startup = "diagonal",
        within_set = FALSE, seed = 4551)
  )
}

# Epsilon of a configuration: its own value, or for NA the calibrated value
# stored in the file of its `calibration` column. NA for the M&S rule.
mtd_eps <- function(cfg, dir = "results") {
  if (cfg$rule != "cheung") return(NA_real_)
  if (!is.na(cfg$eps)) return(cfg$eps)
  f <- file.path(dir, cfg$calibration)
  if (!file.exists(f)) stop("Run calibrate.R first (see README): ", f, " is missing.", call. = FALSE)
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
  med <- switch(cfg$prior_form,
                flat = rep(cfg$prior_median, K),
                scenario1 = c(mtd_truth), scenario2 = mtd_truth_s2, scenario3 = mtd_truth_s3,
                stop("Unknown prior_form: ", cfg$prior_form, call. = FALSE))
  pr <- gPIPE:::beta_prior_from_median(med, cfg$prior_strength, method = cfg$prior_method)
  cheung_design(grid, theta = mtd_theta, eps = eps, target = "MTD", outcome = "binary",
                taper = cfg$taper, weights = cfg$weights, selection = cfg$selection,
                shift = cfg$shift, startup = cfg$startup, prior_a = pr$a, prior_b = pr$b, N = cfg$N,
                cohort = cfg$cohort, start = start,
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

# Start-up phase (#51), one value per trial from a run with startup =
# "diagonal". The phase follows the diagonal from combination 1 until every
# drug is at its top level, so on a full grid it has at most max(levels)
# cohorts (5 here). `n_startup` is the number of cohorts treated in it. A
# trial whose phase is shorter than that was ended by a toxicity.
mtd_startup_length <- function(grid = mtd_grid()) max(grid@levels)

# Summary of the start-up phase of one simulation, one row: mean number of
# cohorts in the phase (with its Monte Carlo SE) and the percentage of trials
# whose phase ended by a toxicity before the end of the diagonal.
mtd_startup_summarise <- function(id, sim) {
  n_startup <- sim@summary$n_startup
  ended <- 100 * mean(n_startup < mtd_startup_length(sim@design@grid))
  data.frame(id = id, mean_n_startup = mean(n_startup),
             se_n_startup = stats::sd(n_startup) / sqrt(length(n_startup)),
             pct_phase_ended_by_toxicity = ended, stringsAsFactors = FALSE)
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
