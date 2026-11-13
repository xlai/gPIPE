# Print the stored results: one row per configuration with the main metrics,
# grouped by factor with the baseline repeated at the top of each group, the
# difference from the baseline in PCS, FDR and TPR (with its SE and z), the
# start-up phase pairs (without and with the phase), the selection grid of one
# configuration, a table of all configurations sorted by PCS, and Riviere's
# reference line.
# Usage: Rscript report.R [id = baseline]
# Run from this folder.
library(gPIPE)
source("common.R")
options(width = 250)
args <- commandArgs(trailingOnly = TRUE)
show_id <- if (length(args) >= 1) args[1] else "baseline"
res <- utils::read.csv("results/mtd_sensitivity.csv",
                       colClasses = c(combo = "character"))
runs <- utils::read.csv("results/mtd_sensitivity_runs.csv")
cfgs <- mtd_configs()
ids <- cfgs$id[cfgs$id %in% res$id]
if (!show_id %in% ids) stop("No results for configuration ", show_id)

# Whole-trial metrics. Percent metrics as percent, fdr, tpr and pcs as
# percent too (x 100), so they sit on the same scale as Riviere's PCS.
get <- function(id, metric, combo = "") res[res$id == id & res$metric == metric & res$combo == combo, ]
cell <- function(id, metric, scale = 1, digits = 1) {
  r <- get(id, metric)
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), scale * r$estimate, scale * r$mc_se)
}
cols <- list(
  PCS = c("pcs", 100), FDR = c("fdr", 100), TPR = c("tpr", 100),
  `%above` = c("pct_patients_above", 1), `%band` = c("pct_patients_in_band", 1),
  `%below` = c("pct_patients_below", 1), P.setFalse = c("p_set_false", 100),
  nRec = c("mean_recommended", 1, 2), `%stop` = c("pct_stopped", 1), nPat = c("mean_patients", 1),
  cohToContour = c("first_contour_cohort", 1), `%noContour` = c("pct_never_contour", 1),
  nDistinct = c("n_distinct", 1), `%skip` = c("pct_skip", 1), DLTs = c("n_dlt", 1),
  `%setHasMTD` = c("pct_set_has_mtd", 1), nMTDtried = c("n_mtd_tried", 1, 2)
)
tab <- data.frame(id = ids, factor = cfgs$factor[match(ids, cfgs$id)],
                  shift = ifelse(cfgs$shift[match(ids, cfgs$id)] == 0, "-",
                                 sprintf("%g", cfgs$shift[match(ids, cfgs$id)])),
                  eps = ifelse(is.na(runs$eps[match(ids, runs$id)]), "-",
                               sprintf("%.3f", runs$eps[match(ids, runs$id)])),
                  stringsAsFactors = FALSE)
for (nm in names(cols)) {
  a <- cols[[nm]]
  tab[[nm]] <- vapply(ids, function(i) cell(i, a[1], as.numeric(a[2]), if (length(a) > 2) as.integer(a[3]) else 1L),
                      "")
}
cat(sprintf("MTD sensitivity study: %d trials per configuration, theta = %g, band = %g.\n",
            runs$nsim[1], mtd_theta, mtd_band))
cat("Estimate (Monte Carlo SE). PCS, FDR, TPR and P.setFalse are in percent.\n")
cat("Rows other than the eps, calibration and band shift rows use the baseline eps.\n")
cat("shift: shift of the decision threshold (theta + shift); the summary band is 0.1 in every row.\n")
cat("The combination rows are exploratory: each changes more than one thing from the baseline.\n")

# One block per factor, in the order of mtd_configs(); the baseline row is
# repeated at the top of each block (except its own).
groups <- unique(tab$factor[tab$id != "baseline"])
show_groups <- function(t) {
  for (g in groups) {
    cat(sprintf("\n== %s ==\n", g))
    b <- t[t$id == "baseline", , drop = FALSE]
    print(rbind(b, t[t$factor == g, , drop = FALSE]), row.names = FALSE, right = FALSE)
  }
}
show_groups(tab)

# Difference from the baseline in PCS, FDR and TPR, in percentage points. The
# runs are independent, so SE = sqrt(se1^2 + se2^2) and z = difference / SE.
dif <- function(id, metric, ref = "baseline") {
  a <- get(id, metric)
  b <- get(ref, metric)
  d <- 100 * (a$estimate - b$estimate)
  se <- 100 * sqrt(a$mc_se^2 + b$mc_se^2)
  c(d = d, se = se, z = if (se > 0) d / se else NA_real_)
}
dtab <- data.frame(id = ids, factor = tab$factor, shift = tab$shift, eps = tab$eps,
                   stringsAsFactors = FALSE)
for (nm in c("PCS", "FDR", "TPR")) {
  v <- t(vapply(ids, dif, c(d = 0, se = 0, z = 0), metric = tolower(nm)))
  base <- ids == "baseline"
  dtab[[paste0("d", nm, " (SE)")]] <- ifelse(base, "-", sprintf("%+.1f (%.1f)", v[, "d"], v[, "se"]))
  dtab[[paste0("z", nm)]] <- ifelse(base, "-", sprintf("%+.1f", v[, "z"]))
}
cat("\nDifference from the baseline in percentage points, with its SE and z (difference / SE).\n")
show_groups(dtab)

# Start-up phase (#51): each configuration with startup = "diagonal" next to
# the same configuration without it (the id without "_startup"). Percent
# metrics and the difference in PCS are in percent (points). The difference
# is with minus without; the runs are independent, so SE and z are as above.
# nStartup: mean number of cohorts in the phase, from mtd_startup_phase.csv.
with_ids <- ids[cfgs$startup[match(ids, cfgs$id)] == "diagonal" & cfgs$pick[match(ids, cfgs$id)] == "interval"]
if (length(with_ids)) {
  phase <- utils::read.csv("results/mtd_startup_phase.csv")
  pairs <- data.frame(without = sub("_startup$", "", with_ids), with = with_ids, stringsAsFactors = FALSE)
  stopifnot(all(pairs$without %in% ids))
  one <- function(id, partner = NULL) {
    d <- if (is.null(partner)) c(d = NA, se = NA, z = NA) else dif(id, "pcs", partner)
    ph <- phase$mean_n_startup[phase$id == id]
    data.frame(pair = pairs$without[pairs$with == id | pairs$without == id][1],
               startup = if (is.null(partner)) "without" else "with",
               PCS = cell(id, "pcs", 100), FDR = cell(id, "fdr", 100), TPR = cell(id, "tpr", 100),
               `%above` = cell(id, "pct_patients_above"), `%band` = cell(id, "pct_patients_in_band"),
               `%below` = cell(id, "pct_patients_below"), cohToContour = cell(id, "first_contour_cohort"),
               DLTs = cell(id, "n_dlt"), nStartup = if (length(ph)) sprintf("%.2f", ph) else "-",
               `dPCS (SE)` = if (is.null(partner)) "-" else sprintf("%+.1f (%.1f)", d[["d"]], d[["se"]]),
               zPCS = if (is.null(partner)) "-" else sprintf("%+.1f", d[["z"]]),
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  st <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    rbind(one(pairs$without[i]), one(pairs$with[i], pairs$without[i]))
  }))
  cat("\n== Start-up phase: without vs with the diagonal start-up phase ==\n")
  cat("pair: the configuration without the phase. Estimate (Monte Carlo SE). PCS, FDR, TPR in percent.\n")
  cat("dPCS: with minus without, in percentage points.\n")
  cat("nStartup: mean cohorts in the phase. The phase has at most 5 cohorts, and ends at the first toxicity.\n")
  print(st, row.names = FALSE, right = FALSE)
  cat(sprintf("Trials whose phase ended by a toxicity (percent): %.1f to %.1f.\n",
              min(phase$pct_phase_ended_by_toxicity), max(phase$pct_phase_ended_by_toxicity)))
}

# Isotonic pick (#54): each configuration with pick = "isotonic" next to its
# interval-pick partner (the id without "_iso": the row with within_set =
# FALSE where there is one, else the within-set row, so the partner may differ
# in within_set; see wset). Percent metrics. The pick classes are computed
# from the stored pct_pick of each combination, so they exist for the partners
# too: over = pick at true toxicity >= 0.45, under = <= 0.15. dPCS: isotonic
# minus partner, in points, with its SE and z (independent runs).
iso_ids <- ids[cfgs$pick[match(ids, cfgs$id)] == "isotonic"]
if (length(iso_ids)) {
  pick_pct <- function(id) {
    v <- res[res$id == id & res$metric == "pct_pick", ]
    stopifnot(identical(v$combo, combo_labels(mtd_grid())))
    c(over = sum(v$estimate[c(mtd_truth) >= 0.45 - 1e-8]), under = sum(v$estimate[c(mtd_truth) <= 0.15 + 1e-8]))
  }
  iso_row <- function(id, partner = NULL) {
    d <- if (is.null(partner)) c(d = NA, se = NA, z = NA) else dif(id, "pcs", partner)
    pk <- pick_pct(id)
    data.frame(pair = sub("_iso$", "", if (is.null(partner)) id else partner),
               pick = if (is.null(partner)) "interval" else "isotonic",
               wset = cfgs$within_set[cfgs$id == id],
               PCS = cell(id, "pcs", 100), FDR = cell(id, "fdr", 100), TPR = cell(id, "tpr", 100),
               `%above` = cell(id, "pct_patients_above"), `%band` = cell(id, "pct_patients_in_band"),
               `%below` = cell(id, "pct_patients_below"), `%pickOver` = sprintf("%.1f", pk[["over"]]),
               `%pickUnder` = sprintf("%.1f", pk[["under"]]),
               cohToContour = cell(id, "first_contour_cohort"), DLTs = cell(id, "n_dlt"),
               `dPCS (SE)` = if (is.null(partner)) "-" else sprintf("%+.1f (%.1f)", d[["d"]], d[["se"]]),
               zPCS = if (is.null(partner)) "-" else sprintf("%+.1f", d[["z"]]),
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  iso_pairs <- sub("_iso$", "", iso_ids)
  stopifnot(all(iso_pairs %in% ids))
  io <- do.call(rbind, lapply(seq_along(iso_ids), function(i) {
    rbind(iso_row(iso_pairs[i]), iso_row(iso_ids[i], iso_pairs[i]))
  }))
  cat("\n== Isotonic pick: Riviere's interval pick vs the isotonic pick ==\n")
  cat("pair: the interval-pick configuration. wset: within_set of the row. Estimate (Monte Carlo SE);\n")
  cat("PCS, FDR, TPR in percent. %pickOver, %pickUnder: picks at true toxicity >= 0.45 and <= 0.15.\n")
  cat("dPCS: isotonic minus interval, in percentage points.\n")
  print(io, row.names = FALSE, right = FALSE)
  cls <- utils::read.csv("results/mtd_pick_classes.csv")
  stopifnot(setequal(cls$id, iso_ids))
  chk <- vapply(cls$id, function(i) {
    stored <- unlist(cls[cls$id == i, c("pct_pick_over", "pct_pick_under")])
    max(abs(unlist(pick_pct(i)) - stored))
  }, 0)
  stopifnot(all(chk < 1e-4))
}

# All configurations sorted by PCS, best first.
pcs_all <- vapply(ids, function(i) get(i, "pcs")$estimate, 0)
best <- tab[order(-pcs_all), c("id", "shift", "eps", "PCS", "FDR", "TPR", "%above", "%band")]
cat("\nBest so far: all configurations sorted by PCS (percent; %above and %band are the share of patients\n")
cat("at toxicity >= 0.40 and between 0.20 and 0.40).\n")
print(best, row.names = FALSE, right = FALSE)
ref <- utils::read.csv("results/riviere2014_scenario1_pcs.csv")
cat(sprintf("Riviere et al. (2014), published PCS of seven designs: %.1f to %.1f (not targets).\n",
            min(ref$pcs), max(ref$pcs)))

grid_of <- function(metric, id) {
  v <- res[res$id == id & res$metric == metric, ]
  m <- matrix(v$estimate, nrow(mtd_truth), ncol(mtd_truth),
              dimnames = list(paste0("A", seq_len(nrow(mtd_truth))), paste0("B", seq_len(ncol(mtd_truth)))))
  stopifnot(identical(v$combo, combo_labels(mtd_grid())))
  m
}
cat(sprintf("\nTrue toxicity (rows = drug A level, columns = drug B level):\n"))
print(mtd_truth, digits = 2)
cat(sprintf("\n%s: %% of trials with the combination in the final set:\n", show_id))
print(round(grid_of("pct_selected", show_id), 1))
cat(sprintf("\n%s: %% of trials whose single pick is the combination (true MTDs: sum = PCS):\n", show_id))
print(round(grid_of("pct_pick", show_id), 1))
cat(sprintf("Sum over the true MTDs: %.1f; PCS: %.1f\n",
            sum(grid_of("pct_pick", show_id)[mtd_is_mtd(c(mtd_truth))]), 100 * get(show_id, "pcs")$estimate))
cat(sprintf("\n%s: mean patients per combination:\n", show_id))
print(round(grid_of("mean_n", show_id), 1))

cat("\nRiviere et al. (2014), Table II, scenario 1: published PCS (%, 2000 replications):\n")
print(ref[c("design", "pcs")], row.names = FALSE)
cat(sprintf("Range %.1f to %.1f. For reference, not targets: the designs differ.\n", min(ref$pcs), max(ref$pcs)))
