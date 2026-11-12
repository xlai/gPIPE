# Print the stored results: one row per configuration with the main metrics,
# the selection grid of one configuration, and Riviere's reference line.
# Usage: Rscript report.R [id = baseline]
# Run from this folder.
library(gPIPE)
source("common.R")
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
cat("Estimate (Monte Carlo SE). PCS, FDR, TPR and P.setFalse are in percent.\n\n")
print(tab, row.names = FALSE, right = FALSE)

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

ref <- utils::read.csv("results/riviere2014_scenario1_pcs.csv")
cat("\nRiviere et al. (2014), Table II, scenario 1: published PCS (%, 2000 replications):\n")
print(ref[c("design", "pcs")], row.names = FALSE)
cat(sprintf("Range %.1f to %.1f. For reference, not targets: the designs differ.\n", min(ref$pcs), max(ref$pcs)))
