# Calibrate epsilon for the baseline design (alpha = 0.05, band = 0.1) and
# write every epsilon tried to results/calibration_baseline.csv.
# Usage: Rscript calibrate.R [seed = 1] [out = results/calibration_baseline.csv]
# Run from this folder. One process (calibrate_eps() is serial).
library(gPIPE)
source("common.R")
args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args) >= 1) as.integer(args[1]) else 1L
out <- if (length(args) >= 2) args[2] else "results/calibration_baseline.csv"

# Settings. Each evaluation runs nsim trials of 60 patients in each scenario.
alpha <- 0.05
scenarios <- mtd_upsets()  # all 56 monotone configurations, so no ensemble is drawn
n_scenarios <- length(scenarios)
nsim <- 600L
tol <- 0.005
lower <- 0.01
upper <- 0.5

base <- mtd_configs()
base <- base[base$id == "baseline", ]
design <- mtd_design(base, eps = upper)  # eps is replaced by the calibration
t0 <- Sys.time()
cal <- calibrate_eps(design, scenarios = scenarios, nsim = nsim, alpha = alpha, band = mtd_band,
                     lower = lower, upper = upper, tol = tol, seed = seed)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
print(cal)
cat(sprintf("Run time: %.0f seconds\n", secs))

tr <- cal@trace
res <- data.frame(eps = tr$eps, fdr = tr$fdr, fdr_se = tr$fdr_se, tpr = tr$tpr,
                  selected = !is.na(cal@eps) & abs(tr$eps - cal@eps) < 1e-12,
                  alpha = alpha, band = mtd_band, n_scenarios = n_scenarios, nsim = nsim,
                  tol = tol, lower = lower, upper = upper, seed = seed, seconds = round(secs))
utils::write.csv(res, out, row.names = FALSE)
