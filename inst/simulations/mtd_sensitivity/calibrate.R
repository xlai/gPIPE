# Calibrate epsilon for the baseline design and write every epsilon tried to
# a file in results/.
# Usage: Rscript calibrate.R [seed = 1] [alpha = 0.05] [out = results/calibration_baseline.csv]
#                            [band = 0.1] [shift = 0]
# `band` is the distance of the simulated truths from theta, and `shift` the
# shift of the design's decision threshold (theta + shift). The summaries of
# run.R use band 0.1 whatever the calibration band.
# For example: Rscript calibrate.R 1 0.10 results/calibration_alpha_0.10.csv
#              Rscript calibrate.R 1 0.05 results/calibration_shift_0.05.csv 0.05 0.05
# Run from this folder. One process (calibrate_eps() is serial).
library(gPIPE)
source("common.R")
args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args) >= 1) as.integer(args[1]) else 1L
alpha <- if (length(args) >= 2) as.numeric(args[2]) else 0.05
out <- if (length(args) >= 3) args[3] else "results/calibration_baseline.csv"
band <- if (length(args) >= 4) as.numeric(args[4]) else mtd_band
shift <- if (length(args) >= 5) as.numeric(args[5]) else 0

# Settings. Each evaluation runs nsim trials of 60 patients in each scenario.
scenarios <- mtd_upsets()  # all 56 monotone configurations, so no ensemble is drawn
n_scenarios <- length(scenarios)
nsim <- 600L
tol <- 0.005
lower <- 0.01
upper <- 0.5

base <- mtd_configs()
base <- base[base$id == "baseline", ]
base$shift <- shift
design <- mtd_design(base, eps = upper)  # eps is replaced by the calibration
t0 <- Sys.time()
cal <- calibrate_eps(design, scenarios = scenarios, nsim = nsim, alpha = alpha, band = band,
                     lower = lower, upper = upper, tol = tol, seed = seed)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
print(cal)
cat(sprintf("Run time: %.0f seconds\n", secs))

tr <- cal@trace
res <- data.frame(eps = tr$eps, fdr = tr$fdr, fdr_se = tr$fdr_se, tpr = tr$tpr,
                  selected = !is.na(cal@eps) & abs(tr$eps - cal@eps) < 1e-12,
                  alpha = alpha, band = band, shift = shift, n_scenarios = n_scenarios, nsim = nsim,
                  tol = tol, lower = lower, upper = upper, seed = seed, seconds = round(secs))
utils::write.csv(res, out, row.names = FALSE)
