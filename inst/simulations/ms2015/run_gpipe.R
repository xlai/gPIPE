# Mander & Sweeting (2015) Tables III and IV with gPIPE.
# Usage: Rscript run_gpipe.R [reps = 4] [cores = 2] [out = results/ms2015_gpipe.csv]
# Run from this folder. Trials per configuration: reps x 1000 (Table III) or
# reps x 2000 (Table IV).
library(gPIPE)
source("common.R")
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 4L
cores <- if (length(args) >= 2) as.integer(args[2]) else 2L
out <- if (length(args) >= 3) args[3] else "results/ms2015_gpipe.csv"

cfgs <- ms2015_configs()
run_one <- function(i) {
  cfg <- cfgs[i, ]
  truth <- ms2015_truth(cfg)
  grid <- DoseGrid(levels = c(A = nrow(truth), B = ncol(truth)))
  design <- ms_design(grid, theta = cfg$theta,
                      prior_median = grid_from_array(grid, ms2015_prior_median(cfg)),
                      prior_strength = 1 / length(truth), N = cfg$N, cohort = cfg$cohort,
                      admis = cfg$admis, selection = cfg$selection, safety = 0.8)
  sim <- simulate(design, nsim = ms2015_nsim(cfg, reps), seed = 3100 + i, truth = truth)
  ms2015_summary(cfg, sim@n, sim@final_set, sim@summary$stopped, sim@truth)
}
res <- parallel::mclapply(seq_len(nrow(cfgs)), run_one, mc.cores = cores, mc.preschedule = FALSE)
ms2015_write(res, out)
