# Mander & Sweeting (2015) Tables III and IV with their own package,
# pipe.design 0.5.1. It is archived on CRAN: install ggplot2, gtools and
# xtable, then the source package from
# https://cloud.r-project.org/src/contrib/Archive/pipe.design/pipe.design_0.5.1.tar.gz
# Usage: Rscript run_pipe_design.R [reps = 4] [cores = 2] [out = results/ms2015_pipe_design.csv]
# Run from this folder. Slow: about 135 s per 1000 trials of Table III.
suppressMessages(library(pipe.design))
source("common.R")
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 4L
cores <- if (length(args) >= 2) as.integer(args[2]) else 2L
out <- if (length(args) >= 3) args[3] else "results/ms2015_pipe_design.csv"

cfgs <- ms2015_configs()
run_one <- function(i) {
  cfg <- cfgs[i, ]
  truth <- ms2015_truth(cfg)
  S <- ms2015_nsim(cfg, reps)
  # pipe.design prints its progress; capture.output() keeps it quiet.
  utils::capture.output(
    fit <- pipe.design(N = cfg$N, S = S, c = cfg$cohort, theta = cfg$theta, pi = truth,
                       prior.med = ms2015_prior_median(cfg),
                       prior.ss = matrix(1 / length(truth), nrow(truth), ncol(truth)),
                       strategy = c(min_S = "ss", inverse_S = "ss-random")[[cfg$selection]],
                       admis = cfg$admis, constraint = "neighbouring", epsilon = 0.8,
                       seed = 4100 + i)
  )
  # Trials x combinations, in column-major order of the A x B matrix.
  n <- t(vapply(fit$n.sim, c, numeric(length(truth))))
  rec <- t(vapply(fit$rpII.list, function(rp) {
    m <- matrix(0, nrow(truth), ncol(truth))
    for (h in seq_len(ncol(rp))) m[rp[1, h], rp[2, h]] <- 1
    c(m)
  }, numeric(length(truth))))
  # A stopped trial has dose 0 in its last cohort.
  stopped <- apply(fit$rec.i.sim[, , 1, drop = FALSE] == 0, 1, any)
  ms2015_summary(cfg, n, rec, stopped, c(truth))
}
res <- parallel::mclapply(seq_len(nrow(cfgs)), run_one, mc.cores = cores, mc.preschedule = FALSE)
res <- do.call(rbind, res)
res$estimate <- signif(res$estimate, 6)
res$se <- signif(res$se, 4)
write.csv(res, out, row.names = FALSE)
