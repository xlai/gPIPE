# Simulate the configurations of mtd_configs(): 10000 trials each, then
# summarise them (common.R) and store the summary.
# Usage: Rscript run.R [cores = 4] [ids ...]
# Run from this folder. With no ids, every configuration runs. With ids
# (e.g. Rscript run.R 4 baseline eps_0.5), only those configurations run, and
# only their rows in results/mtd_sensitivity.csv and
# results/mtd_sensitivity_runs.csv are replaced. The rows with startup =
# "diagonal" also write results/mtd_startup_phase.csv, and the rows with pick =
# "isotonic" also write results/mtd_pick_classes.csv.
#
# Seeds: each configuration has its own `seed` in mtd_configs() and runs its
# 10000 trials in one process with set.seed(seed) (Mersenne-Twister). The
# parallel loop is over configurations only, so the numbers do not depend on
# the number of cores, nor on which other configurations run.
library(gPIPE)
source("common.R")
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 4L
cfgs <- mtd_configs()
stopifnot(!anyDuplicated(cfgs$id), !anyDuplicated(cfgs$seed))
ids <- if (length(args) >= 2) args[-1] else cfgs$id
if (!all(ids %in% cfgs$id)) stop("Unknown configuration: ", paste(setdiff(ids, cfgs$id), collapse = ", "))
dir <- "results"
eps <- vapply(ids, function(i) mtd_eps(cfgs[cfgs$id == i, ], dir), numeric(1))  # fails early without calibration

run_one <- function(id) {
  cfg <- cfgs[cfgs$id == id, ]
  t0 <- Sys.time()
  sim <- simulate(mtd_design(cfg, eps[[id]]), nsim = mtd_nsim, seed = cfg$seed, truth = mtd_truth,
                  keep_trials = TRUE)
  res <- mtd_summarise(id, sim)
  phase <- if (cfg$startup == "diagonal") mtd_startup_summarise(id, sim)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  run <- data.frame(id = id, nsim = mtd_nsim, seed = cfg$seed, eps = eps[[id]], seconds = round(secs))
  classes <- if (cfg$pick == "isotonic") mtd_pick_classes(id, sim, mtd_truth)
  list(res = res, phase = phase, classes = classes, run = run)
}
out <- parallel::mclapply(ids, run_one, mc.cores = min(cores, length(ids)), mc.preschedule = FALSE)
bad <- !vapply(out, is.list, logical(1))
if (any(bad)) stop("These configurations failed: ", paste(ids[bad], collapse = ", "))

# Replace the rows of the configurations that ran; keep the others, in the
# order of mtd_configs().
store <- function(new, file) {
  new$estimate <- signif(new$estimate, 7)
  if ("mc_se" %in% names(new)) new$mc_se <- signif(new$mc_se, 5)
  if (file.exists(file)) {
    old <- utils::read.csv(file, colClasses = c(combo = "character"))
    new <- rbind(old[!old$id %in% new$id, names(new)], new)
  }
  new <- new[order(match(new$id, cfgs$id), seq_len(nrow(new))), ]
  utils::write.csv(new, file, row.names = FALSE)
}
res <- do.call(rbind, lapply(out, `[[`, "res"))
run <- do.call(rbind, lapply(out, `[[`, "run"))
store(res, file.path(dir, "mtd_sensitivity.csv"))
f <- file.path(dir, "mtd_sensitivity_runs.csv")
if (file.exists(f)) {
  old <- utils::read.csv(f)
  run <- rbind(old[!old$id %in% run$id, names(run)], run)
}
run <- run[order(match(run$id, cfgs$id)), ]
utils::write.csv(run, f, row.names = FALSE)
print(run, row.names = FALSE)

# Mechanism of the start-up rows (#51): kept apart from the long format so
# that the other rows and their metrics stay as they were. One row per
# configuration with startup = "diagonal".
phase <- do.call(rbind, lapply(out, `[[`, "phase"))
if (!is.null(phase)) {
  phase$mean_n_startup <- signif(phase$mean_n_startup, 7)
  phase$se_n_startup <- signif(phase$se_n_startup, 5)
  phase$pct_phase_ended_by_toxicity <- signif(phase$pct_phase_ended_by_toxicity, 7)
  f <- file.path(dir, "mtd_startup_phase.csv")
  if (file.exists(f)) {
    old <- utils::read.csv(f)
    phase <- rbind(old[!old$id %in% phase$id, names(phase)], phase)
  }
  phase <- phase[order(match(phase$id, cfgs$id)), ]
  utils::write.csv(phase, f, row.names = FALSE)
  print(phase, row.names = FALSE)
}

# Classes of the single pick of the isotonic rows (#54): kept apart from the
# long format so that the stored rows of the other configurations stay as
# they were. One row per configuration with pick = "isotonic".
classes <- do.call(rbind, lapply(out, `[[`, "classes"))
if (!is.null(classes)) {
  classes[-1] <- lapply(classes[-1], signif, 7)
  f <- file.path(dir, "mtd_pick_classes.csv")
  if (file.exists(f)) {
    old <- utils::read.csv(f)
    classes <- rbind(old[!old$id %in% classes$id, names(classes)], classes)
  }
  classes <- classes[order(match(classes$id, cfgs$id)), ]
  utils::write.csv(classes, f, row.names = FALSE)
  print(classes[1:5], row.names = FALSE)
}
