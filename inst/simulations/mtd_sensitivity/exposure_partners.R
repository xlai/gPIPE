# Exposure to over-toxic doses (#57) for the partner rows of the closest rows.
# Re-simulates each partner with its stored seed and settings (10000 trials),
# checks that the PCS equals the stored PCS exactly (so the seeds reproduce),
# and adds its four exposure numbers to results/mtd_exposure.csv. The stored
# results of the partners (mtd_sensitivity*.csv) are not touched.
# Usage: Rscript exposure_partners.R [cores = 4] [ids ...]
# Run from this folder. With no ids, the partners of the closest rows.
library(gPIPE)
source("common.R")
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 4L
cfgs <- mtd_configs()
partners <- unname(mtd_closest_partners())
ids <- if (length(args) >= 2) args[-1] else partners
stopifnot(all(ids %in% cfgs$id))
dir <- "results"
res <- utils::read.csv(file.path(dir, "mtd_sensitivity.csv"), colClasses = c(combo = "character"))

run_one <- function(id) {
  cfg <- cfgs[cfgs$id == id, ]
  sim <- simulate(mtd_design(cfg, mtd_eps(cfg, dir)), nsim = mtd_nsim, seed = cfg$seed, truth = mtd_truth,
                  keep_trials = TRUE)
  pcs <- mean(sim@summary$pick %in% which(mtd_is_mtd(c(mtd_truth))))
  list(exposure = mtd_exposure(id, sim, mtd_truth), pcs = pcs)
}
out <- parallel::mclapply(ids, run_one, mc.cores = min(cores, length(ids)), mc.preschedule = FALSE)
bad <- !vapply(out, is.list, logical(1))
if (any(bad)) stop("These configurations failed: ", paste(ids[bad], collapse = ", "))
for (i in seq_along(ids)) {
  stored <- res$estimate[res$id == ids[i] & res$metric == "pcs"]
  ok <- identical(signif(out[[i]]$pcs, 7), stored)
  cat(sprintf("%-30s PCS re-simulated %.7g, stored %.7g: %s\n", ids[i], out[[i]]$pcs, stored,
              if (ok) "equal" else "DIFFERENT"))
  if (!ok) stop("The re-simulation of ", ids[i], " does not reproduce its stored PCS.", call. = FALSE)
}
exposure <- do.call(rbind, lapply(out, `[[`, "exposure"))
exposure[-1] <- lapply(exposure[-1], signif, 7)
f <- file.path(dir, "mtd_exposure.csv")
if (file.exists(f)) {
  old <- utils::read.csv(f)
  exposure <- rbind(old[!old$id %in% exposure$id, names(exposure)], exposure)
}
exposure <- exposure[order(match(exposure$id, cfgs$id)), ]
utils::write.csv(exposure, f, row.names = FALSE)
print(exposure, row.names = FALSE)
