# Cheung, Chandereng and Diaz (2022), Table 2: five PIPE designs on the 252
# plateau scenarios of the sedentary break study. See README.md for settings.
#
#   Rscript inst/simulations/cheung-table2.R main [reps]   # default 400 per scenario
#   Rscript inst/simulations/cheung-table2.R calib         # soft epsilon check
#
# Results go to inst/simulations/results/. Seeds are fixed per design and
# scenario, so results do not depend on the number of cores.

suppressPackageStartupMessages(library(gPIPE))
if (!exists("here")) here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
source(file.path(here, "cheung-table2-report.R"))

K_LEVELS <- 5L
N_TRIALS <- 125L
SD_TRUE <- 0.16
MU_EFFECTIVE <- 0.15

# Every 0/1 up-set of a k x k grid (252 for k = 5), as matrices with rows =
# drug A levels. Row j is effective from column t_j on; t never increases.
cheung_upsets <- function(k = K_LEVELS) {
  t <- as.matrix(expand.grid(rep(list(seq_len(k + 1L)), k)))
  t <- t[apply(t, 1, function(r) all(diff(r) <= 0)), , drop = FALSE]
  lapply(seq_len(nrow(t)), function(i) {
    outer(seq_len(k), seq_len(k), function(j, kk) as.integer(kk >= t[i, j]))
  })
}

cheung_grid <- function() DoseGrid(levels = c(A = K_LEVELS, B = K_LEVELS))

cheung_designs <- function(grid = cheung_grid()) {
  mk <- function(eps, taper, selection) {
    cheung_design(grid, theta = 0, eps = eps, target = "MED", outcome = "normal",
                  sigma0 = SD_TRUE, N = N_TRIALS, taper = taper, selection = selection)
  }
  list(
    "1-EQ" = mk(0.05, FALSE, "equal"), "1-PP" = mk(0.0625, FALSE, "prob"),
    "2-EQ" = mk(0.0425, TRUE, "equal"), "2-PP" = mk(0.0525, TRUE, "prob")
  )
}

# Per-trial FDR (eq. 7), TPR (eq. 8) and n_eff, from the final estimate. `s`
# is the true gamma (0/1 over the combinations); `n`, `gamma`, `set` are
# trials x combinations. TPR is NA when no combination is effective.
cheung_metrics <- function(s, n, gamma, set) {
  right <- s == 1
  list(
    fdr = rowSums(set[, !right, drop = FALSE]) / pmax(1, rowSums(set)),
    tpr = if (any(right)) rowSums(gamma[, right, drop = FALSE] == 1) / sum(right) else rep(NA_real_, nrow(n)),
    n_eff = rowSums(n[, right, drop = FALSE])
  )
}

scen_row <- function(design, eps, scenario, s, m) {
  one <- function(x) c(mean = if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE),
                       var = if (sum(!is.na(x)) > 1) stats::var(x, na.rm = TRUE) else NA_real_)
  f <- one(m$fdr); t <- one(m$tpr); e <- one(m$n_eff)
  data.frame(design = design, eps = eps, scenario = scenario, n_trials = length(m$fdr),
             n_right = sum(s == 1), fdr_mean = f[["mean"]], fdr_var = f[["var"]],
             tpr_mean = t[["mean"]], tpr_var = t[["var"]],
             n_eff_mean = e[["mean"]], n_eff_var = e[["var"]])
}

truth_for <- function(s) list(mean = MU_EFFECTIVE * s, sd = SD_TRUE)

# Adaptive designs: full trials through simulate(), for scenario `i`.
run_adaptive <- function(des, eps, name, scenarios, i, reps, seed0) {
  des <- S7::set_props(des, eps = eps)
  s <- grid_from_array(des@grid, scenarios[[i]])
  sim <- simulate(des, nsim = reps, seed = seed0 + i, truth = truth_for(scenarios[[i]]))
  scen_row(name, eps, i, s, cheung_metrics(s, sim@n, sim@gamma, sim@final_set))
}

# Non-adaptive design: 5 patients on every combination; the estimate at each
# epsilon in `eps` uses the same data.
run_nonadaptive <- function(scenarios, reps, eps, seed0, grid = cheung_grid()) {
  des <- cheung_design(grid, theta = 0, eps = eps[1], target = "MED", outcome = "normal",
                       sigma0 = SD_TRUE, N = N_TRIALS)
  d <- gPIPE:::design_list(des)
  K <- d$g$K
  per <- N_TRIALS / K
  n <- rep(as.integer(per), K)
  lapply(seq_along(scenarios), function(i) {
    set.seed(seed0 + i)
    s <- grid_from_array(grid, scenarios[[i]])
    mu <- MU_EFFECTIVE * s
    y <- matrix(stats::rnorm(reps * K * per, rep(mu, each = per), SD_TRUE), reps, K * per, byrow = TRUE)
    id <- rep(seq_len(K), each = per)
    s1 <- t(rowsum(t(y), id))
    s2 <- t(rowsum(t(y^2), id))
    lapply(eps, function(e) {
      est <- lapply(seq_len(reps), function(r) {
        p <- d$model$exceedance(n, s1[r, ], s2[r, ])
        gamma <- gPIPE:::trial_estimate(d, p, n, e)
        list(gamma = gamma, set = seq_len(K) %in% gPIPE:::pipe_final_set(d$g, gamma, "MED"))
      })
      gam <- do.call(rbind, lapply(est, `[[`, "gamma"))
      set <- do.call(rbind, lapply(est, `[[`, "set"))
      scen_row("non-adaptive", e, i, s, cheung_metrics(s, matrix(n, reps, K, byrow = TRUE), gam, set))
    })
  })
}

cores <- function() max(1L, as.integer(Sys.getenv("GPIPE_CORES", parallel::detectCores())))

args <- commandArgs(TRUE)
stage <- if (length(args)) args[1] else "main"
out_dir <- file.path(here, "results")
dir.create(out_dir, showWarnings = FALSE)

save_scen <- function(x) {
  utils::write.csv(x, file.path(out_dir, "cheung_table2_scenarios.csv"), row.names = FALSE)
}
show <- function(x) {
  print(cheung_table2_report(x, reps_cheung = 40)[, c("design", "metric", "published", "gpipe", "z", "pass")],
        digits = 3)
}

if (stage == "main") {
  reps <- if (length(args) > 1) as.integer(args[2]) else 400L
  scenarios <- cheung_upsets()
  stopifnot(length(scenarios) == 252L)
  designs <- cheung_designs()
  pub <- cheung_table2_published
  t0 <- Sys.time()
  tick <- function(what) cat(sprintf("%s done (%.0f s)\n", what, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  # Non-adaptive first: no allocation rule, so it checks the posterior and the metrics alone.
  na <- run_nonadaptive(scenarios, reps, pub$eps[pub$design == "non-adaptive"], seed0 = 1000)
  res <- do.call(rbind, unlist(na, recursive = FALSE))
  save_scen(res)
  tick("non-adaptive")
  show(res)
  for (nm in names(designs)) {
    rows <- parallel::mclapply(seq_along(scenarios), function(i) {
      run_adaptive(designs[[nm]], pub$eps[pub$design == nm], nm, scenarios, i, reps,
                   seed0 = 100000L * match(nm, names(designs)))
    }, mc.cores = cores())
    res <- rbind(res, do.call(rbind, rows))
    save_scen(res)
    tick(nm)
    show(res[res$design == nm, ])
  }
  utils::write.csv(cheung_table2_report(res), file.path(out_dir, "cheung_table2_cells.csv"), row.names = FALSE)
  show(res)
}
