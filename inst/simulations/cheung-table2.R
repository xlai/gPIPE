# nolint start: object_usage_linter. Objects from cheung-table2-report.R are sourced at run time.
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
source(file.path(here, "cheung-table2-report.R"), local = TRUE)

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
  one <- function(x) {
    c(mean = if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE),
      var = if (sum(!is.na(x)) > 1) stats::var(x, na.rm = TRUE) else NA_real_)
  }
  f <- one(m$fdr)
  t <- one(m$tpr)
  e <- one(m$n_eff)
  data.frame(design = design, eps = eps, scenario = scenario, n_trials = length(m$fdr),
             n_right = sum(s == 1), fdr_mean = f[["mean"]], fdr_var = f[["var"]],
             tpr_mean = t[["mean"]], tpr_var = t[["var"]],
             n_eff_mean = e[["mean"]], n_eff_var = e[["var"]])
}

# Posterior variants for the normal outcome (README, setting 1).
#  "pooled": gPIPE's model, one pooled variance (outcome_normal()).
#  "cell":   Cheung's Sec. 2.2 model read per combination: each cell has its own
#            variance, flat prior on the mean and 1/sigma on the SD, so
#            P(mu > 0) = t_{n-1}(ybar / (s / sqrt(n))). With fewer than 2
#            patients the variance is unknown, so the known SD is used.
cell_exceedance <- function(theta, sigma0) {
  function(n, s1, s2) {
    p <- rep(0.5, length(n))
    tried <- n > 0
    ybar <- s1 / pmax(n, 1)
    ss <- pmax(s2 - s1^2 / pmax(n, 1), 0)
    ok <- n >= 2 & ss > 0
    z <- ifelse(ok, stats::pt((ybar - theta) / (sqrt(ss / pmax(n - 1, 1)) / sqrt(pmax(n, 1))), pmax(n - 1, 1)),
                stats::pnorm((ybar - theta) / (sigma0 / sqrt(pmax(n, 1)))))
    p[tried] <- z[tried]
    p
  }
}

# The plain design list the trial loop runs on, with the chosen posterior.
design_with <- function(des, posterior) {
  d <- gPIPE:::design_list(des)
  if (posterior == "cell") d$model$exceedance <- cell_exceedance(des@theta, SD_TRUE)
  d
}

truth_for <- function(s) list(mean = MU_EFFECTIVE * s, sd = SD_TRUE)

# Adaptive designs: full trials through simulate(), for scenario `i`.
run_adaptive <- function(des, eps, name, scenarios, i, reps, seed0, posterior = "pooled") {
  des <- S7::set_props(des, eps = eps)
  s <- grid_from_array(des@grid, scenarios[[i]])
  if (posterior == "pooled") {
    sim <- simulate(des, nsim = reps, seed = seed0 + i, truth = truth_for(scenarios[[i]]))
    sim <- list(n = sim@n, gamma = sim@gamma, final_set = sim@final_set)
  } else {
    tr <- truth_for(scenarios[[i]])
    tr$mean <- c(grid_from_array(des@grid, tr$mean))
    sim <- gPIPE:::pipe_simulate(design_with(des, posterior), tr, reps, seed = seed0 + i, keep_trials = FALSE)
  }
  scen_row(name, eps, i, s, cheung_metrics(s, sim$n, sim$gamma, sim$final_set))
}

# Non-adaptive design: 5 patients on every combination; the estimate at each
# epsilon in `eps` uses the same data.
run_nonadaptive <- function(scenarios, reps, eps, seed0, grid = cheung_grid(), posterior = "pooled") {
  des <- cheung_design(grid, theta = 0, eps = eps[1], target = "MED", outcome = "normal",
                       sigma0 = SD_TRUE, N = N_TRIALS)
  d <- design_with(des, posterior)
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
stage <- if (isTRUE(getOption("gpipe.sim.source_only"))) "none" else if (length(args)) args[1] else "main"
out_dir <- file.path(here, "results")
dir.create(out_dir, showWarnings = FALSE)

posterior <- if (length(args) > 2) args[3] else "pooled"
stopifnot(posterior %in% c("pooled", "cell"))
suffix <- if (posterior == "pooled") "" else paste0("_", posterior)
save_scen <- function(x) {
  utils::write.csv(x, file.path(out_dir, paste0("cheung_table2_scenarios", suffix, ".csv")), row.names = FALSE)
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
  na <- run_nonadaptive(scenarios, reps, pub$eps[pub$design == "non-adaptive"], seed0 = 1000,
                        posterior = posterior)
  res <- do.call(rbind, unlist(na, recursive = FALSE))
  save_scen(res)
  tick("non-adaptive")
  show(res[res$design == "non-adaptive", ])
  for (nm in names(designs)) {
    rows <- parallel::mclapply(seq_along(scenarios), function(i) {
      run_adaptive(designs[[nm]], pub$eps[pub$design == nm], nm, scenarios, i, reps,
                   seed0 = 100000L * match(nm, names(designs)), posterior = posterior)
    }, mc.cores = cores())
    res <- rbind(res, do.call(rbind, rows))
    save_scen(res)
    tick(nm)
    show(res[res$design == nm, ])
  }
  utils::write.csv(cheung_table2_report(res), file.path(out_dir, paste0("cheung_table2_cells", suffix, ".csv")),
                   row.names = FALSE)
  show(res)
}

# Soft check: the epsilon that keeps FDR <= 0.05 (reported, not a gate).
#  - non-adaptive: FDR and TPR on Cheung's grid 0.01, 0.0125, ..., 0.10, both posteriors;
#  - adaptive: calibrate_eps() on the 252 scenarios (pooled posterior), 40 trials per scenario (as Cheung).
if (stage == "calib") {
  scenarios <- cheung_upsets()
  grid_eps <- seq(0.01, 0.10, by = 0.0025)
  for (post in c("pooled", "cell")) {
    na <- run_nonadaptive(scenarios, 400L, grid_eps, seed0 = 1000, posterior = post)
    df <- do.call(rbind, unlist(na, recursive = FALSE))
    tab <- do.call(rbind, lapply(split(df, df$eps), function(x) {
      data.frame(posterior = post, eps = x$eps[1], fdr = mean(x$fdr_mean),
                 fdr_se = sqrt(sum(x$fdr_var / x$n_trials)) / nrow(x), tpr = mean(x$tpr_mean, na.rm = TRUE))
    }))
    utils::write.csv(tab, file.path(out_dir, paste0("cheung_table2_calibration_nonadaptive_", post, ".csv")),
                     row.names = FALSE)
    cat(post, ": largest grid eps with FDR <= 0.05 =", max(tab$eps[tab$fdr <= 0.05]), "\n")
  }
  designs <- cheung_designs()
  cal <- parallel::mclapply(names(designs), function(nm) {
    r <- calibrate_eps(designs[[nm]], scenarios = scenarios, nsim = 40, band = MU_EFFECTIVE, sd = SD_TRUE,
                       lower = 0.01, upper = 0.10, tol = 0.0025, seed = 7)
    data.frame(design = nm, eps = r@eps, fdr = r@fdr, fdr_se = r@fdr_se, tpr = r@tpr,
               printed_eps = cheung_table2_published$eps[cheung_table2_published$design == nm])
  }, mc.cores = cores())
  cal <- do.call(rbind, cal)
  utils::write.csv(cal, file.path(out_dir, "cheung_table2_calibration_adaptive.csv"), row.names = FALSE)
  print(cal)
}
# nolint end
