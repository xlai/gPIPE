# Mander & Sweeting (2015) Tables III and IV: scenarios, designs and the
# summary shared by the gPIPE run (run_gpipe.R) and the reference run with
# pipe.design 0.5.1 (run_pipe_design.R).
#
# Every matrix is rows = drug A level, columns = drug B level. The paper
# prints drug B in rows (Table I: highest level on top; Table II: level 1 on
# top) and drug A in columns.

from_v5 <- function(x) t(matrix(x, 6, 6, byrow = TRUE)[6:1, ])
from_v6 <- function(x) t(matrix(x, 4, 4, byrow = TRUE)) / 100

# Table I: 6 x 6, theta = 0.30 (Table III).
v5_scen <- list(
  "1" = from_v5(c(.23, .27, .30, .36, .44, .49, .18, .21, .26, .30, .39, .44, .11, .14, .18, .23, .30, .36,
                  .06, .09, .14, .18, .26, .30, .03, .05, .09, .13, .21, .27, .02, .03, .06, .10, .18, .23)),
  "2" = from_v5(c(.23, .30, .45, .50, .55, .60, .18, .21, .30, .45, .50, .55, .11, .14, .18, .30, .45, .50,
                  .06, .09, .14, .18, .30, .45, .03, .05, .09, .13, .21, .30, .02, .03, .06, .10, .18, .23)),
  "3" = from_v5(c(.30, .38, .48, .58, .68, .78, .22, .30, .40, .50, .60, .70, .17, .25, .35, .45, .50, .60,
                  .12, .20, .30, .40, .45, .55, .06, .14, .24, .34, .39, .49, .02, .10, .20, .30, .35, .45)),
  "4" = from_v5(c(.265, .295, .325, .355, .385, .415, .250, .280, .310, .340, .370, .400,
                  .235, .265, .295, .325, .355, .385, .220, .250, .280, .310, .340, .370,
                  .205, .235, .265, .295, .325, .355, .190, .220, .250, .280, .310, .340))
)

# Table II: 4 x 4, theta = 0.20 (Table IV).
v6_scen <- list(
  A = from_v6(c(4, 8, 12, 16, 10, 14, 18, 22, 16, 20, 24, 28, 22, 26, 30, 34)),
  B = from_v6(c(2, 4, 6, 8, 5, 7, 9, 11, 8, 10, 12, 14, 11, 13, 15, 17)),
  C = from_v6(c(10, 20, 30, 40, 25, 35, 45, 55, 40, 50, 60, 70, 55, 65, 75, 85)),
  D = from_v6(c(44, 48, 52, 56, 50, 54, 58, 62, 56, 60, 64, 68, 62, 66, 70, 74)),
  E = from_v6(c(8, 18, 28, 29, 9, 19, 29, 30, 10, 20, 30, 31, 11, 21, 31, 41)),
  F = from_v6(c(12, 13, 14, 15, 16, 18, 20, 22, 44, 45, 46, 47, 50, 52, 54, 55)),
  G = from_v6(c(1, 2, 3, 4, 4, 10, 15, 20, 6, 15, 30, 45, 10, 30, 50, 80))
)

# One row per simulated configuration. `selection` uses gPIPE's names:
# "min_S" is pipe.design's strategy "ss", "inverse_S" is "ss-random".
# The weak prior has strength 1 / K per dose and its medians are scenario 1
# (Table III) or scenario A (Table IV).
ms2015_configs <- function() {
  v5 <- expand.grid(admis = c("closest", "adjacent"), selection = c("min_S", "inverse_S"),
                    scenario = names(v5_scen), stringsAsFactors = FALSE)
  v6 <- data.frame(admis = "closest", selection = "min_S", scenario = names(v6_scen),
                   stringsAsFactors = FALSE)
  rbind(cbind(table = "III", theta = 0.3, N = 40, cohort = 2, v5),
        cbind(table = "IV", theta = 0.2, N = 50, cohort = 1, v6))
}
ms2015_truth <- function(cfg) if (cfg$table == "III") v5_scen[[cfg$scenario]] else v6_scen[[cfg$scenario]]
ms2015_prior_median <- function(cfg) if (cfg$table == "III") v5_scen[["1"]] else v6_scen[["A"]]

# True-toxicity bands. Table III: toxicity in percent, [0,15), [15,25),
# [25,35), [35,46), 46 and over. Table IV: distance from theta in percentage
# points: 0, 1 to 10, more than 10 (either side).
ms2015_band <- function(table, truth, theta) {
  if (table == "III") {
    cut(round(100 * truth, 1), c(-Inf, 14.99, 24.99, 34.99, 45.99, Inf),
        labels = c("0-14", "15-24", "25-34", "35-45", "46+"))
  } else {
    d <- abs(round(100 * (truth - theta)))
    factor(ifelse(d == 0, "at", ifelse(d <= 10, "1-10", ">10")), levels = c("at", "1-10", ">10"))
  }
}

# Estimate and Monte Carlo standard error of a ratio of totals,
# sum(x) / sum(d), over trials (delta method).
ratio_se <- function(x, d) {
  m <- length(x)
  est <- if (sum(d) > 0) sum(x) / sum(d) else NA_real_
  se <- if (sum(d) > 0) sqrt(m * stats::var(x - est * d)) / sum(d) else NA_real_
  c(estimate = est, se = se)
}

# M&S's summary of one configuration, as pipe.design computes it.
#   n, rec:  trials x combinations: patients treated, and 1 where the
#            combination is a recommended phase II dose;
#   stopped: trials stopped by the safety rule;
#   truth:   true toxicity per combination (same order as the columns).
# Experimentation %: patients per band over all patients. Table IV counts
# the patients not treated in stopped trials as "None" (denominator N per
# trial). Recommendation %: recommended doses per band, pooled over trials.
# Table IV counts each stopped trial once as "None".
# Returns one row per cell with the estimate (percent, or doses for
# mean_n_rec) and its Monte Carlo standard error.
ms2015_summary <- function(cfg, n, rec, stopped, truth) {
  band <- ms2015_band(cfg$table, c(truth), cfg$theta)
  by_band <- function(x) vapply(levels(band), function(b) rowSums(x[, band == b, drop = FALSE]), numeric(nrow(x)))
  n_b <- by_band(n)
  rec_b <- by_band(rec)
  if (cfg$table == "IV") {
    n_b <- cbind(n_b, None = cfg$N - rowSums(n))
    rec_b <- cbind(rec_b, None = as.numeric(stopped))
  }
  cell <- function(metric, x) {
    out <- vapply(colnames(x), function(b) 100 * ratio_se(x[, b], rowSums(x)), numeric(2))
    data.frame(metric = metric, band = colnames(x), estimate = out[1, ], se = out[2, ])
  }
  n_rec <- rowSums(rec)
  out <- rbind(cell("exp", n_b), cell("rec", rec_b),
               data.frame(metric = "mean_n_rec", band = "", estimate = mean(n_rec),
                          se = stats::sd(n_rec) / sqrt(length(n_rec))))
  cbind(cfg[c("table", "scenario", "admis", "selection")], nsim = nrow(n), out, row.names = NULL)
}

# Trials per configuration: REPS times the paper's 1000 (III) or 2000 (IV).
ms2015_nsim <- function(cfg, reps) reps * if (cfg$table == "III") 1000 else 2000
