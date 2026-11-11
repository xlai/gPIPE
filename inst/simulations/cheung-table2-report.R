# The 15-cell comparison with Cheung et al. (2022), Table 2.
# Input: one row per design x scenario (see cheung-table2.R).

cheung_table2_published <- data.frame(
  design = c("non-adaptive", "1-EQ", "1-PP", "2-EQ", "2-PP"),
  eps = c(0.0525, 0.05, 0.0625, 0.0425, 0.0525),
  fdr = c(0.046, 0.048, 0.048, 0.046, 0.049),
  tpr = c(0.68, 0.79, 0.78, 0.88, 0.91),
  n_eff = c(63, 103, 111, 66, 91),
  stringsAsFactors = FALSE
)

# Printed unit of each metric, for the rounding error.
cheung_table2_units <- c(fdr = 0.001, tpr = 0.01, n_eff = 1)

# `scen`: design, scenario, n_trials, n_right and the per-scenario mean and
# variance of each metric (fdr_mean, fdr_var, ...). `tpr0`: what to do with
# the all-ineffective scenario, where TPR is 0/0: "exclude", "zero" or "one".
cheung_table2_report <- function(scen, reps_cheung = 40, tpr0 = c("exclude", "zero", "one")) {
  tpr0 <- match.arg(tpr0)
  rows <- lapply(seq_len(nrow(cheung_table2_published)), function(i) {
    pub <- cheung_table2_published[i, ]
    s <- scen[scen$design == pub$design, ]
    lapply(c("fdr", "tpr", "n_eff"), function(m) {
      mu <- s[[paste0(m, "_mean")]]
      v <- s[[paste0(m, "_var")]]
      if (m == "tpr") {
        undefined <- is.na(mu)
        if (tpr0 == "zero") mu[undefined] <- 0
        if (tpr0 == "one") mu[undefined] <- 1
        v[undefined] <- 0
        keep <- !is.na(mu)
        mu <- mu[keep]
        v <- v[keep]
        n <- s$n_trials[keep]
      } else {
        n <- s$n_trials
      }
      S <- length(mu)
      est <- mean(mu)
      se_g <- sqrt(sum(v / n)) / S
      se_c <- sqrt(sum(v / reps_cheung)) / S
      se_r <- cheung_table2_units[[m]] / sqrt(12)
      se <- sqrt(se_g^2 + se_c^2 + se_r^2)
      d <- est - pub[[m]]
      data.frame(design = pub$design, metric = m, published = pub[[m]], gpipe = est,
                 diff = d, se_cheung = se_c, se_gpipe = se_g, se_round = se_r, se = se,
                 z = d / se, margin = 3 * se, pass = abs(d) <= 3 * se)
    })
  })
  out <- do.call(rbind, unlist(rows, recursive = FALSE))
  rownames(out) <- NULL
  out
}
