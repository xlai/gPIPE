# Stored results of the MTD component sensitivity study (#45, steps 1 and 2).
#
# The results are stored in inst/simulations/mtd_sensitivity/results: the
# calibration of epsilon for the baseline design, and 10000 trials of each
# configuration (baseline, epsilon = 0.5, and the Mander & Sweeting design)
# on Riviere et al.'s scenario 1. The scripts beside them make the runs
# again. This checks the stored files only, so it does not change when the
# package code changes.

mtd_dir <- function() system.file("simulations", "mtd_sensitivity", package = "gPIPE")

test_that("the stored MTD sensitivity results are complete and consistent", {
  skip_on_cran()
  skip_if(mtd_dir() == "", "simulation results are not installed")
  rdir <- file.path(mtd_dir(), "results")
  read <- function(f, ...) utils::read.csv(file.path(rdir, f), ...)
  for (f in c("mtd_sensitivity.csv", "mtd_sensitivity_runs.csv", "calibration_baseline.csv",
              "riviere2014_scenario1_pcs.csv")) {
    expect_true(file.exists(file.path(rdir, f)), info = f)
  }
  res <- read("mtd_sensitivity.csv", colClasses = c(combo = "character"))
  runs <- read("mtd_sensitivity_runs.csv")
  cal <- read("calibration_baseline.csv")
  ref <- read("riviere2014_scenario1_pcs.csv")

  ids <- c("baseline", "eps_0.5", "ref_ms")
  expect_true(all(ids %in% res$id))
  expect_true(all(ids %in% runs$id))
  expect_equal(runs$nsim[match(ids, runs$id)], rep(10000, 3))

  # Every metric is there, for every configuration, and the combinations are
  # the 15 of the 5 x 3 grid.
  overall <- c("fdr", "tpr", "pcs", "pct_patients_above", "pct_patients_in_band",
               "pct_patients_below", "p_set_false", "mean_recommended", "pct_stopped",
               "mean_patients", "first_contour_cohort", "pct_never_contour", "n_distinct",
               "pct_skip", "n_dlt", "pct_set_has_mtd", "n_mtd_tried")
  by_combo <- c("pct_selected", "pct_pick", "mean_n")
  labels <- paste0("A", rep(1:5, 3), ".B", rep(1:3, each = 5))
  for (id in ids) {
    r <- res[res$id == id, ]
    expect_setequal(r$metric[r$combo == ""], overall)
    for (m in by_combo) expect_identical(r$combo[r$metric == m], labels)
    expect_equal(nrow(r), length(overall) + 3 * 15, info = id)
  }

  # Estimates and standard errors are in range. Proportions lie in [0, 1],
  # percentages in [0, 100], and no estimate or SE is missing.
  expect_false(anyNA(res$estimate))
  expect_true(all(res$mc_se >= 0 | is.na(res$mc_se)))
  is_prop <- res$metric %in% c("fdr", "tpr", "pcs", "p_set_false")
  is_pct <- grepl("^pct_", res$metric)
  expect_true(all(res$estimate[is_prop] >= 0 & res$estimate[is_prop] <= 1))
  expect_true(all(res$estimate[is_pct] >= 0 & res$estimate[is_pct] <= 100))
  expect_true(all(res$estimate[res$metric %in% c("mean_patients", "mean_n")] <= 60))
  # A trial has 20 cohorts, so the first cohort on the contour is 1 to 20.
  fc <- res$estimate[res$metric == "first_contour_cohort"]
  expect_true(all(fc >= 1 & fc <= 20))
  expect_true(all(res$estimate[res$metric == "n_mtd_tried"] <= 3))
  expect_true(all(res$estimate[res$metric == "n_distinct"] <= 15))
  # Patients per combination add up to the mean number enrolled.
  for (id in ids) {
    r <- res[res$id == id, ]
    expect_equal(sum(r$estimate[r$metric == "mean_n"]), r$estimate[r$metric == "mean_patients"],
                 tolerance = 1e-5, info = id)
  }

  # PCS is the sum of the single-pick percentages over the true MTDs: the
  # combinations with toxicity 0.30, A4.B1, A3.B2 and A2.B3.
  mtd <- c("A4.B1", "A3.B2", "A2.B3")
  for (id in ids) {
    r <- res[res$id == id, ]
    pick <- r[r$metric == "pct_pick", ]
    expect_equal(sum(pick$estimate[pick$combo %in% mtd]) / 100, r$estimate[r$metric == "pcs"],
                 tolerance = 1e-5, info = id)
    expect_lte(sum(pick$estimate), 100 + 1e-6)
  }

  # The calibrated epsilon is the one the baseline row used, and it is the
  # selected row of the calibration path, whose false discovery rate is
  # within alpha.
  expect_equal(sum(cal$selected), 1)
  eps_cal <- cal$eps[cal$selected]
  expect_equal(runs$eps[runs$id == "baseline"], eps_cal)
  expect_lte(cal$fdr[cal$selected], cal$alpha[1])
  expect_true(is.na(runs$eps[runs$id == "ref_ms"]))
  expect_equal(runs$eps[runs$id == "eps_0.5"], 0.5)

  # Riviere's published PCS for scenario 1 (percent): seven designs.
  expect_equal(nrow(ref), 7)
  expect_true(all(ref$pcs > 60 & ref$pcs < 80))
})
