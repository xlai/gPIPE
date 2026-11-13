# Stored results of the MTD component sensitivity study (#45).
#
# The results are stored in inst/simulations/mtd_sensitivity/results: the
# calibrations of epsilon (three values of alpha, four band shifts), and
# 10000 trials of each configuration (the baseline, the one-factor-at-a-time
# rows, the selection, band shift and exploratory combination rows, the rows with
# a start-up phase, the isotonic-pick rows, and the Mander & Sweeting design) on
# Riviere et al.'s scenario 1. The scripts
# beside them make the runs again. This checks the stored files only, so it does not change when the
# package code changes.

mtd_dir <- function() system.file("simulations", "mtd_sensitivity", package = "gPIPE")

test_that("the stored MTD sensitivity results are complete and consistent", {
  skip_on_cran()
  skip_if(mtd_dir() == "", "simulation results are not installed")
  rdir <- file.path(mtd_dir(), "results")
  read <- function(f, ...) utils::read.csv(file.path(rdir, f), ...)
  shifts <- c("0.1", "0.075", "0.05", "0.025")
  for (f in c("mtd_sensitivity.csv", "mtd_sensitivity_runs.csv", "calibration_baseline.csv",
              "calibration_alpha_0.10.csv", "calibration_alpha_0.20.csv",
              "riviere2014_scenario1_pcs.csv", sprintf("calibration_shift_%s.csv", shifts))) {
    expect_true(file.exists(file.path(rdir, f)), info = f)
  }
  res <- read("mtd_sensitivity.csv", colClasses = c(combo = "character"))
  runs <- read("mtd_sensitivity_runs.csv")
  cal <- read("calibration_baseline.csv")
  ref <- read("riviere2014_scenario1_pcs.csv")

  # The start-up rows (#51): each is the row of the same name without the suffix.
  with_startup <- c("baseline_startup", "single_pick_off_startup", "eps_0.3_startup",
                    "eps_0.4_startup", "eps_0.3_pick_all_startup", "eps_0.4_pick_all_startup",
                    "selection_fewest_startup", "selection_inverse_n_startup",
                    "shift_0.1_pick_all_fewest_startup", "shift_0.1_eps_fixed_startup",
                    "shift_0.1_eps_fixed_pick_all_startup")
  # The isotonic-pick rows (#54, step 2): the id of the interval-pick partner
  # plus "_iso", seeds 4552 to 4569.
  with_iso <- c("single_pick_off_iso", "baseline_iso", "single_pick_off_startup_iso",
                "eps_0.3_pick_all_iso", "eps_0.4_pick_all_iso", "eps_0.5_iso",
                "eps_0.4_pick_all_startup_iso", "shift_0.1_eps_fixed_pick_all_iso",
                "shift_0.1_eps_fixed_pick_all_startup_iso", "selection_fewest_iso",
                "selection_inverse_n_iso", "taper_off_iso", "weights_1_iso",
                "prior_strength_1_iso", "cohort_1_iso", "N_120_iso", "N_30_iso", "ref_ms_iso")
  ids <- c("baseline", "eps_0.5", "ref_ms",
           "eps_0.4", "eps_0.3", "eps_0.15", "eps_0.05", "alpha_0.10", "alpha_0.20",
           "taper_off", "weights_1", "prior_strength_1", "prior_strength_4",
           "prior_beta_1_1", "prior_beta_0.5_0.5", "prior_correct", "prior_too_low",
           "prior_too_high", "prior_correct_strength_1", "selection_prob", "cohort_1",
           "N_30", "N_120", "start_top", "single_pick_off",
           "selection_fewest", "selection_inverse_n", "selection_interval",
           "shift_0.1", "shift_0.075", "shift_0.05", "shift_0.025", "shift_0.1_eps_fixed",
           "eps_0.4_pick_all", "eps_0.3_pick_all", "shift_0.1_pick_all",
           "shift_0.1_pick_all_fewest", "shift_0.1_pick_all_interval", "shift_0.05_pick_all",
           "shift_0.1_eps_fixed_pick_all", with_startup, with_iso)
  expect_setequal(res$id, ids)
  expect_setequal(runs$id, ids)
  expect_equal(runs$nsim, rep(10000, length(ids)))
  expect_false(anyDuplicated(runs$seed) > 0)
  expect_equal(runs$seed[match(c("shift_0.1_eps_fixed_pick_all", with_startup), runs$id)], 4540:4551)
  expect_equal(runs$seed[match(with_iso, runs$id)], 4552:4569)

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
    expect_false(anyNA(r$estimate), info = id)
  }

  # Estimates and standard errors are in range. Proportions lie in [0, 1],
  # percentages in [0, 100], and no estimate or SE is missing.
  expect_false(anyNA(res$estimate))
  expect_true(all(res$mc_se >= 0 | is.na(res$mc_se)))
  is_prop <- res$metric %in% c("fdr", "tpr", "pcs", "p_set_false")
  is_pct <- grepl("^pct_", res$metric)
  expect_true(all(res$estimate[is_prop] >= 0 & res$estimate[is_prop] <= 1))
  expect_true(all(res$estimate[is_pct] >= 0 & res$estimate[is_pct] <= 100))
  # N is 30, 60 or 120 patients, in cohorts of 3 (cohort_1: of 1).
  n_max <- c(N_30 = 30, N_120 = 120, N_30_iso = 30, N_120_iso = 120)
  n_cap <- ifelse(res$id %in% names(n_max), n_max[res$id], 60)
  is_n <- res$metric %in% c("mean_patients", "mean_n")
  expect_true(all(res$estimate[is_n] <= n_cap[is_n]))
  # The first cohort on the contour is 1 to the number of cohorts: 10, 20,
  # 40 or, with cohorts of 1, 60.
  max_cohorts <- rep(20, nrow(res))
  max_cohorts[res$id %in% c("N_30", "N_30_iso")] <- 10
  max_cohorts[res$id %in% c("N_120", "N_120_iso")] <- 40
  max_cohorts[res$id %in% c("cohort_1", "cohort_1_iso")] <- 60
  fc <- res$metric == "first_contour_cohort"
  expect_true(all(res$estimate[fc] >= 1 & res$estimate[fc] <= max_cohorts[fc]))
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

  # Each calibration file has one selected epsilon, whose false discovery
  # rate is within its alpha. Alpha grows with the three alpha files: 0.05,
  # 0.10, 0.20, and so does the selected epsilon.
  cals <- list(baseline = cal, alpha_0.10 = read("calibration_alpha_0.10.csv"),
               alpha_0.20 = read("calibration_alpha_0.20.csv"))
  for (nm in names(cals)) {
    cl <- cals[[nm]]
    expect_equal(sum(cl$selected), 1, info = nm)
    expect_lte(cl$fdr[cl$selected], cl$alpha[1])
  }
  expect_equal(vapply(cals, function(x) x$alpha[1], 0, USE.NAMES = FALSE), c(0.05, 0.10, 0.20))
  sel <- vapply(cals, function(x) x$eps[x$selected], 0, USE.NAMES = FALSE)
  expect_true(all(diff(sel) > 0))

  # The shift calibrations: alpha 0.05, and the band and the shift of the
  # file name (the shift is the band). The older files have no shift column
  # and shift 0, with band 0.1. The selected epsilon is a tried value with
  # an FDR within alpha, and it is one value in the file.
  scals <- lapply(setNames(shifts, shifts), function(b) read(sprintf("calibration_shift_%s.csv", b)))
  for (b in shifts) {
    cl <- scals[[b]]
    expect_equal(unique(cl$band), as.numeric(b), info = b)
    expect_equal(unique(cl$shift), as.numeric(b), info = b)
    expect_equal(unique(cl$alpha), 0.05, info = b)
    expect_equal(sum(cl$selected), 1, info = b)
    expect_lte(cl$fdr[cl$selected], 0.05)
    expect_equal(unique(cl$n_scenarios), 56, info = b)
  }
  for (cl in cals) {
    expect_equal(unique(cl$band), 0.1)
    expect_false("shift" %in% names(cl))
  }
  sel_shift <- vapply(scals, function(x) x$eps[x$selected], 0)

  # Rows that use a calibrated epsilon used the one selected in their file:
  # the baseline and every row that keeps its epsilon use the baseline
  # calibration, the two alpha rows use their own, and each shift row uses
  # the calibration of its shift. Fixed-epsilon rows have the value in their
  # id, the M&S row has none, and the fixed-epsilon shift row keeps the
  # baseline epsilon.
  eps_of <- function(id) runs$eps[runs$id == id]
  fixed <- c(eps_0.5 = 0.5, eps_0.4 = 0.4, eps_0.3 = 0.3, eps_0.15 = 0.15, eps_0.05 = 0.05,
             eps_0.4_pick_all = 0.4, eps_0.3_pick_all = 0.3, eps_0.3_startup = 0.3,
             eps_0.4_startup = 0.4, eps_0.3_pick_all_startup = 0.3, eps_0.4_pick_all_startup = 0.4,
             eps_0.3_pick_all_iso = 0.3, eps_0.4_pick_all_iso = 0.4, eps_0.5_iso = 0.5,
             eps_0.4_pick_all_startup_iso = 0.4)
  for (id in names(fixed)) expect_equal(eps_of(id), fixed[[id]], info = id)
  expect_true(is.na(eps_of("ref_ms")))
  expect_true(is.na(eps_of("ref_ms_iso")))
  expect_equal(eps_of("alpha_0.10"), sel[2])
  expect_equal(eps_of("alpha_0.20"), sel[3])
  shift_of <- c(shift_0.1 = "0.1", shift_0.075 = "0.075", shift_0.05 = "0.05", shift_0.025 = "0.025",
                shift_0.1_pick_all = "0.1", shift_0.1_pick_all_fewest = "0.1",
                shift_0.1_pick_all_interval = "0.1", shift_0.05_pick_all = "0.05",
                shift_0.1_pick_all_fewest_startup = "0.1")
  for (id in names(shift_of)) expect_equal(eps_of(id), sel_shift[[shift_of[[id]]]], info = id)
  own <- c(names(fixed), "ref_ms", "ref_ms_iso", "alpha_0.10", "alpha_0.20", names(shift_of))
  for (id in setdiff(ids, own)) expect_equal(eps_of(id), sel[1], info = id)

  # Mechanism of the shift: with a shift of 0.1 and the same epsilon, a dose
  # with truth 0.30 is called safe more often, so the true MTDs are in the
  # final set more often than at the baseline. (With the smaller calibrated
  # epsilon the two effects cancel, so only the fixed-epsilon row is checked.)
  p_sel <- function(id) {
    r <- res[res$id == id & res$metric == "pct_selected" & res$combo %in% mtd, ]
    sum(r$estimate)
  }
  expect_gt(p_sel("shift_0.1_eps_fixed"), p_sel("baseline"))

  # Start-up phase (#51): the mechanism file has one row for each start-up
  # configuration. The phase has 1 to 5 cohorts (the diagonal of the 5 x 3
  # grid), and nearly every trial's phase ends by a toxicity. The partner of
  # each start-up row has the same epsilon and exists in the results.
  ph <- read("mtd_startup_phase.csv")
  expect_setequal(ph$id, c(with_startup, grep("startup_iso$", with_iso, value = TRUE)))
  expect_false(anyNA(ph))
  expect_true(all(ph$mean_n_startup >= 1 & ph$mean_n_startup <= 5))
  expect_true(all(ph$pct_phase_ended_by_toxicity > 90 & ph$pct_phase_ended_by_toxicity <= 100))
  expect_true(all(ph$se_n_startup > 0))
  partner <- sub("_startup$", "", with_startup)
  expect_true(all(partner %in% ids))
  expect_equal(runs$eps[match(with_startup, runs$id)], runs$eps[match(partner, runs$id)])
  # The phase moves patients up the diagonal early: the first cohort on the
  # contour comes sooner than in the partner row.
  fcc <- function(id) {
    vapply(id, function(i) res$estimate[res$id == i & res$metric == "first_contour_cohort"], 0)
  }
  expect_true(all(fcc(with_startup) < fcc(partner)))

  # Isotonic pick (#54): the pick classes are stored for the isotonic rows
  # only. They agree with the stored pct_pick of the row (PCS, over, under),
  # the classes cannot add up to more than 100, and the 15 combinations add up
  # to 100 minus the share without a pick. The partner of each row is the id
  # without "_iso" and has the same epsilon.
  pc <- read("mtd_pick_classes.csv")
  expect_setequal(pc$id, with_iso)
  expect_false(anyNA(pc))
  combo_cols <- paste0("pct_pick_", labels)
  expect_true(all(combo_cols %in% names(pc)))
  expect_equal(ncol(pc), 1 + 4 + 15)
  truth <- c(0.05, 0.10, 0.15, 0.30, 0.45, 0.10, 0.15, 0.30, 0.45, 0.55,
             0.15, 0.30, 0.45, 0.50, 0.60)
  for (id in with_iso) {
    p <- pc[pc$id == id, ]
    pick <- res[res$id == id & res$metric == "pct_pick", ]
    expect_equal(unlist(p[combo_cols]), setNames(pick$estimate, combo_cols), tolerance = 1e-5, info = id)
    expect_equal(p$pct_pick_mtd, 100 * res$estimate[res$id == id & res$metric == "pcs"], tolerance = 1e-5,
                 info = id)
    expect_equal(p$pct_pick_over, sum(pick$estimate[truth >= 0.45 - 1e-8]), tolerance = 1e-5, info = id)
    expect_equal(p$pct_pick_under, sum(pick$estimate[truth <= 0.15 + 1e-8]), tolerance = 1e-5, info = id)
    expect_equal(sum(unlist(p[combo_cols])) + p$pct_no_pick, 100, tolerance = 1e-4, info = id)
    expect_lte(p$pct_pick_mtd + p$pct_pick_over + p$pct_pick_under, 100 + 1e-6)
  }
  iso_partner <- sub("_iso$", "", with_iso)
  expect_true(all(iso_partner %in% ids))
  expect_equal(runs$eps[match(with_iso, runs$id)], runs$eps[match(iso_partner, runs$id)])
  # The pick changes the single dose only: the trial conduct and the final set
  # are the partner's, so FDR and TPR agree with the partner's within Monte
  # Carlo error (differences below 1.5 points here).
  for (m in c("fdr", "tpr")) {
    gap <- vapply(seq_along(with_iso), function(i) {
      est <- function(id) res$estimate[res$id == id & res$metric == m]
      abs(est(with_iso[i]) - est(iso_partner[i]))
    }, 0)
    expect_true(all(gap < 0.015), info = m)
  }

  # Riviere's published PCS for scenario 1 (percent): seven designs.
  expect_equal(nrow(ref), 7)
  expect_true(all(ref$pcs > 60 & ref$pcs < 80))
})
