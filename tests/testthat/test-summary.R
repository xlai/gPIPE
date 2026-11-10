g22 <- DoseGrid(levels = c(A = 2, B = 2))   # combinations 1..4

# A hand-built simulation: 4 trials, 4 combinations, truth known.
hand_sim <- function(target = "MTD") {
  des <- if (target == "MTD") {
    ms_design(g22, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 4, N = 12)
  } else {
    cheung_design(g22, theta = 0.5, eps = 0.1, N = 12)
  }
  fs <- rbind(c(1, 1, 0, 0), c(1, 0, 0, 0), c(0, 1, 1, 0), c(0, 0, 0, 0)) == 1
  PipeSimulation(
    design = des,
    summary = data.frame(trial = 1:4, n_enrolled = c(12, 12, 12, 6),
                         stopped = c(FALSE, FALSE, FALSE, TRUE), pick = c(2L, 1L, 3L, NA)),
    n = rbind(c(6, 6, 0, 0), c(4, 4, 4, 0), c(0, 6, 3, 3), c(3, 3, 0, 0)),
    s1 = matrix(0, 4, 4),
    gamma = rbind(c(0, 0, 1, 1), c(0, 1, 1, 1), c(0, 0, 0, 1), c(1, 1, 1, 1)),
    final_set = fs, trials = list(),
    truth = c(0.1, 0.3, 0.45, 0.6)
  )
}

test_that("summary() matches hand-computed tables for an MTD design", {
  s <- summary(hand_sim(), band = 0.075)  # safe <= 0.225, toxic >= 0.375
  expect_s3_class(s, "gPIPE::PipeSummary")
  expect_equal(as.character(s@by_combo$class), c("favourable", "indifferent", "unfavourable", "unfavourable"))
  expect_equal(s@by_combo$pct_selected, c(50, 50, 25, 0))
  expect_equal(s@by_combo$mean_n, c(3.25, 4.75, 1.75, 0.75))
  expect_equal(s@by_combo$pct_patients, 100 * c(13, 19, 7, 3) / 42)
  expect_equal(s@by_combo$pct_trials_tried, c(75, 100, 50, 25))
  expect_equal(s@trials$mean_n, 10.5)
  expect_equal(s@trials$pct_stopped, 25)
  m <- setNames(s@metrics$estimate, s@metrics$metric)
  # Declared safe (gamma = 0): {1,2}, {1}, {1,2,4}->{1,2}... trial 3 is c(0,0,0,1): {1,2,3}; trial 4 stopped.
  # Trial 1: declared {1,2}: 0 unfavourable -> FDR 0, TPR 1.
  # Trial 2: declared {1}: FDR 0, TPR 1.
  # Trial 3: declared {1,2,3}: one unfavourable (3) of 3 -> FDR 1/3, TPR 1.
  # Trial 4: stopped -> nothing declared: FDR 0, TPR 0.
  expect_equal(m[["fdr"]], (0 + 0 + 1 / 3 + 0) / 4)
  expect_equal(m[["tpr"]], (1 + 1 + 1 + 0) / 4)
  # Only combination 2 is exactly at theta, picks 2, 1, 3, NA.
  expect_equal(m[["pcs"]], 0.25)
  expect_equal(m[["mean_recommended"]], (2 + 1 + 2 + 0) / 4)
  # Share of patients at unfavourable combinations 3, 4.
  expect_equal(m[["pct_patients_unfavourable"]], mean(100 * c(0, 4, 6, 0) / c(12, 12, 12, 6)))
  expect_equal(m[["p_set_unfavourable"]], 0.25)
  expect_true(all(is.finite(s@metrics$mc_se)))
})

test_that("the Mander & Sweeting recommendation table pools final sets by band", {
  s <- summary(hand_sim())
  # d = (-20, 0, 15, 30) points: more than 10 below, at theta, more than 10 above x2.
  expect_equal(s@bands$band, c("more than 10 below", "1-10 below", "at theta", "1-10 above",
                               "more than 10 above"))
  expect_equal(s@bands$recommended, c(2, 0, 2, 0, 1))
  expect_equal(s@bands$pct_recommended, 100 * c(2, 0, 2, 0, 1) / 5)
  s2 <- summary(hand_sim(), breaks = c(0, 0.25, 0.4, 1))
  expect_equal(s2@bands$recommended, c(2, 2, 1))
})

test_that("MED designs declare gamma = 1 and swap the favourable side", {
  s <- summary(hand_sim("MED"))  # theta 0.5, band 0: only combination 4 (0.6) is effective
  expect_equal(as.character(s@by_combo$class), c("unfavourable", "unfavourable", "unfavourable", "favourable"))
  expect_false("pcs" %in% s@metrics$metric)
  m <- setNames(s@metrics$estimate, s@metrics$metric)
  # Declared (gamma = 1): {3,4}, {2,3,4}, {4}, stopped -> none.
  expect_equal(m[["fdr"]], (1 / 2 + 2 / 3 + 0 + 0) / 4)
  expect_equal(m[["tpr"]], (1 + 1 + 1 + 0) / 4)
})

test_that("summary() needs the truth and a valid band", {
  sim <- hand_sim()
  sim@truth <- numeric()
  expect_error(summary(sim), "truth")
  expect_error(summary(hand_sim(), band = -1), "band")
})

test_that("simulate() stores the truth, from a vector, an array or a normal scenario", {
  des <- ms_design(g22, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 4, N = 6)
  tr <- c(0.1, 0.2, 0.3, 0.5)
  sim <- simulate(des, nsim = 3, seed = 1, truth = tr)
  expect_equal(sim@truth, tr)
  s <- summary(sim)
  expect_equal(nrow(s@by_combo), 4L)
  expect_output(print(s), "Operating characteristics")
  ch <- cheung_design(g22, theta = 0.5, eps = 0.1, N = 8)
  sim2 <- simulate(ch, nsim = 3, seed = 1, truth = list(mean = c(0, 0.5, 1, 1), sd = 1))
  expect_equal(sim2@truth, c(0, 0.5, 1, 1))
  expect_s3_class(summary(sim2), "gPIPE::PipeSummary")
})
