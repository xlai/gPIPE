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
  s <- summary(hand_sim(), band = 0.075)  # below <= 0.225, above >= 0.375
  expect_s3_class(s, "gPIPE::PipeSummary")
  expect_equal(as.character(s@by_combo$class), c("below", "indifferent", "above", "above"))
  expect_equal(s@by_combo$pct_selected, c(50, 50, 25, 0))
  expect_equal(s@by_combo$mean_n, c(3.25, 4.75, 1.75, 0.75))
  expect_equal(s@by_combo$pct_patients, 100 * c(13, 19, 7, 3) / 42)
  expect_equal(s@by_combo$pct_trials_tried, c(75, 100, 50, 25))
  expect_equal(s@trials$mean_n, 10.5)
  expect_equal(s@trials$pct_stopped, 25)
  m <- setNames(s@metrics$estimate, s@metrics$metric)
  # FDR is over the final set. Final sets: {1,2}, {1}, {2,3}, {} (stopped); combination 3 is
  # above (true gamma = 1, wrong for MTD), so trial 3 gives 1/2 and the rest 0.
  # TPR is over the combinations with gamma = 0 in the final estimate: {1,2}, {1}, {1,2,3}, none
  # (stopped); combination 1 is the only one with true gamma = 0.
  expect_equal(m[["fdr"]], (0 + 0 + 1 / 2 + 0) / 4)
  expect_equal(m[["tpr"]], (1 + 1 + 1 + 0) / 4)
  # Only combination 2 is exactly at theta, picks 2, 1, 3, NA.
  expect_equal(m[["pcs"]], 0.25)
  expect_equal(m[["mean_recommended"]], (2 + 1 + 2 + 0) / 4)
  # Share of patients at the above combinations 3, 4, below combination 1, and in band combination 2.
  expect_equal(m[["pct_patients_above"]], mean(100 * c(0, 4, 6, 0) / c(12, 12, 12, 6)))
  expect_equal(m[["pct_patients_below"]], mean(100 * c(6, 4, 0, 3) / c(12, 12, 12, 6)))
  expect_equal(m[["pct_patients_in_band"]], mean(100 * c(6, 4, 6, 3) / c(12, 12, 12, 6)))
  expect_equal(m[["p_set_false"]], 0.25)
  expect_true(all(is.finite(s@metrics$mc_se)))
})

test_that("the Mander & Sweeting recommendation table pools final sets by band", {
  s <- summary(hand_sim())
  # d = (-20, 0, 15, 30) points: more than 10 below, at theta, more than 10 above x2.
  expect_equal(s@bands$band, c("more than 10 below", "within 10 below", "at theta", "within 10 above",
                               "more than 10 above", "None"))
  # Trial 4 recommended nothing and counts once as None.
  expect_equal(s@bands$recommended, c(2, 0, 2, 0, 1, 1))
  expect_equal(s@bands$pct_recommended, 100 * c(2, 0, 2, 0, 1, 1) / 6)
  s2 <- summary(hand_sim(), breaks = c(0, 0.25, 0.4, 1))
  expect_equal(s2@bands$recommended, c(2, 2, 1, 1))
})

test_that("MED designs look for gamma = 1, MTD for gamma = 0", {
  s <- summary(hand_sim("MED"))  # theta 0.5, normal outcome so band 0: only combination 4 (0.6) is above
  expect_equal(as.character(s@by_combo$class), c("below", "below", "below", "above"))
  expect_false("pcs" %in% s@metrics$metric)
  m <- setNames(s@metrics$estimate, s@metrics$metric)
  # FDR is over the final sets {1,2}, {1}, {2,3}, {}: all below (true gamma = 0, wrong for MED), so 1, 1, 1, 0.
  expect_equal(m[["fdr"]], 0.75)
  # TPR is over the declared sets (gamma = 1): {3,4}, {2,3,4}, {4}, stopped.
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

test_that("band edges are robust to floating point", {
  # Exactly 10 points from theta is in "within 10", not "more than 10".
  b <- ms_band_table(c(1, 1, 1, 1), c(0.2, 0.4, 0.5, 0.3), 0.3, NULL, 0)
  expect_equal(b$recommended[1:5], c(0, 1, 1, 1, 1))
  b2 <- ms_band_table(c(1, 1), c(0.3 + 0.05, 0.1 + 0.2), 0.3, c(0, 0.35, 1), 0)
  expect_equal(b2$recommended, c(1, 1, 0))
  cls <- classify_truth(0.3, c(0.1 + 0.2, 0.5), 0)
  expect_equal(as.character(cls), c("below", "above"))
  # A value exactly at the band edge is outside the interval.
  expect_equal(as.character(classify_truth(0.3, c(0.3 - 0.075, 0.3 + 0.075, 0.3), 0.075)),
               c("below", "above", "indifferent"))
})

test_that("summary() checks `breaks` and `band`", {
  sim <- hand_sim()  # truth 0.1 to 0.6
  expect_error(summary(sim, breaks = c(0, 0.25, 0.35)), "cover")
  expect_error(summary(sim, breaks = c(0.5, 0.1)), "increasing")
  expect_error(summary(sim, breaks = 0.5), "two")
  expect_error(summary(sim, band = NA_real_), "band")
  expect_error(summary(sim, band = "a"), "band")
})

test_that("one band sets the classes and the recommendation table", {
  sim <- hand_sim()  # theta 0.3, truth 0.1, 0.3, 0.45, 0.6
  # The default is 0.1 for a binary outcome: 0.1 is below, 0.3 in band, 0.45 and 0.6 above.
  s <- summary(sim)
  expect_equal(s@band, 0.1)
  expect_equal(as.character(s@by_combo$class), c("below", "indifferent", "above", "above"))
  expect_equal(summary(hand_sim("MED"))@band, 0)  # normal outcome
  # The labels follow the width.
  s5 <- summary(sim, band = 0.05)
  expect_equal(s5@bands$band[c(1, 2, 4, 5)],
               c("more than 5 below", "within 5 below", "within 5 above", "more than 5 above"))
  # 0.1 is 0.2 below theta, 0.45 is 0.15 above, 0.6 is 0.3 above.
  expect_equal(s5@bands$recommended[c(1, 3, 5)], c(2, 2, 1))
  expect_equal(summary(sim, band = 0.2)@bands$band[1], "more than 20 below")
  # mtd_tol only changes pcs: truth 0.3 is a true MTD only within mtd_tol of theta.
  pcs <- function(x) setNames(x@metrics$estimate, x@metrics$metric)[["pcs"]]
  expect_equal(pcs(summary(sim)), 0.25)
  # Picks are 2, 1, 3 and none; within 0.2 of theta are combinations 1, 2 and 3.
  expect_equal(pcs(summary(sim, mtd_tol = 0.2)), 0.75)
  expect_equal(pcs(summary(sim, band = 0.2)), 0.25)
})

test_that("FDR is 0 when no combination has the wrong true gamma", {
  sim <- hand_sim()
  sim@truth <- rep(0.1, 4)  # all below, none above (wrong for MTD)
  s <- summary(sim)
  m <- setNames(s@metrics$estimate, s@metrics$metric)
  expect_equal(m[["fdr"]], 0)
  expect_true(is.na(m[["p_set_false"]]))
})
