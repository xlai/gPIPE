skip_if_not_installed("ggplot2")

test_that("plot() on a simulation draws a tile map for two drugs and bars otherwise", {
  g <- DoseGrid(levels = c(A = 2, B = 3))
  des <- ms_design(g, theta = 0.3, prior_median = 0.3, prior_strength = 1 / 6, N = 9)
  sim <- simulate(des, nsim = 4, seed = 1, truth = seq(0.1, 0.5, length.out = 6))
  p <- plot(sim)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(nrow(b$data[[1]]), 6L)
  expect_equal(sort(b$data[[2]]$label), sort(sprintf("%.1f", 100 * colMeans(sim@final_set))))
  p2 <- plot(sim, type = "patients")
  expect_equal(sort(ggplot2::ggplot_build(p2)$data[[2]]$label), sort(sprintf("%.1f", colMeans(sim@n))))
  expect_error(plot(sim, type = "x"))

  g3 <- DoseGrid(levels = c(A = 2, B = 2, C = 2))
  des3 <- cheung_design(g3, theta = 0.3, eps = 0.2, target = "MTD", outcome = "binary", N = 6)
  sim3 <- simulate(des3, nsim = 3, seed = 1, truth = seq(0.1, 0.5, length.out = 8))
  b3 <- ggplot2::ggplot_build(plot(sim3))
  expect_equal(b3$data[[1]]$y, 100 * colMeans(sim3@final_set))
})

test_that("plot() on a calibration shows FDR and TPR with alpha and the chosen eps", {
  g <- DoseGrid(levels = c(A = 2, B = 2))
  des <- cheung_design(g, theta = 0.3, eps = 0.5, target = "MTD", outcome = "binary",
                       prior_a = 0.3, prior_b = 0.7, N = 8)
  cal <- calibrate_eps(des, n_scenarios = 3, nsim = 10, tol = 0.2, seed = 1)
  p <- plot(cal)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(nrow(b$data[[1]]), 2L * nrow(cal@trace))
  expect_equal(b$data[[3]]$yintercept, cal@alpha)
  expect_equal(b$data[[4]]$xintercept, cal@eps)
})
