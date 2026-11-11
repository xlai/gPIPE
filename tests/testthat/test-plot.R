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
  # positions: x is the drug B level, y the drug A level, row by row
  d <- b$data[[1]]
  expect_equal(d$x, match(g@combos[, 2], sort(unique(g@combos[, 2]))))
  expect_equal(d$y, match(g@combos[, 1], sort(unique(g@combos[, 1]))))
  expect_equal(b$data[[2]]$label, sprintf("%.1f", 100 * colMeans(sim@final_set)))
  expect_equal(b$plot$data$value, 100 * colMeans(sim@final_set))
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
  pts <- b$data[[2]]
  n <- nrow(cal@trace)
  col <- unique(pts$colour)
  expect_length(col, 2L)
  fdr_col <- pts$colour[1]  # FDR rows come first
  expect_equal(pts$y[pts$colour == fdr_col], cal@trace$fdr)
  expect_equal(pts$y[pts$colour != fdr_col], cal@trace$tpr)
  expect_equal(pts$x[1:n], cal@trace$eps)
})

test_that("plot() on a calibration with no feasible eps has no vline", {
  g <- DoseGrid(levels = c(A = 2, B = 2))
  des <- cheung_design(g, theta = 0.3, eps = 0.5, target = "MTD", outcome = "binary",
                       prior_a = 0.3, prior_b = 0.7, N = 8)
  cal <- suppressWarnings(calibrate_eps(des, n_scenarios = 2, nsim = 5, tol = 0.2,
                                        alpha = 1e-6, lower = 0.3, seed = 1))
  expect_true(is.na(cal@eps))
  p <- plot(cal)
  expect_s3_class(p, "ggplot")
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomVline"), logical(1))))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot() on a calibration with one eps or non-finite values is quiet", {
  g <- DoseGrid(levels = c(A = 2, B = 2))
  des <- cheung_design(g, theta = 0.3, eps = 0.5, target = "MTD", outcome = "binary",
                       prior_a = 0.3, prior_b = 0.7, N = 8)
  cal <- calibrate_eps(des, n_scenarios = 3, nsim = 10, tol = 0.2, alpha = 0.99, seed = 1)
  expect_no_message(ggplot2::ggplotGrob(plot(cal)))
  expect_no_warning(ggplot2::ggplotGrob(plot(cal)))
  cal@trace$tpr <- NaN
  expect_no_warning(ggplot2::ggplotGrob(plot(cal)))
  one <- cal
  one@trace <- cal@trace[1, ]
  expect_no_message(ggplot2::ggplotGrob(plot(one)))
})

test_that("plot() on a one-drug grid draws bars in dose order", {
  g <- DoseGrid(levels = c(A = 4))
  des <- cheung_design(g, theta = 0.3, eps = 0.2, target = "MTD", outcome = "binary", N = 6)
  sim <- simulate(des, nsim = 3, seed = 1, truth = seq(0.1, 0.5, length.out = 4))
  b <- ggplot2::ggplot_build(plot(sim))
  expect_equal(b$data[[1]]$x, 1:4)
  expect_equal(b$data[[1]]$y, 100 * colMeans(sim@final_set))
  expect_equal(levels(b$plot$data$combo), combo_labels(g))
})

test_that("plot() errors clearly when ggplot2 is missing", {
  local_mocked_bindings(requireNamespace = function(...) FALSE, .package = "base")
  g <- DoseGrid(levels = c(A = 2))
  des <- cheung_design(g, theta = 0.3, eps = 0.2, target = "MTD", outcome = "binary", N = 4)
  sim <- simulate(des, nsim = 2, seed = 1, truth = c(0.1, 0.3))
  expect_error(plot(sim), "needs the ggplot2 package")
})
