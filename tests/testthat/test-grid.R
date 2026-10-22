test_that("full 2D grid numbers combinations in array order", {
  g <- DoseGrid(levels = c(A = 3, B = 2))
  expect_equal(n_combos(g), 6L)
  expect_equal(unname(g@combos[, "A"]), c(1L, 2L, 3L, 1L, 2L, 3L))
  expect_equal(unname(g@combos[, "B"]), c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_equal(combo_labels(g)[1:2], c("A1.B1", "A2.B1"))

  m <- matrix(1:6, 3, 2)
  expect_equal(grid_from_array(g, m), 1:6)
  expect_equal(grid_as_array(g, 1:6), m)
})

test_that("combo_index maps levels to indices and rejects off-grid levels", {
  g <- DoseGrid(levels = c(4, 3))
  expect_equal(combo_index(g, c(2, 3)), 10L)
  expect_equal(combo_index(g, rbind(c(1, 1), c(5, 1), c(4, 3))), c(1L, NA, 12L))
})

test_that("cover edges on a full grid are the unit steps", {
  g <- DoseGrid(levels = c(3, 3))
  # 2 * 3 * (3 - 1) unit steps
  expect_equal(nrow(g@edges), 12L)
  steps <- g@combos[g@edges[, "upper"], ] - g@combos[g@edges[, "lower"], ]
  expect_true(all(rowSums(steps) == 1L))
  expect_true(all(steps >= 0L))
})

test_that("cover edges on a subset skip the missing combination", {
  keep <- matrix(TRUE, 3, 3)
  keep[2, 2] <- FALSE
  g <- DoseGrid(levels = c(3, 3), keep = keep)
  expect_equal(n_combos(g), 8L)
  lbl <- combo_labels(g)
  e <- paste(lbl[g@edges[, "lower"]], lbl[g@edges[, "upper"]])
  # with (2,2) gone, (1,2) is covered by (3,2) and (2,1) by (2,3)
  expect_true("A1.B2 A3.B2" %in% e)
  expect_true("A2.B1 A2.B3" %in% e)
  expect_false("A1.B2 A2.B3" %in% e)
  expect_false("A1.B1 A2.B3" %in% e)
  expect_true(is.na(combo_index(g, c(2, 2))))
})

test_that("subset cover edges match the full-grid rule when nothing is dropped", {
  g_full <- DoseGrid(levels = c(3, 2, 2))
  g_sub <- DoseGrid(levels = c(3, 2, 2), keep = rep(TRUE, 12))
  sub_edges <- gPIPE:::grid_cover_edges(g_sub@combos, g_sub@lookup, FALSE)
  expect_equal(sub_edges, g_full@edges)
})

test_that("up-set check and extreme elements follow the partial order", {
  g <- DoseGrid(levels = c(2, 2))
  expect_true(gPIPE:::is_upset(g, c(0, 0, 1, 1)))
  expect_false(gPIPE:::is_upset(g, c(1, 0, 0, 0)))
  inside <- c(FALSE, TRUE, TRUE, TRUE)
  expect_equal(sort(gPIPE:::minimal_elements(g, inside)), c(2L, 3L))
  expect_equal(gPIPE:::maximal_elements(g, inside), 4L)
})

test_that("neighbourhood is a Chebyshev ball of radius 1", {
  g <- DoseGrid(levels = c(4, 4))
  k <- combo_index(g, c(1, 1))
  expect_equal(sort(gPIPE:::grid_neighbourhood(g, k)), combo_index(g, rbind(c(1, 1), c(2, 1), c(1, 2), c(2, 2))))
  expect_length(gPIPE:::grid_neighbourhood(g, combo_index(g, c(2, 2))), 9L)
})

test_that("invalid grids are rejected", {
  expect_error(DoseGrid(levels = c(0, 2)))
  expect_error(DoseGrid(levels = c(2, 2), keep = rep(FALSE, 4)))
  expect_error(DoseGrid(levels = c(2, 2), keep = TRUE))
  expect_error(DoseGrid(levels = c(2.5, 2)), "whole numbers")
  expect_error(DoseGrid(levels = c(NA, 2)), "whole numbers")
})

test_that("the validator catches inconsistent properties", {
  g <- DoseGrid(levels = c(3, 2))
  expect_error(g@levels <- c(A = 2L, B = 2L), "dimensions of the full grid")
  expect_error(g@combos <- g@combos[, 1, drop = FALSE], "one column per drug")
  expect_error(g@edges <- rbind(g@edges, c(1L, 99L)), "pairs of combination indices")
})

test_that("a three-drug grid round-trips between array and vector views", {
  g <- DoseGrid(levels = c(A = 2, B = 3, C = 2))
  a <- array(seq_len(12) * 10, dim = c(2, 3, 2))
  expect_equal(grid_from_array(g, a), as.vector(a))
  expect_equal(grid_as_array(g, grid_from_array(g, a)), a)
  expect_equal(combo_index(g, c(2, 3, 2)), 12L)
  expect_equal(nrow(g@edges), 1 * 3 * 2 + 2 * 2 * 2 + 2 * 3 * 1)
})
