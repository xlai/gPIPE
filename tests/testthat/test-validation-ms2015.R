# Validation against Mander & Sweeting (2015), Tables III and IV (#31, step 6).
#
# The results are stored in inst/simulations/ms2015/results: gPIPE and
# pipe.design 0.5.1 at 4000 (Table III) or 8000 (Table IV) trials per
# configuration, and the printed tables. The scripts beside them make the
# runs again. Here the stored results are compared cell by cell: experimentation
# and recommendation percentages by true-toxicity band, and the mean number
# of recommended doses. z is the difference in Monte Carlo standard errors.


# Pass margin, provisional (#31). With 239 cells, about one cell in two
# runs passes 3 SE by chance alone, so the limit is 4 SE on every cell.
MS2015_MAX_Z <- 4

ms2015_dir <- function() system.file("simulations", "ms2015", package = "gPIPE")

test_that("gPIPE reproduces pipe.design on M&S Tables III and IV", {
  skip_on_cran()
  skip_if(ms2015_dir() == "", "simulation results are not installed")
  source(file.path(ms2015_dir(), "compare.R"), local = TRUE)
  cmp <- ms2015_compare(file.path(ms2015_dir(), "results"))

  # Every configuration and cell is there: 16 x 11 (III) + 7 x 9 (IV).
  expect_equal(nrow(cmp$pipe_design), 16 * 11 + 7 * 9)
  expect_equal(nrow(cmp$paper), 16 * 11 + 7 * 8)

  expect_lt(max(abs(cmp$pipe_design$z)), MS2015_MAX_Z)
  expect_lt(max(abs(cmp$paper$z)), MS2015_MAX_Z)
})
