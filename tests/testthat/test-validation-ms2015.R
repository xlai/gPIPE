# Validation against Mander & Sweeting (2015), Tables III and IV (#31, step 6).
#
# The results are stored in inst/simulations/ms2015/results: gPIPE and
# pipe.design 0.5.1 at 4000 (Table III) or 8000 (Table IV) trials per
# configuration, with the Monte Carlo covariance of each configuration's
# cells, and the printed tables. The scripts beside them make the runs again.
#
# The cells are the experimentation and recommendation percentages by
# true-toxicity band and the mean number of recommended doses. The criterion
# is one chi-square test over all configurations (see compare.R), not a limit
# on each cell.
#
# This checks the stored runs, so it does not change when the package code
# changes; the golden tests check the code. Make the gPIPE run again after a
# change to the M&S rules.

ms2015_dir <- function() system.file("simulations", "ms2015", package = "gPIPE")

test_that("gPIPE reproduces pipe.design on M&S Tables III and IV", {
  skip_on_cran()
  skip_if(ms2015_dir() == "", "simulation results are not installed")
  source(file.path(ms2015_dir(), "compare.R"), local = TRUE)
  cmp <- ms2015_compare(file.path(ms2015_dir(), "results"))

  # Every configuration and cell is there: 16 x 11 (III) + 7 x 9 (IV); the
  # paper has no mean number recommended for Table IV.
  expect_equal(nrow(cmp$pipe_design), 16 * 11 + 7 * 9)
  expect_equal(nrow(cmp$paper), 16 * 11 + 7 * 8)
  expect_equal(nrow(cmp$by_config), 23)

  # Against pipe.design: both runs have a covariance, so the test is exact up
  # to the normal approximation.
  expect_gt(cmp$test_pipe_design[["p"]], 0.05)
  # Against the printed tables the covariance is borrowed from the gPIPE run,
  # so the test is approximate (pipe.design itself gets p = 0.06 against the
  # paper this way) and the level is 0.01.
  expect_gt(cmp$test_paper[["p"]], 0.01)
})
